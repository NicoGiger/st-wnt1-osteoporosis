#!/usr/bin/env python3
"""Reproducible SpatialDM pipeline for Visium Wnt sections.

This script converts the notebook workflow into a command-line pipeline.
It supports two ligand-receptor modes:
1) cellchat: canonical mouse LR pairs via SpatialDM/CellChatDB
2) wnt1_all: custom Wnt1-versus-all-genes screen

Example
-------
python spatialdm_pipeline.py \
  --visium-dir /path/to/Wnt \
  --counts-file raw_feature_bc_matrix.h5 \
  --group-csv /path/to/spa.csv \
  --group-column group \
  --section Wt \
  --lr-mode wnt1_all \
  --output-dir results/spatialdm_wildtype
"""

from __future__ import annotations

import argparse
import json
import logging
import os
import sys
import time
from pathlib import Path
from typing import Optional

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import requests
import scanpy as sc
import squidpy as sq
import spatialdm as sdm
import spatialdm.plottings as sdm_pl


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Run a reproducible SpatialDM pipeline on a Visium section.")
    p.add_argument("--visium-dir", required=True, help="Directory passed to squidpy.read.visium().")
    p.add_argument("--counts-file", default="raw_feature_bc_matrix.h5", help="Counts file inside visium-dir.")
    p.add_argument("--group-csv", default=None, help="Optional CSV with per-spot section labels/clusters.")
    p.add_argument("--group-column", default="group", help="Column in group CSV containing the section label.")
    p.add_argument("--section", default=None, help="Subset to one section/group, e.g. 'Wt' or 'WntOverExp'.")
    p.add_argument("--index-col", default=0, type=int, help="Index column for the group CSV.")
    p.add_argument("--filter-in-tissue", action="store_true", help="Keep only spots with obs['in_tissue'] == 1.")
    p.add_argument("--min-gene-counts", type=int, default=3, help="Filter genes with total counts below this threshold.")
    p.add_argument("--target-sum", type=float, default=1e4, help="Library-size normalization target.")
    p.add_argument("--skip-normalization", action="store_true", help="Assume adata.X is already normalized/log-transformed.")
    p.add_argument("--lr-mode", choices=["cellchat", "wnt1_all"], default="cellchat", help="Ligand-receptor set.")
    p.add_argument("--species", default="mouse", help="Species for ligand-receptor database selection.")
    p.add_argument("--lr-db-csv", default=None, help="Optional local ligand-receptor CSV with ligand/receptor columns. Recommended for reproducibility.")
    p.add_argument("--lr-db-url", default=None, help="Optional URL to a ligand-receptor CSV. Used only if --lr-db-csv is not provided.")
    p.add_argument("--ligand-gene", default="Wnt1", help="Ligand gene for --lr-mode wnt1_all.")
    p.add_argument("--min-cell", type=int, default=3, help="Min expressing spots for CellChatDB extraction.")
    p.add_argument("--l", type=float, default=100.0, help="SpatialDM RBF kernel length scale.")
    p.add_argument("--n-neighbors", type=int, default=None, help="Optional n_neighbors for weight_matrix.")
    p.add_argument("--single-cell", action="store_true", help="Pass single_cell=True to SpatialDM.")
    p.add_argument("--method", default="z-score", help="SpatialDM inference method.")
    p.add_argument("--n-perm", type=int, default=1000, help="Number of permutations for global/local tests.")
    p.add_argument("--nproc", type=int, default=1, help="Number of CPU processes.")
    p.add_argument("--pair-fdr", type=float, default=0.1, help="FDR threshold for significant LR pairs.")
    p.add_argument("--spot-p", type=float, default=0.1, help="P-value threshold for significant local spots.")
    p.add_argument("--top-n", type=int, default=10, help="Number of top pairs to export/plot.")
    p.add_argument("--output-dir", required=True, help="Output directory.")
    return p.parse_args()


def setup_logging(outdir: Path) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    log_file = outdir / "pipeline.log"
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
        handlers=[logging.FileHandler(log_file), logging.StreamHandler(sys.stdout)],
    )


def load_visium(args: argparse.Namespace):
    logging.info("Loading Visium data from %s", args.visium_dir)
    adata = sq.read.visium(path=args.visium_dir, counts_file=args.counts_file)
    adata.var_names_make_unique()

    if args.filter_in_tissue:
        if "in_tissue" not in adata.obs.columns:
            raise KeyError("--filter-in-tissue was set, but obs['in_tissue'] is missing.")
        adata = adata[adata.obs["in_tissue"] == 1].copy()

    if args.group_csv is not None:
        group_df = pd.read_csv(args.group_csv, index_col=args.index_col)
        if args.group_column not in group_df.columns:
            raise KeyError(f"Column '{args.group_column}' not found in {args.group_csv}.")
        adata.obs = adata.obs.join(group_df[[args.group_column]], how="left")
        logging.info("Joined group annotations from %s", args.group_csv)
        logging.info("Group counts:\n%s", adata.obs[args.group_column].value_counts(dropna=False).to_string())

    if args.section is not None:
        if args.group_column not in adata.obs.columns:
            raise KeyError("--section was provided, but section labels are missing from adata.obs.")
        adata = adata[adata.obs[args.group_column] == args.section].copy()
        logging.info("Subset section '%s': %d spots, %d genes", args.section, adata.n_obs, adata.n_vars)

    if adata.n_obs == 0:
        raise ValueError("No spots left after filtering/subsetting.")

    return adata


def normalize_adata(adata, args: argparse.Namespace):
    adata.raw = adata.copy()
    if args.skip_normalization:
        logging.info("Skipping normalization because --skip-normalization was set.")
        return adata

    logging.info("Filtering genes with min_counts >= %d", args.min_gene_counts)
    sc.pp.filter_genes(adata, min_counts=args.min_gene_counts)
    logging.info("Normalizing to target_sum=%s and log1p transforming", args.target_sum)
    sc.pp.normalize_total(adata, target_sum=args.target_sum)
    sc.pp.log1p(adata)
    return adata


def extract_lr_wnt1_all(adata, ligand_gene: str = "Wnt1", mean: str = "algebra"):
    if ligand_gene not in adata.var_names:
        raise ValueError(f"Ligand gene '{ligand_gene}' not found in adata.var_names.")

    receptors = adata.var_names.tolist()
    gene_inter = pd.DataFrame(
        {
            "interaction_name": [f"{ligand_gene}_{r}".upper() for r in receptors],
            "ligand": [ligand_gene] * len(receptors),
            "receptor": receptors,
            "annotation": "Secreted Signaling",
        }
    ).set_index("interaction_name")
    gene_inter.index.name = None

    adata.uns["mean"] = mean
    adata.uns["geneInter"] = gene_inter
    adata.uns["ligand"] = pd.DataFrame({"Ligand": gene_inter["ligand"]}, index=gene_inter.index)
    adata.uns["receptor"] = pd.DataFrame({"Receptor": gene_inter["receptor"]}, index=gene_inter.index)
    adata.uns["num_pairs"] = len(gene_inter)
    logging.info("Created %d custom %s-vs-all LR pairs.", len(gene_inter), ligand_gene)
    return adata



def _read_lr_csv(path_or_url: str) -> pd.DataFrame:
    if str(path_or_url).startswith(("http://", "https://")):
        r = requests.get(path_or_url, timeout=60)
        r.raise_for_status()
        if not r.text.strip():
            raise ValueError(f"Empty response from {path_or_url}")
        from io import StringIO
        return pd.read_csv(StringIO(r.text))
    return pd.read_csv(path_or_url)


def _normalize_lr_table(df: pd.DataFrame) -> pd.DataFrame:
    rename_map = {}
    for c in df.columns:
        lc = c.lower()
        if lc in {"ligand", "ligand_symbol", "source", "gene_a"}:
            rename_map[c] = "ligand"
        elif lc in {"receptor", "receptor_symbol", "target", "gene_b"}:
            rename_map[c] = "receptor"
        elif lc in {"interaction_name", "interaction", "pair", "lr_pair"}:
            rename_map[c] = "interaction_name"
        elif lc in {"annotation", "pathway_name", "category"}:
            rename_map[c] = "annotation"
    df = df.rename(columns=rename_map)

    required = {"ligand", "receptor"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(
            "Ligand-receptor CSV must contain ligand and receptor columns "
            f"(missing: {sorted(missing)}). Columns seen: {list(df.columns)}"
        )

    df = df.copy()
    df["ligand"] = df["ligand"].astype(str)
    df["receptor"] = df["receptor"].astype(str)
    if "interaction_name" not in df.columns:
        df["interaction_name"] = (df["ligand"] + "_" + df["receptor"]).str.upper()
    if "annotation" not in df.columns:
        df["annotation"] = "CellChatDB"
    return df[["interaction_name", "ligand", "receptor", "annotation"]]


def extract_lr_from_table(adata, lr_df: pd.DataFrame, min_cell: int = 3, mean: str = "algebra"):
    lr_df = _normalize_lr_table(lr_df)

    genes = pd.Index(adata.var_names)
    keep = lr_df["ligand"].isin(genes) & lr_df["receptor"].isin(genes)
    lr_df = lr_df.loc[keep].copy()

    if lr_df.empty:
        raise RuntimeError("No LR pairs overlap with adata.var_names.")

    X = adata.X
    valid_rows = []
    for _, row in lr_df.iterrows():
        l_expr = np.asarray(adata[:, row["ligand"]].X.mean(axis=1)).ravel()
        r_expr = np.asarray(adata[:, row["receptor"]].X.mean(axis=1)).ravel()
        if (l_expr > 0).sum() >= min_cell and (r_expr > 0).sum() >= min_cell:
            valid_rows.append(row)

    if not valid_rows:
        raise RuntimeError("No LR pairs passed the min_cell filter after overlap with the expression matrix.")

    gene_inter = pd.DataFrame(valid_rows).drop_duplicates(subset=["interaction_name"]).set_index("interaction_name")
    gene_inter.index.name = None

    adata.uns["mean"] = mean
    adata.uns["geneInter"] = gene_inter
    adata.uns["ligand"] = pd.DataFrame({"Ligand": gene_inter["ligand"]}, index=gene_inter.index)
    adata.uns["receptor"] = pd.DataFrame({"Receptor": gene_inter["receptor"]}, index=gene_inter.index)
    adata.uns["num_pairs"] = len(gene_inter)
    logging.info("Loaded %d LR pairs from explicit table input.", len(gene_inter))
    return adata

def prepare_lr_pairs(adata, args: argparse.Namespace):
    if args.lr_mode == "cellchat":
        if args.lr_db_csv is not None:
            logging.info("Loading LR pairs from local CSV: %s", args.lr_db_csv)
            lr_df = _read_lr_csv(args.lr_db_csv)
            extract_lr_from_table(adata, lr_df, min_cell=args.min_cell)
        elif args.lr_db_url is not None:
            logging.info("Loading LR pairs from URL: %s", args.lr_db_url)
            lr_df = _read_lr_csv(args.lr_db_url)
            extract_lr_from_table(adata, lr_df, min_cell=args.min_cell)
        else:
            logging.info("Extracting LR pairs from SpatialDM built-in CellChatDB route for species=%s", args.species)
            try:
                sdm.extract_lr(adata, args.species, min_cell=args.min_cell)
            except Exception as e:
                raise RuntimeError(
                    "SpatialDM extract_lr() failed. This often happens because the hard-coded "
                    "remote CellChatDB download returned an empty response. Re-run with "
                    "--lr-db-csv <local_db.csv> for a reproducible local database."
                ) from e
    else:
        extract_lr_wnt1_all(adata, ligand_gene=args.ligand_gene)

    if "geneInter" not in adata.uns or len(adata.uns["geneInter"]) == 0:
        raise RuntimeError("No ligand-receptor pairs available after extraction.")
    logging.info("Number of LR pairs: %d", len(adata.uns["geneInter"]))
    return adata

def build_weight_matrix(adata, args: argparse.Namespace):
    logging.info("Building weight matrix with l=%s, n_neighbors=%s, single_cell=%s",
                 args.l, args.n_neighbors, args.single_cell)
    kwargs = {"l": args.l, "single_cell": args.single_cell}
    if args.n_neighbors is not None:
        kwargs["n_neighbors"] = args.n_neighbors
    sdm.weight_matrix(adata, **kwargs)

    if "weight" not in adata.obsp or adata.obsp["weight"].nnz == 0:
        raise RuntimeError("SpatialDM weight matrix is empty. Adjust --l or --n-neighbors.")
    return adata


def run_spatialdm(adata, args: argparse.Namespace):
    logging.info("Running global selection")
    t0 = time.time()
    sdm.spatialdm_global(adata, n_perm=args.n_perm, specified_ind=None, method=args.method, nproc=args.nproc)
    sdm.sig_pairs(adata, method=args.method, fdr=True, threshold=args.pair_fdr)
    logging.info("Global selection finished in %.2f s", time.time() - t0)

    global_res = adata.uns.get("global_res")
    if global_res is None or global_res.empty:
        raise RuntimeError("SpatialDM returned an empty global_res table.")

    logging.info("Running local selection")
    t0 = time.time()
    sdm.spatialdm_local(adata, n_perm=args.n_perm, method=args.method, specified_ind=None, nproc=args.nproc)
    sdm.sig_spots(adata, method=args.method, fdr=False, threshold=args.spot_p)
    logging.info("Local selection finished in %.2f s", time.time() - t0)
    return adata


def save_summary(adata, args: argparse.Namespace, outdir: Path):
    global_res = adata.uns["global_res"].sort_values("fdr").copy()
    global_csv = outdir / "global_results.csv"
    global_res.to_csv(global_csv)

    selected = global_res.query("selected") if "selected" in global_res.columns else global_res.iloc[0:0]
    selected_csv = outdir / "significant_pairs.csv"
    selected.to_csv(selected_csv)

    top_pairs = global_res.index[: min(args.top_n, len(global_res))].tolist()
    pd.Series(top_pairs, name="pair").to_csv(outdir / "top_pairs.csv", index=False)

    meta = {
        "spatialdm_version": getattr(sdm, "__version__", "unknown"),
        "n_obs": int(adata.n_obs),
        "n_vars": int(adata.n_vars),
        "lr_mode": args.lr_mode,
        "section": args.section,
        "n_lr_pairs": int(len(adata.uns.get("geneInter", []))),
        "n_significant_pairs": int(selected.shape[0]),
        "top_pairs": top_pairs,
        "parameters": vars(args),
    }
    with open(outdir / "run_metadata.json", "w") as f:
        json.dump(meta, f, indent=2)

    logging.info("Saved %s and %s", global_csv.name, selected_csv.name)
    logging.info("Significant pairs: %d", selected.shape[0])
    return global_res, top_pairs


def save_plots(adata, global_res: pd.DataFrame, top_pairs: list[str], args: argparse.Namespace, outdir: Path):
    if len(top_pairs) == 0:
        logging.warning("No top pairs available for plotting.")
        return

    if args.ligand_gene in adata.var_names:
        fig = sq.pl.spatial_scatter(adata, color=args.ligand_gene, return_fig=True, show=False)
        fig.savefig(outdir / f"spatial_{args.ligand_gene}.png", dpi=200, bbox_inches="tight")
        plt.close(fig)

    fig = sdm_pl.global_plot(adata, pairs=top_pairs, figsize=(6, 5), vmin=-1.5, vmax=2, show=False)
    if fig is None:
        fig = plt.gcf()
    fig.savefig(outdir / "global_plot_top_pairs.png", dpi=200, bbox_inches="tight")
    plt.close(fig)

    fig = sdm_pl.plot_pairs(
        adata,
        pairs_to_plot=top_pairs,
        marker="s",
        cmap="RdPu",
        cmap_l="RdPu",
        cmap_r="RdPu",
        s=3,
        vmin=-0.1,
        show=False,
    )
    if fig is None:
        fig = plt.gcf()
    fig.savefig(outdir / "local_plot_top_pairs.png", dpi=200, bbox_inches="tight")
    plt.close(fig)

    if adata.n_obs > 50 and "weight" in adata.obsp:
        plt.figure(figsize=(4, 3))
        plt.scatter(
            adata.obsm["spatial"][:, 0],
            adata.obsm["spatial"][:, 1],
            c=np.asarray(adata.obsp["weight"][50].todense()).ravel(),
            s=3,
            cmap="RdPu",
        )
        plt.title("SpatialDM kernel influence (spot 50)")
        plt.tight_layout()
        plt.savefig(outdir / "weight_matrix_example.png", dpi=200)
        plt.close()


def main() -> int:
    args = parse_args()
    outdir = Path(args.output_dir)
    setup_logging(outdir)

    logging.info("SpatialDM version: %s", getattr(sdm, "__version__", "unknown"))
    adata = load_visium(args)
    adata = normalize_adata(adata, args)
    adata = prepare_lr_pairs(adata, args)
    adata = build_weight_matrix(adata, args)
    adata = run_spatialdm(adata, args)
    global_res, top_pairs = save_summary(adata, args, outdir)
    save_plots(adata, global_res, top_pairs[: min(args.top_n, len(top_pairs))], args, outdir)

    adata.write(outdir / "spatialdm_processed.h5ad")
    logging.info("Saved processed AnnData to %s", outdir / "spatialdm_processed.h5ad")
    logging.info("Pipeline completed successfully.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
