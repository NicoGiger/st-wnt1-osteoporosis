# Libraries ####

library(splines)
library(ggplot2)
library(pheatmap)
library(ggplot2)
library(msigdbr)

# database ####

# Reactome (CP:REACTOME)
db_reactome <- msigdbr(species = "Mus musculus") %>%
  filter(gs_subcollection == "CP:REACTOME") %>%   # NOTE: some msigdbr versions use gs_subcat; others use gs_subcollection
  transmute(term = gs_name, gene = gene_symbol) %>%
  distinct()
db_go <- msigdbr(species = "Mus musculus") %>%
  filter(gs_subcollection == "GO:BP") %>%   # NOTE: some msigdbr versions use gs_subcat; others use gs_subcollection
  transmute(term = gs_name, gene = gene_symbol) %>%
  distinct()
db_kegg <- msigdbr(species = "Mus musculus") %>%
  filter(gs_subcollection == "CP:KEGG") %>%
  transmute(term = gs_name, gene = gene_symbol) %>%
  distinct()
db_kegg <- msigdbr(species = "Mus musculus", category = "C2", subcollection = "CP:KEGG_LEGACY") %>%
  transmute(term = gs_name, gene = gene_symbol) %>%
  distinct()
# Hallmark (H)
db_hallmark <- msigdbr(species = "Mus musculus", category = "H") %>%
  transmute(term = gs_name, gene = gene_symbol) %>%
  distinct()

# Functions ####

plot_heatmap_original_rings <- function(res, top_n = 40, scale_rows = TRUE) {

  stopifnot(!is.null(res$fit), !is.null(res$design), !is.null(res$meta))

  # --- Top genes by shape ---
  tt <- res$shape[order(res$shape$F, decreasing = TRUE), , drop = FALSE]
  genes <- rownames(tt)[seq_len(min(top_n, nrow(tt)))]

  meta <- res$meta

  # Extract ring numbers and order them 0 → max
  rings <- sort(unique(meta$dist))

  # Get WT and OE design rows ordered by ring
  ord_wt <- unlist(lapply(rings, function(r)
    which(meta$group == "WT" & meta$dist == r)
  ))

  ord_oe <- unlist(lapply(rings, function(r)
    which(meta$group == "OE" & meta$dist == r)
  ))

  design_wt <- res$design[ord_wt, , drop = FALSE]
  design_oe <- res$design[ord_oe, , drop = FALSE]

  cn <- colnames(res$design)

  genes <- intersect(genes, rownames(res$fit$coefficients))
  beta  <- res$fit$coefficients[genes, cn, drop = FALSE]

  # Predicted fitted values
  WT   <- beta %*% t(design_wt)
  OE   <- beta %*% t(design_oe)
  DIFF <- OE - WT

  # --- Label columns exactly with ring numbers ---
  colnames(WT)   <- rings
  colnames(OE)   <- rings
  colnames(DIFF) <- rings

  # Label rows with gene names
  rownames(WT)   <- genes
  rownames(OE)   <- genes
  rownames(DIFF) <- genes

  sc <- if (scale_rows) "row" else "none"

  pheatmap(WT,
           scale = sc,
           main = "WT fitted (SCT residuals)",
           cluster_cols = FALSE, cluster_rows = FALSE)

  pheatmap(OE,
           scale = sc,
           main = "OE fitted (SCT residuals)",
           cluster_cols = FALSE, cluster_rows = FALSE)

  pheatmap(DIFF,
           scale = sc,
           main = "OE − WT difference",
           cluster_cols = FALSE, clustering_distance_rows = "correlation",
           clustering_method = "average")
}

volcano_shift <- function(res, GOIs = character(), filename = NULL,
                          top_pct = 0.05,            # top 5% |logFC| defines vertical cutoffs
                          pval_sig = 0.01,           # "UP/DOWN" adjusted p value threshold,
                          pval_sug = 0.05,           # "SUGGESTIVE" pvalue (not adjusted) threshold (dashed line shown)
                          use_adj_for_y = FALSE,     # use adj.P.Val on y-axis if TRUE
                          force = 1.5,
                          force_pull = 0.5,
                          max.overlaps = 10,
                          n_label = 10,
                          lim = FALSE,
                          xlim_abs = 3,
                          ylim = 100) {

  tt <- res$shift
  if (is.null(tt) || nrow(tt) == 0) stop("res$shift is missing or empty.")
  tt <- tibble::as_tibble(tt, rownames = "gene")

  pcol <- if (use_adj_for_y) "adj.P.Val" else "P.Value"
  if (!(pcol %in% colnames(tt))) stop("Column missing in res$shift: ", pcol)

  tt <- tt %>%
    mutate(
      pvalue = .data[[pcol]],
      y = -log10(pmax(pvalue, .Machine$double.xmin))
    )

  # Effect cutoff from data (top X% absolute)
  thr <- unname(stats::quantile(abs(tt$logFC), probs = 1 - top_pct, na.rm = TRUE))

  # DEG categories aligned with the dashed lines:
  # - UP/DOWN: beyond effect cutoff AND p < pval_sig
  # - SUGGESTIVE: beyond effect cutoff AND p < pval_sug (but not < pval_sig)
  tt <- tt %>%
    mutate(
      DEG = case_when(
        logFC >=  thr & pvalue < pval_sig ~ "UP",
        logFC <= -thr & pvalue < pval_sig ~ "DOWN",
        logFC >=  thr & pvalue < pval_sug ~ "SUGGESTIVE UP",
        logFC <= -thr & pvalue < pval_sug ~ "SUGGESTIVE DOWN",
        TRUE ~ "NO"
      ),
      DEG = factor(DEG, levels = c("DOWN","SUGGESTIVE DOWN","NO","SUGGESTIVE UP","UP"))
    )

  # label: top 5 UP and DOWN (by pvalue) + GOIs
  genes_up <- tt %>% filter(DEG == "UP") %>% slice_min(order_by = pvalue, n = n_label) %>% pull(gene)
  genes_down <- tt %>% filter(DEG == "DOWN") %>% slice_min(order_by = pvalue, n = n_label) %>% pull(gene)
  genes_to_label <- unique(c(genes_up, genes_down, GOIs))

  plt <- ggplot(tt, aes(x = logFC, y = y, col = DEG)) +
    geom_point() +
    scale_color_manual(values = c(
      "DOWN" = "purple",
      "NO" = "grey",
      "UP" = "orange",
      "SUGGESTIVE DOWN" = "#DB7093",
      "SUGGESTIVE UP" = "#FFD700"
    )) +
    theme(text = element_text(size = 20)) +
    geom_vline(xintercept = c(-thr, thr), linetype = "dashed", col = "black") +
    geom_hline(yintercept = -log10(pval_sug), linetype = "dashed", col = "black") +
    geom_label_repel(
      data = tt %>% filter(gene %in% genes_to_label),
      aes(label = gene),
      label.size = 0,
      force = force,
      force_pull = force_pull,
      max.overlaps = max.overlaps
    ) +
    xlab("Effect size (SCT residual units)") +
    ylab(paste0("p-value (-log10) [", pcol, "]")) +
    theme_bw()

  if (isTRUE(lim)) {
    plt <- plt + xlim(c(-xlim_abs, xlim_abs)) + ylim(c(0, ylim))
  }

  if (!is.null(filename)) {
    ggsave(filename, plt, width = 7, height = 6, dpi = 300)
  }

  return(plt)
}

predict_gene_curves <- function(res, gene, ngrid = 200) {
  if (is.null(res$meta) || nrow(res$meta) == 0) stop("res$meta is missing/empty.")
  if (!gene %in% rownames(res$fit$coefficients)) stop("Gene not in fit: ", gene)

  dist_train <- res$meta$dist
  grid <- seq(min(dist_train), max(dist_train), length.out = ngrid)

  B <- ns(grid, knots = res$spline$knots, Boundary.knots = res$spline$bknots)
  colnames(B) <- paste0("b", seq_len(ncol(B)))

  new_WT <- cbind(data.frame(group = factor("WT", levels = levels(res$meta$group)), dist = grid), B)
  new_OE <- cbind(data.frame(group = factor("OE", levels = levels(res$meta$group)), dist = grid), B)

  form <- as.formula(paste0("~ group * (", paste(colnames(B), collapse = "+"), ")"))
  X_WT <- model.matrix(form, data = new_WT)
  X_OE <- model.matrix(form, data = new_OE)

  cn <- colnames(res$design)
  X_WT <- X_WT[, cn, drop = FALSE]
  X_OE <- X_OE[, cn, drop = FALSE]

  beta <- res$fit$coefficients[gene, cn]
  data.frame(
    dist = grid,
    WT = as.numeric(X_WT %*% beta),
    OE = as.numeric(X_OE %*% beta)
  )
}

plot_spline_curves <- function(curves, gene) {
  df1 <- rbind(
    data.frame(dist = curves$dist, group = "WT", y = curves$WT),
    data.frame(dist = curves$dist, group = "OE", y = curves$OE)
  )

  p_curves <- ggplot(df1, aes(dist, y, linetype = group)) +
    geom_line(linewidth = 0.9) +
    theme_classic() +
    labs(title = gene, x = "Absolute ring distance", y = "Fitted expression (spline)")

  df2 <- data.frame(dist = curves$dist, diff = curves$OE - curves$WT)
  p_diff <- ggplot(df2, aes(dist, diff)) +
    geom_hline(yintercept = 0, linewidth = 0.4) +
    geom_line(linewidth = 0.9) +
    theme_classic() +
    labs(title = paste0(gene, " (OE - WT)"), x = "Absolute ring distance", y = "Difference")

  list(curves = p_curves, diff = p_diff)
}



#gsea #####

## Signed GSEA for "shape" (spline interaction) using:
##   sign = mean fitted (OE - WT) across distance
##   magnitude = sqrt(F) from the shape F-test
##
## Requirements:
## - Your `res` must come from the spline function that stores:
##     res$fit, res$design, res$meta (with group levels), res$spline$knots, res$spline$bknots
## - limma fit was done with design ~ group * (b1 + b2 + ...)


# ---- 0) Build gene lists ----

# SHIFT ranking (signed)
gl_shift <- res$shift$t
names(gl_shift) <- rownames(res$shift)
gl_shift <- gl_shift[is.finite(gl_shift)]
gl_shift <- sort(gl_shift, decreasing = TRUE)
gl_shift %>% head()
#SHAPE ranking (unsigned)
gl_shape_unsigned <- res$shape$F
names(gl_shape_unsigned) <- rownames(res$shape)
gl_shape_unsigned <- gl_shape_unsigned[is.finite(gl_shape_unsigned)]
gl_shape_unsigned <- sort(gl_shape_unsigned, decreasing = TRUE)

# ---- 1) TERM2GENE: Reactome + Hallmark (mouse) ----



# Choose one mapping to run (set to db_reactome or db_hallmark)
TERM2GENE_use <- db_reactome   # change to db_hallmark if you want Hallmark

# ---- 2) Run GSEA twice ----

gsea_shape <- GSEA(
  geneList = gl_shape_unsigned,
  TERM2GENE = TERM2GENE_use,
  pvalueCutoff = 0.05,
  seed = TRUE,
  nPermSimple = 10000,
  scoreType = "pos"
)

gsea_shift <- GSEA(
  geneList = gl_shift,
  TERM2GENE = TERM2GENE_use,
  pvalueCutoff = 0.05,
  seed = TRUE,
  nPermSimple = 10000
)


# ---- 3) Dotplots ----
gsea_shape@result %>% filter(p.adjust < 0.05) %>% pull(Description)
p_shape <- dotplot(gsea_shape , showCategory = 10) +
  theme_classic() +
  labs(title = "GSEA (SHAPE): spatial pattern changes")

cats <- c("REACTOME_ERYTHROCYTES_TAKE_UP_CARBON_DIOXIDE_AND_RELEASE_OXYGEN",
  "REACTOME_CELL_CYCLE_MITOTIC", "REACTOME_CELLULAR_SENESCENCE",
  "REACTOME_SIGNALING_BY_NOTCH", "REACTOME_SIGNALING_BY_WNT")


p_shift <- dotplot(gsea_shift, showCategory = 10, split = ".sign") + # split=".sign" separates positive vs negative NES.
  facet_grid(. ~ .sign) +
  theme_classic() +
  labs(title = "GSEA (SHIFT): overall expression difference")
gsea_shift@result %>% filter(p.adjust < 0.05) %>% pull(Description)

p_shape
p_shift

pathway_name <- "REACTOME_SIGNALING_BY_WNT"  # adjust to your exact term

le_string <- gsea_shift@result$core_enrichment[
  gsea_shape@result$Description == pathway_name
]
genes <- strsplit(le_string, "/")[[1]]

# Shift Volcano ####

gois <- c("Wnt1","Col1a2", "Col1a1", "Runx2", "Sost", "Dmp1", "Col1a1", "Mmp9",
  "Bglap", 'Lrp5')
osteoblast_genes <- c(
  "Runx2",
  "Sp7",
  "Alpl",
  "Col1a1",
  "Ibsp",
  "Bglap",
  "Spp1",
  "Dmp1",
  "Sost"
)
shift_volcano <- volcano_shift(res, GOIs = union(osteoblast_genes, gois), n_label = 5,
                               top_pct = 0.05, pval_sug = 0.05, lim = F,
                               force_pull = 10, max.overlaps = 20,)
shift_volcano
# Shape heatmap ####
plot_heatmap_original_rings(res, top_n = 75, scale_rows = TRUE)

# DistPlots (specific genes) ####
gene <- "Mmp9"
res$shape[gene,]
res$shift[gene,]
SpatialFeaturePlot(dd.norm, features = gene, alpha = c(0,3))
cur <- predict_gene_curves(res, gene)
pp <- plot_spline_curves(cur, gene)
pp$curves
pp$diff

