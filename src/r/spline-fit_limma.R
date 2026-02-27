# Libraries ####

library(Seurat)
library(limma)
library(splines)

get_logdata <- function(seu, assay = NULL, layer = "data") {
  if (is.null(assay)) assay <- DefaultAssay(seu)
  GetAssayData(seu, assay = assay, layer = layer)
}


ring_spline_limma_abs <- function(
    seu_wt,
    seu_oe,
    ring_col  = "ring",
    assay     = "SCT",
    layer     = "data",
    nbins     = NULL,     # NULL = per-ring means; or set e.g. 10/20
    spline_df = 6,
    perc_filter = NULL,    # e.g. 0.01 for 1% detection filter using RNA counts; NULL = no filter
    max_rings = NULL
) {
  suppressPackageStartupMessages({
    library(limma)
    library(splines)
    library(Matrix)
  })

  # --- Optional gene prefilter: detected in >= perc_filter spots in either sample (RNA counts) ---
  if (!is.null(perc_filter)) {
    Xw_cnt <- GetAssayData(seu_wt, assay = assay, layer = layer)
    Xo_cnt <- GetAssayData(seu_oe, assay = assay, layer = layer)
    g0 <- intersect(rownames(Xw_cnt), rownames(Xo_cnt))
    Xw_cnt <- Xw_cnt[g0, , drop = FALSE]
    Xo_cnt <- Xo_cnt[g0, , drop = FALSE]
    pct_wt <- Matrix::rowMeans(Xw_cnt > 0)
    pct_oe <- Matrix::rowMeans(Xo_cnt > 0)
    keep_detect <- (pct_wt >= perc_filter) | (pct_oe >= perc_filter)
    keep_genes <- g0[keep_detect]
  } else {
    keep_genes <- NULL
  }

  # --- Expression used for modeling ---
  X_wt <- GetAssayData(seu_wt, assay = assay, layer = layer)
  X_oe <- GetAssayData(seu_oe, assay = assay, layer = layer)

  genes <- intersect(rownames(X_wt), rownames(X_oe))
  if (!is.null(keep_genes)) genes <- intersect(genes, keep_genes)

  X_wt <- X_wt[genes, , drop = FALSE]
  X_oe <- X_oe[genes, , drop = FALSE]

  # --- Rings / absolute distance ---
  r_wt <- as.integer(as.character(seu_wt[[ring_col]][, 1]))
  r_oe <- as.integer(as.character(seu_oe[[ring_col]][, 1]))

  keep_wt <- !is.na(r_wt)
  keep_oe <- !is.na(r_oe)

  X_wt <- X_wt[, keep_wt, drop = FALSE]
  X_oe <- X_oe[, keep_oe, drop = FALSE]
  r_wt <- r_wt[keep_wt]
  r_oe <- r_oe[keep_oe]

  if (is.null(max_rings)){
    max_shared <- min(max(r_wt, na.rm = TRUE), max(r_oe, na.rm = TRUE))
  }
  else{
    max_shared <- max_rings
  }
  keep_wt2 <- r_wt <= max_shared
  keep_oe2 <- r_oe <= max_shared


  X_wt <- X_wt[, keep_wt2, drop = FALSE]
  X_oe <- X_oe[, keep_oe2, drop = FALSE]
  r_wt <- r_wt[keep_wt2]
  r_oe <- r_oe[keep_oe2]

  # --- Aggregate to bins or per-ring + compute weights (bin sizes) ---
  if (!is.null(nbins)) {
    brks <- seq(0, max_shared, length.out = nbins + 1)
    bin_wt <- cut(r_wt, breaks = brks, include.lowest = TRUE, labels = FALSE)
    bin_oe <- cut(r_oe, breaks = brks, include.lowest = TRUE, labels = FALSE)

    n_wt <- sapply(seq_len(nbins), function(b) sum(bin_wt == b))
    n_oe <- sapply(seq_len(nbins), function(b) sum(bin_oe == b))

    Y_wt <- sapply(seq_len(nbins), function(b) {
      cols <- which(bin_wt == b)
      if (length(cols) == 0) return(rep(NA_real_, nrow(X_wt)))
      rowMeans(X_wt[, cols, drop = FALSE])
    })
    Y_oe <- sapply(seq_len(nbins), function(b) {
      cols <- which(bin_oe == b)
      if (length(cols) == 0) return(rep(NA_real_, nrow(X_oe)))
      rowMeans(X_oe[, cols, drop = FALSE])
    })

    Y_wt <- t(Y_wt); Y_wt <- t(Y_wt)  # genes x bins
    Y_oe <- t(Y_oe); Y_oe <- t(Y_oe)

    dist_axis <- seq(0, max_shared, length.out = nbins)
    w_wt <- n_wt
    w_oe <- n_oe

  } else {
    rings <- 0:max_shared

    n_wt <- sapply(rings, function(rr) sum(r_wt == rr))
    n_oe <- sapply(rings, function(rr) sum(r_oe == rr))

    Y_wt <- sapply(rings, function(rr) {
      cols <- which(r_wt == rr)
      if (length(cols) == 0) return(rep(NA_real_, nrow(X_wt)))
      rowMeans(X_wt[, cols, drop = FALSE])
    })
    Y_oe <- sapply(rings, function(rr) {
      cols <- which(r_oe == rr)
      if (length(cols) == 0) return(rep(NA_real_, nrow(X_oe)))
      rowMeans(X_oe[, cols, drop = FALSE])
    })

    Y_wt <- t(Y_wt); Y_wt <- t(Y_wt)  # genes x rings
    Y_oe <- t(Y_oe); Y_oe <- t(Y_oe)

    dist_axis <- rings
    w_wt <- n_wt
    w_oe <- n_oe
  }

  # Keep only distance points present in BOTH groups (and thus weights > 0)
  ok <- which(colSums(is.na(rbind(Y_wt, Y_oe))) == 0)
  if (length(ok) < 4) stop("Too few shared distance points after filtering.")

  Y_wt <- Y_wt[, ok, drop = FALSE]
  Y_oe <- Y_oe[, ok, drop = FALSE]
  Y <- cbind(Y_wt, Y_oe)

  # weights per observation (bin/ring mean), same ordering as columns in Y
  w <- c(w_wt[ok], w_oe[ok])
  w <- as.numeric(w)
  # avoid zero weights (shouldn't happen if ok chosen correctly, but safe)
  w[w <= 0] <- NA_real_

  meta <- data.frame(
    group = relevel(factor(c(rep("WT", ncol(Y_wt)), rep("OE", ncol(Y_oe)))), ref = "WT"),
    dist  = c(dist_axis[ok], dist_axis[ok]),
    weight = w
  )

  # --- Spline basis ---
  B_train <- ns(meta$dist, df = spline_df)
  knots  <- attr(B_train, "knots")
  bknots <- attr(B_train, "Boundary.knots")

  B <- ns(meta$dist, knots = knots, Boundary.knots = bknots)
  colnames(B) <- paste0("b", seq_len(ncol(B)))
  metaB <- cbind(meta, B)

  form <- as.formula(paste0("~ group * (", paste(colnames(B), collapse = "+"), ")"))
  design <- model.matrix(form, data = metaB)

  # --- Weighted limma fit ---
  fit <- eBayes(lmFit(Y, design, weights = metaB$weight))

  cn <- colnames(design)
  int_terms  <- grep("^group.*:b", cn, value = TRUE)
  main_group <- grep("^group", cn, value = TRUE)
  main_group <- main_group[!grepl(":", main_group)]

  tt_shape <- topTable(fit, coef = match(int_terms, cn), number = Inf, sort.by = "F")
  tt_shift <- topTable(fit, coef = main_group, number = Inf, sort.by = "P")

  list(
    fit = fit,
    design = design,
    meta = metaB,  # includes weights
    spline = list(df = spline_df, knots = knots, bknots = bknots),
    shape = tt_shape,
    shift = tt_shift
  )
}




# Example
dd.norm_split <- SplitObject(dd.norm, split.by = "group")
seu_wt <- dd.norm_split$Wt
seu_oe <- dd.norm_split$WntOverExp
colnames(seu_oe) %>% length() #519
colnames(seu_wt) %>% length() #351
res <- ring_spline_limma_abs(seu_wt, seu_oe, ring_col="ring",
                            assay="SCT", layer="data", nbins=NULL, spline_df=14,
                            max_rings = 16, perc_filter = 0.05)
perc1 <- res$shape %>% rownames()
Genes_in <- c('Wnt1', 'Sost', 'Runx2', 'Col1a1', 'Col1a2', 'Fzd10', 'Fzd4', 'Fzd1', 'Fzd7', 'Fzd9', 'Lrp5', 'Lrp6')
Genes_in %in% perc1
#temp <- FindSpatiallyVariableFeatures(seu_wt, assay = 'SCT', selection.method = 'moransi', nfeatures = 3000)
#temp_oe <- FindSpatiallyVariableFeatures(seu_oe, assay = 'SCT', selection.method = 'moransi', nfeatures = 3000)
#temp_union <- union(SpatiallyVariableFeatures(temp), SpatiallyVariableFeatures(temp_oe))
#genes_filtered <- union(perc1, temp_union)
#
#Genes_in %in% perc1
#Genes_in %in% temp_union
#Genes_in %in% genes_filtered
#
#
#
