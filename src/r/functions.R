make_signed_shape_stat <- function(res, ngrid = 200) {

  # Basic checks
  stopifnot(!is.null(res$fit))
  stopifnot(!is.null(res$design))
  stopifnot(!is.null(res$meta))
  stopifnot(!is.null(res$spline))
  stopifnot(!is.null(res$shape))

  cn <- colnames(res$design)
  genes <- rownames(res$fit$coefficients)

  # Distance grid across fitted range
  dist_train <- res$meta$dist
  grid <- seq(min(dist_train), max(dist_train), length.out = ngrid)

  # Rebuild spline basis using SAME knots
  B <- ns(grid,
          knots = res$spline$knots,
          Boundary.knots = res$spline$bknots)

  colnames(B) <- paste0("b", seq_len(ncol(B)))

  # Build design rows for WT and OE
  lvl <- levels(res$meta$group)

  new_WT <- cbind(
    data.frame(group = factor("WT", levels = lvl), dist = grid),
    B
  )

  new_OE <- cbind(
    data.frame(group = factor("OE", levels = lvl), dist = grid),
    B
  )

  form <- as.formula(paste0("~ group * (", paste(colnames(B), collapse = "+"), ")"))

  X_WT <- model.matrix(form, data = new_WT)
  X_OE <- model.matrix(form, data = new_OE)

  # Align to fitted design
  X_WT <- X_WT[, cn, drop = FALSE]
  X_OE <- X_OE[, cn, drop = FALSE]

  # Average contrast vector (OE − WT across distance)
  contrast_vec <- colMeans(X_OE - X_WT)

  # Gene-specific directional effect
  beta <- res$fit$coefficients[, cn, drop = FALSE]
  effect <- as.numeric(beta %*% contrast_vec)
  names(effect) <- genes

  # F-statistic from shape test
  Fstat <- res$shape$F
  names(Fstat) <- rownames(res$shape)

  # Align genes
  common <- intersect(names(effect), names(Fstat))

  signed_stat <- sign(effect[common]) * sqrt(Fstat[common])

  # Clean + sort for GSEA
  signed_stat <- signed_stat[is.finite(signed_stat)]
  signed_stat <- sort(signed_stat, decreasing = TRUE)

  signed_stat
}



# ----------------------------
# Helpers
# ----------------------------
top_shape_genes <- function(res, n = 20) {
  tt <- res$shape
  tt <- tt[order(tt$F, decreasing = TRUE), , drop = FALSE]
  rownames(tt)[seq_len(min(n, nrow(tt)))]
}

# Mean expression per (absolute) ring for one Seurat object and one gene (RNA log-norm)
ring_profile_gene <- function(seu, gene, ring_col = "ring", assay = "RNA", layer = "data") {
  x <- GetAssayData(seu, assay = assay, layer = layer)[gene, ]
  r <- as.integer(as.character(seu[[ring_col]][, 1]))
  df <- data.frame(ring = r, expr = as.numeric(x))
  df <- df[!is.na(df$ring), ]
  prof <- aggregate(expr ~ ring, df, mean)
  prof[order(prof$ring), ]
}

# Overlay curves WT vs OE for one gene (RNA log-norm, averaged per ring)
plot_ring_overlay <- function(seu_wt, seu_oe, gene,
                              ring_col = "ring",
                              assay = "Spatial", layer = "data",
                              max_shared = NULL) {
  p1 <- ring_profile_gene(seu_wt, gene, ring_col, assay, layer) %>% mutate(group = "WT")
  p2 <- ring_profile_gene(seu_oe, gene, ring_col, assay, layer) %>% mutate(group = "OE")
  df <- bind_rows(p1, p2)

  if (is.null(max_shared)) {
    max_shared <- min(max(p1$ring), max(p2$ring))
  }
  df <- df[df$ring <= max_shared, ]

  ggplot(df, aes(ring, expr, linetype = group)) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2) +
    theme_classic() +
    labs(title = gene, x = "Absolute ring distance", y = paste0(assay, " ", layer, " (mean per ring)"))
}

plot_ring_overlay(seu_wt = seu_wt, seu_oe = seu_oe, assay = "SCT")
