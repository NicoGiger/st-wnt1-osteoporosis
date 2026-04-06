shape_heatmap <- function(res, top_n = 40, scale_rows = TRUE) {

  stopifnot(!is.null(res$fit), !is.null(res$design), !is.null(res$meta))

  # --- Top genes by shape ---
  tt <- res$shape[order(res$shape$F, decreasing = TRUE), , drop = FALSE]
  genes <- rownames(tt)[seq_len(min(top_n, nrow(tt)))]

  meta <- res$meta

  # Extract dist numbers and order them
  dists <- sort(unique(meta$dist))

  # Use bin labels if present, otherwise fall back to numeric dist
  dist_labels <- vapply(
    dists,
    function(r) {
      lab <- unique(meta$dist_label[meta$dist == r])
      if (length(lab) > 0 && !all(is.na(lab))) {
        as.character(lab[1])
      } else {
        as.character(r)
      }
    },
    character(1)
  )

  # Get WT and OE design rows ordered by dist
  ord_wt <- unlist(lapply(dists, function(r) {
    which(meta$group == "Control OVX" & meta$dist == r)
  }))

  ord_oe <- unlist(lapply(dists, function(r) {
    which(meta$group == "Wnt1tg OVX" & meta$dist == r)
  }))

  design_wt <- res$design[ord_wt, , drop = FALSE]
  design_oe <- res$design[ord_oe, , drop = FALSE]

  cn <- colnames(res$design)

  genes <- intersect(genes, rownames(res$fit$coefficients))
  beta <- res$fit$coefficients[genes, cn, drop = FALSE]

  # Predicted fitted values
  WT <- beta %*% t(design_wt)
  OE <- beta %*% t(design_oe)
  DIFF <- OE - WT

  # Label columns with bin labels
  colnames(WT) <- dist_labels
  colnames(OE) <- dist_labels
  colnames(DIFF) <- dist_labels

  # Label rows with gene names
  rownames(WT) <- genes
  rownames(OE) <- genes
  rownames(DIFF) <- genes

  sc <- if (scale_rows) "row" else "none"

  pheatmap::pheatmap(
    WT,
    scale = sc,
    main = "Control OVX fitted (SCT residuals)",
    cluster_cols = FALSE,
    cluster_rows = FALSE,
    angle_col = 0
  )

  pheatmap::pheatmap(
    OE,
    scale = sc,
    main = "Wnt1tg OVX fitted (SCT residuals)",
    cluster_cols = FALSE,
    cluster_rows = FALSE,
    angle_col = 0
  )

  pheatmap::pheatmap(
    DIFF,
    scale = sc,
    main = "Wnt1 OVX - ctl OVX difference",
    cluster_cols = FALSE,
    clustering_distance_rows = "correlation",
    clustering_method = "average",
    angle_col = 0
  )
}
