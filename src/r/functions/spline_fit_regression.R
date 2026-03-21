# ----------------------------
# Detection-filter sets + tables
# ----------------------------
make_filter_gene_sets <- function(
    seu_wt,
    seu_oe,
    assay_model = "SCT",
    layer_model = "data",
    filter_assay = "Spatial",
    filter_layer = "counts",
    perc_filter = 0.01
) {

  Xw_mod <- Seurat::GetAssayData(seu_wt, assay = assay_model, layer = layer_model)
  Xo_mod <- Seurat::GetAssayData(seu_oe, assay = assay_model, layer = layer_model)

  genes_model <- intersect(rownames(Xw_mod), rownames(Xo_mod))
  if (length(genes_model) == 0) {
    stop("No shared model genes between WT and OE.")
  }

  Xw_cnt <- Seurat::GetAssayData(seu_wt, assay = filter_assay, layer = filter_layer)
  Xo_cnt <- Seurat::GetAssayData(seu_oe, assay = filter_assay, layer = filter_layer)

  genes_cnt <- intersect(rownames(Xw_cnt), rownames(Xo_cnt))
  genes_base <- intersect(genes_model, genes_cnt)

  if (length(genes_base) == 0) {
    stop("No shared genes between model assay/layer and filter assay/layer.")
  }

  Xw_cnt <- Xw_cnt[genes_base, , drop = FALSE]
  Xo_cnt <- Xo_cnt[genes_base, , drop = FALSE]

  pct_wt <- Matrix::rowMeans(Xw_cnt > 0)
  pct_oe <- Matrix::rowMeans(Xo_cnt > 0)

  det_wt <- names(pct_wt)[pct_wt >= perc_filter]
  det_oe <- names(pct_oe)[pct_oe >= perc_filter]

  kept <- union(det_wt, det_oe)
  kept <- intersect(kept, genes_model)

  list(
    perc_filter = perc_filter,
    genes_model = genes_model,
    genes_base = genes_base,
    pct_wt = pct_wt,
    pct_oe = pct_oe,
    detected_wt = det_wt,
    detected_oe = det_oe,
    kept = kept,
    wt_only = setdiff(det_wt, det_oe),
    oe_only = setdiff(det_oe, det_wt),
    overlap = intersect(det_wt, det_oe)
  )
}

make_filter_tables <- function(sets) {
  summary_counts <- data.frame(
    metric = c(
      "model_genes_intersection",
      "genes_used_for_detection_calc",
      "detected_WT",
      "detected_OE",
      "overlap_detected",
      "WT_only_detected",
      "OE_only_detected",
      "kept_OR_rule"
    ),
    n = c(
      length(sets$genes_model),
      length(sets$genes_base),
      length(sets$detected_wt),
      length(sets$detected_oe),
      length(sets$overlap),
      length(sets$wt_only),
      length(sets$oe_only),
      length(sets$kept)
    )
  )

  all_genes <- sort(unique(sets$genes_base))
  membership_table <- data.frame(
    gene = all_genes,
    pct_WT = as.numeric(sets$pct_wt[all_genes]),
    pct_OE = as.numeric(sets$pct_oe[all_genes]),
    detected_WT = all_genes %in% sets$detected_wt,
    detected_OE = all_genes %in% sets$detected_oe,
    kept_OR_rule = all_genes %in% sets$kept,
    stringsAsFactors = FALSE
  )

  list(
    summary_counts = summary_counts,
    membership_table = membership_table
  )
}



# ----------------------------
# Main: distance spline limma
# ----------------------------
spline_limma_abs <- function(
    seu_wt,
    seu_oe,
    dist_col = "graphdist_to_TrabBone",
    assay = "SCT",
    layer = "data",
    nbins = NULL,
    spline_df = 6,
    perc_filter = NULL,
    max_dist = NULL,
    filter_assay = "Spatial",
    filter_layer = "counts",
    min_shared_points = 4,
    auto_cap_df = TRUE,
    return_filter_outputs = TRUE
) {

  filter_out <- NULL
  keep_genes <- NULL

  if (!is.null(perc_filter)) {
    sets <- make_filter_gene_sets(
      seu_wt,
      seu_oe,
      assay_model = assay,
      layer_model = layer,
      filter_assay = filter_assay,
      filter_layer = filter_layer,
      perc_filter = perc_filter
    )

    keep_genes <- sets$kept

    if (return_filter_outputs) {
      tabs <- make_filter_tables(sets)
      filter_out <- list(
        sets = sets,
        summary_counts = tabs$summary_counts,
        membership_table = tabs$membership_table
      )
    }
  }

  # Expression for modeling
  X_wt <- Seurat::GetAssayData(seu_wt, assay = assay, layer = layer)
  X_oe <- Seurat::GetAssayData(seu_oe, assay = assay, layer = layer)

  genes <- intersect(rownames(X_wt), rownames(X_oe))
  if (length(genes) == 0) {
    stop("No shared genes between WT and OE in modeling assay/layer.")
  }

  if (!is.null(keep_genes)) {
    genes <- intersect(genes, keep_genes)
  }

  if (length(genes) == 0) {
    stop("No genes left after detection filtering.")
  }

  X_wt <- X_wt[genes, , drop = FALSE]
  X_oe <- X_oe[genes, , drop = FALSE]

  # Distances
  r_wt <- suppressWarnings(as.integer(as.character(seu_wt[[dist_col]][, 1])))
  r_oe <- suppressWarnings(as.integer(as.character(seu_oe[[dist_col]][, 1])))

  keep_wt <- !is.na(r_wt)
  keep_oe <- !is.na(r_oe)

  X_wt <- X_wt[, keep_wt, drop = FALSE]
  X_oe <- X_oe[, keep_oe, drop = FALSE]
  r_wt <- r_wt[keep_wt]
  r_oe <- r_oe[keep_oe]

  if (length(r_wt) == 0 || length(r_oe) == 0) {
    stop("No spots left after removing NA distances.")
  }

  max_shared <- if (is.null(max_dist)) {
    min(max(r_wt, na.rm = TRUE), max(r_oe, na.rm = TRUE))
  } else {
    as.integer(max_dist)
  }

  keep_wt2 <- r_wt <= max_shared
  keep_oe2 <- r_oe <= max_shared

  X_wt <- X_wt[, keep_wt2, drop = FALSE]
  X_oe <- X_oe[, keep_oe2, drop = FALSE]
  r_wt <- r_wt[keep_wt2]
  r_oe <- r_oe[keep_oe2]

  # Aggregate means by distance / bin
  aggregate_means <- function(X, r, dist_levels) {
    n <- vapply(dist_levels, function(d) sum(r == d), numeric(1))

    Y <- vapply(
      dist_levels,
      function(d) {
        cols <- which(r == d)
        if (length(cols) == 0) {
          rep(NA_real_, nrow(X))
        } else {
          Matrix::rowMeans(X[, cols, drop = FALSE])
        }
      },
      FUN.VALUE = numeric(nrow(X))
    )

    Y <- t(Y)
    Y <- t(Y)

    list(Y = Y, n = n)
  }

  if (!is.null(nbins)) {
    nbins <- as.integer(nbins)
    if (nbins < 2) {
      stop("nbins must be >= 2.")
    }

    brks <- seq(0, max_shared, length.out = nbins + 1)

    bin_wt <- cut(
      r_wt,
      breaks = brks,
      include.lowest = TRUE,
      labels = FALSE
    )

    bin_oe <- cut(
      r_oe,
      breaks = brks,
      include.lowest = TRUE,
      labels = FALSE
    )

    dist_levels <- seq_len(nbins)

    agg_wt <- aggregate_means(X_wt, bin_wt, dist_levels)
    agg_oe <- aggregate_means(X_oe, bin_oe, dist_levels)

    Y_wt <- agg_wt$Y
    w_wt <- agg_wt$n
    Y_oe <- agg_oe$Y
    w_oe <- agg_oe$n

    # numeric positions for spline fitting (unchanged statistical role)
    dist_axis <- seq(0, max_shared, length.out = nbins)

    # shared modeled bins only
    ok <- which((w_wt > 0) & (w_oe > 0))
    if (length(ok) < min_shared_points) {
      stop(
        "Too few shared distance points after filtering: ",
        length(ok),
        " (need >= ",
        min_shared_points,
        ")."
      )
    }

    # exact integer labels from observed values in the retained shared bins
    bin_labels_ok <- vapply(ok, function(i) {
      vals_wt <- r_wt[bin_wt == i]
      vals_oe <- r_oe[bin_oe == i]
      vals <- sort(unique(c(vals_wt, vals_oe)))

      if (length(vals) == 0) {
        paste0("bin", i)
      } else {
        paste0(min(vals), "-", max(vals))
      }
    }, character(1))

    Y_wt <- Y_wt[, ok, drop = FALSE]
    Y_oe <- Y_oe[, ok, drop = FALSE]
    Y <- cbind(Y_wt, Y_oe)

    w <- as.numeric(c(w_wt[ok], w_oe[ok]))
    w[w <= 0] <- 1

    meta <- data.frame(
      group = stats::relevel(
        factor(c(rep("Control OVX", ncol(Y_wt)), rep("Wnt1tg OVX", ncol(Y_oe)))),
        ref = "Control OVX"
      ),
      dist = c(dist_axis[ok], dist_axis[ok]),
      dist_label = c(bin_labels_ok, bin_labels_ok),
      weight = w
    )

    dist_axis_ok <- dist_axis[ok]
    dist_labels_ok <- bin_labels_ok

  } else {
    dist_axis <- 0:max_shared
    dist_levels <- dist_axis

    agg_wt <- aggregate_means(X_wt, r_wt, dist_levels)
    agg_oe <- aggregate_means(X_oe, r_oe, dist_levels)

    Y_wt <- agg_wt$Y
    w_wt <- agg_wt$n
    Y_oe <- agg_oe$Y
    w_oe <- agg_oe$n

    ok <- which((w_wt > 0) & (w_oe > 0))
    if (length(ok) < min_shared_points) {
      stop(
        "Too few shared distance points after filtering: ",
        length(ok),
        " (need >= ",
        min_shared_points,
        ")."
      )
    }

    Y_wt <- Y_wt[, ok, drop = FALSE]
    Y_oe <- Y_oe[, ok, drop = FALSE]
    Y <- cbind(Y_wt, Y_oe)

    w <- as.numeric(c(w_wt[ok], w_oe[ok]))
    w[w <= 0] <- 1

    dist_labels_ok <- as.character(dist_axis[ok])

    meta <- data.frame(
      group = stats::relevel(
        factor(c(rep("Control OVX", ncol(Y_wt)), rep("Wnt1tg OVX", ncol(Y_oe)))),
        ref = "Control OVX"
      ),
      dist = c(dist_axis[ok], dist_axis[ok]),
      dist_label = c(dist_labels_ok, dist_labels_ok),
      weight = w
    )

    dist_axis_ok <- dist_axis[ok]
  }

  # Cap df
  K <- length(unique(meta$dist))
  if (auto_cap_df) {
    df_cap <- max(
      2L,
      min(
        as.integer(spline_df),
        as.integer(K - 1L),
        as.integer(floor(K / 3))
      )
    )

    if (df_cap != spline_df) {
      message(
        "Adjusted spline_df from ",
        spline_df,
        " to ",
        df_cap,
        " based on K=",
        K,
        " shared distance points."
      )
    }

    spline_df <- df_cap
  } else {
    if (spline_df >= K) {
      stop("spline_df must be < K (K=", K, ").")
    }
  }

  # Spline basis + design
  B_train <- splines::ns(meta$dist, df = spline_df)
  knots <- attr(B_train, "knots")
  bknots <- attr(B_train, "Boundary.knots")

  B <- splines::ns(meta$dist, knots = knots, Boundary.knots = bknots)
  colnames(B) <- paste0("b", seq_len(ncol(B)))
  metaB <- cbind(meta, B)

  form <- stats::as.formula(
    paste0("~ group * (", paste(colnames(B), collapse = "+"), ")")
  )
  design <- stats::model.matrix(form, data = metaB)

  # Fit
  fit <- limma::eBayes(
    limma::lmFit(Y, design, weights = metaB$weight)
  )

  cn <- colnames(design)

  int_terms <- grep("^groupOE:b", cn, value = TRUE)
  if (length(int_terms) == 0) {
    int_terms <- grep("^group.*:b", cn, value = TRUE)
  }
  if (length(int_terms) == 0) {
    stop("Could not find interaction terms for shape test in design matrix.")
  }

  main_group <- grep("^groupOE$", cn, value = TRUE)
  if (length(main_group) == 0) {
    main_group <- grep("^group", cn, value = TRUE)
    main_group <- main_group[!grepl(":", main_group)]
    main_group <- setdiff(main_group, "(Intercept)")
  }
  if (length(main_group) == 0) {
    stop("Could not find main group term for shift test in design matrix.")
  }

  tt_shape <- limma::topTable(
    fit,
    coef = match(int_terms, cn),
    number = Inf,
    sort.by = "F"
  )

  tt_shift <- limma::topTable(
    fit,
    coef = match(main_group[1], cn),
    number = Inf,
    sort.by = "P"
  )

  list(
    fit = fit,
    design = design,
    meta = metaB,
    spline = list(
      df = spline_df,
      knots = knots,
      bknots = bknots
    ),
    shape = tt_shape,
    shift = tt_shift,
    filter = filter_out,
    info = list(
      K_shared = K,
      ok_points = ok,
      dist_axis_ok = dist_axis_ok,
      dist_labels_ok = dist_labels_ok,
      max_shared = max_shared,
      nbins = nbins,
      assay = assay,
      layer = layer,
      filter_assay = filter_assay,
      filter_layer = filter_layer,
      perc_filter = perc_filter
    )
  )
}
