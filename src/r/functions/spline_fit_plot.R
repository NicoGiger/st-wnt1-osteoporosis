observed_dist_means <- function(seu_wt, seu_oe, gene,
                                dist_col = "dist",
                                assay = "SCT", layer = "data") {
  Xw <- Seurat::GetAssayData(seu_wt, assay = assay, layer = layer)
  Xo <- Seurat::GetAssayData(seu_oe, assay = assay, layer = layer)

  if (!gene %in% rownames(Xw) || !gene %in% rownames(Xo)) stop("Gene not found: ", gene)

  df_wt <- data.frame(
    group = "Control OVX",
    dist  = seu_wt[[dist_col]][, 1],
    y     = as.numeric(Xw[gene, ])
  )
  df_oe <- data.frame(
    group = "Wnt1tg OVX",
    dist  = seu_oe[[dist_col]][, 1],
    y     = as.numeric(Xo[gene, ])
  )

  dplyr::bind_rows(df_wt, df_oe) %>%
    dplyr::filter(!is.na(dist)) %>%
    dplyr::group_by(group, dist) %>%
    dplyr::summarise(
      y_obs = mean(y, na.rm = TRUE),
      n     = sum(!is.na(y)),
      .groups = "drop"
    )
}

plot_spline_clean <- function(curves_long, gene, obs = NULL,
                              cols = c(`Control OVX` = "#00BFC4", `Wnt1tg OVX` = "#F8766D")) {
  p <- ggplot2::ggplot(
    curves_long,
    ggplot2::aes(x = dist, y = y, color = group)
  ) +
    ggplot2::geom_line(linewidth = 1.2) +
    ggplot2::scale_color_manual(values = cols) +
    scale_x_continuous(
      limits = c(0, NA),
      breaks = function(x) {

        xmax <- max(x)

        step <- if (xmax <= 12) {
          2
        } else if (xmax <= 25) {
          4
        } else if (xmax <= 50) {
          5
        } else {
          10
        }

        seq(0, floor(xmax / step) * step, by = step)
      },
      expand = c(0, 0)
    ) +
    ggplot2::theme_classic() +
    ggplot2::labs(
      title = gene,
      x = "Minimal graph distance",
      y = "Fitted expression"
    )

  if (!is.null(obs) && nrow(obs) > 0) {
    p <- p +
      ggplot2::geom_point(
        data = obs,
        position = position_dodge(width = 0.4),
        ggplot2::aes(x = dist, y = y_obs, color = group),
        inherit.aes = FALSE,
        size = 2,
        alpha = 0.9
      )
  }
  p
}



predict_gene_curves <- function(res, gene, ngrid = 200) {
  if (is.null(res$meta) || nrow(res$meta) == 0) stop("res$meta is missing/empty.")
  if (!gene %in% rownames(res$fit$coefficients)) stop("Gene not in fit: ", gene)

  dist_train <- res$meta$dist
  grid <- seq(min(dist_train), max(dist_train), length.out = ngrid)

  B <- splines::ns(
    grid,
    knots = res$spline$knots,
    Boundary.knots = res$spline$bknots
  )
  colnames(B) <- paste0("b", seq_len(ncol(B)))

  new_WT <- cbind(
    data.frame(group = factor("Control OVX", levels = levels(res$meta$group)), dist = grid),
    B
  )
  new_OE <- cbind(
    data.frame(group = factor("Wnt1tg OVX", levels = levels(res$meta$group)), dist = grid),
    B
  )

  form <- stats::as.formula(
    paste0("~ group * (", paste(colnames(B), collapse = "+"), ")")
  )

  X_WT <- stats::model.matrix(form, data = new_WT)
  X_OE <- stats::model.matrix(form, data = new_OE)

  cn <- colnames(res$design)
  X_WT <- X_WT[, cn, drop = FALSE]
  X_OE <- X_OE[, cn, drop = FALSE]

  beta <- res$fit$coefficients[gene, cn]

  wide <- data.frame(
    dist = grid,
    "Control OVX" = as.numeric(X_WT %*% beta),
    "Wnt1tg OVX" = as.numeric(X_OE %*% beta),
    check.names = FALSE
  )

  long <- tidyr::pivot_longer(
    wide,
    cols = -dist,
    names_to = "group",
    values_to = "y"
  )

  list(wide = wide, long = long)
}

spline_fit_plt <- function(res, seu_wt, seu_oe, gene,
                                    dist_col = "dist",
                                    assay = "SCT", layer = "data",
                                    ngrid = 200,
                                    cols = c(`Control OVX` = "#00BFC4", `Wnt1tg OVX` = "#F8766D")) {
  pr <- predict_gene_curves(res, gene, ngrid = ngrid)
  obs <- observed_dist_means(seu_wt, seu_oe, gene, dist_col = dist_col, assay = assay, layer = layer)
  # Restrict obs to fitted x-range (avoids points outside model range)
  obs <- obs[obs$dist >= min(pr$wide$dist) & obs$dist <= max(pr$wide$dist), , drop = FALSE]
  curves = plot_spline_clean(pr$long, gene, obs = obs, cols = cols)
  return(curves)
}

# ----------------------------
# ggvenn plot creator
# ----------------------------

make_filter_ggvenn <- function(
    sets,
    title = NULL,
    fill = c(`Control OVX` = "#00BFC4", `Wnt1tg OVX` = "#F8766D"),
    stroke_size = 0.5,
    set_name_size = 4,
    text_size = 4
) {

  if (is.null(title)) {
    title <- paste0("Detection filter (>= ", sets$perc_filter * 100, "% spots)")
  }

  ggvenn::ggvenn(
    list(`Control OVX` = sets$detected_wt, `Wnt1tg OVX` = sets$detected_oe),
    fill_color = unname(fill),
    stroke_size = stroke_size,
    set_name_size = set_name_size,
    text_size = text_size
  ) +
    ggplot2::ggtitle(title) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
}
