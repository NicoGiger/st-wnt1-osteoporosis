shift_volcano <- function(res, GOIs = character(), filename = NULL,
                          top_pct = 0.05,
                          pval_sig = 0.05,
                          use_adj_for_y = FALSE,
                          n_label = 10,
                          lim = FALSE,
                          xlim_abs = 5,
                          ylim = 20) {

  tt <- res$shift
  if (is.null(tt) || nrow(tt) == 0) stop("res$shift is missing or empty.")
  tt <- tibble::as_tibble(tt, rownames = "gene")

  pcol <- if (use_adj_for_y) "adj.P.Val" else "P.Value"

  tt <- tt |>
    dplyr::mutate(
      pvalue = .data[[pcol]],
      y = -log10(pmax(pvalue, .Machine$double.xmin))
    )

  thr <- stats::quantile(abs(tt$logFC), probs = 1 - top_pct, na.rm = TRUE)

  tt <- tt |>
    dplyr::mutate(
      DEG = dplyr::case_when(
        logFC >= thr & pvalue < pval_sig ~ "UP",
        logFC <= -thr & pvalue < pval_sig ~ "DOWN",
        TRUE ~ "NO"
      )
    )

  genes_up <- tt |>
    dplyr::filter(DEG == "UP") |>
    dplyr::slice_min(pvalue, n = n_label) |>
    dplyr::pull(gene)

  genes_down <- tt |>
    dplyr::filter(DEG == "DOWN") |>
    dplyr::slice_min(pvalue, n = n_label) |>
    dplyr::pull(gene)

  genes_to_label <- unique(c(genes_up, genes_down, GOIs))

  tt <- tt |>
    dplyr::mutate(
      plot_group = dplyr::case_when(
        gene %in% GOIs & DEG == "UP" ~ "UP",
        gene %in% GOIs & DEG == "DOWN" ~ "DOWN",
        gene %in% genes_up ~ "UP",
        gene %in% genes_down ~ "DOWN",
        DEG %in% c("UP", "DOWN") ~ "SIG",
        TRUE ~ "NO"
      )
    )

  plt <- ggplot(tt, aes(logFC, y, color = plot_group)) +
    geom_point(size = 1.3) +
    geom_vline(xintercept = c(-thr, thr), linetype = "dashed", color = "grey40") +
    geom_hline(yintercept = -log10(pval_sig), linetype = "dashed", color = "grey40") +

    ggrepel::geom_text_repel(
      data = dplyr::filter(tt, gene %in% genes_to_label),
      aes(label = gene),
      size = 3,
      segment.color = "grey50",
      max.overlaps = Inf
    ) +

    scale_color_manual(values = c(
      UP = "#D9A300",      # orange
      DOWN = "#8A2BE2",    # purple
      SIG = "black",
      NO = "grey75"
    )) +

    xlab("Effect size (SCT residual units)") +
    ylab(paste0("p-value (-log10)")) +

    theme_classic() +
    theme(
      panel.grid.major.x = element_line(color = "grey80", linetype = "dotted"),
      panel.grid.major.y = element_blank(),
      legend.position = "none"
    )

  if (isTRUE(lim)) {
    plt <- plt + coord_cartesian(xlim = c(-xlim_abs, xlim_abs), ylim = c(0, ylim))
  }

  plt
}
