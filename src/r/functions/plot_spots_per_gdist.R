plot_spots_per_gdist <- function(
    seu,
    ring_col = "graphdist_to_TrabBone",
    group_col = "group",
    cols = NULL,
    point_size = 2.2
) {

  meta <- seu@meta.data

  obs <- meta %>%
    dplyr::select(
      dist = dplyr::all_of(ring_col),
      group = dplyr::all_of(group_col)
    ) %>%
    dplyr::mutate(
      dist = as.numeric(dist),
      group = as.character(group)
    )

  obs <- obs %>%
    dplyr::count(group, dist, name = "n") %>%
    dplyr::arrange(group, dist) %>% na.omit()

  ymax <- max(obs$n, na.rm = TRUE)

  p <- ggplot(obs, aes(x = dist, y = n, color = group)) +
    geom_line(linewidth = 0.9) +
    geom_point(size = point_size) +
    scale_y_continuous(
      breaks = pretty(c(0, ymax), n = 8),
      minor_breaks = seq(0, ymax, by = 1)
    ) +
    theme_classic() +
    labs(
      x = "Minimal graph distance",
      y = "Nr of Spots",
      color = "Group"
    )

  if (!is.null(cols)) {
    p <- p + scale_color_manual(values = cols)
  }

  return(p)
}
