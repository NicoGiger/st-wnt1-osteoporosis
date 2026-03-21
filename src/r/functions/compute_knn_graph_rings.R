make_adj_list <- function(dd, k_candidates = 18, thresh_mult = 1.10){
  coords <- GetTissueCoordinates(dd)
  coords_mat <- as.matrix(coords[, c("x","y")])
  n <- nrow(coords_mat)

  nn2 <- get.knn(coords_mat, k = 2)
  spacing <- median(nn2$nn.dist[, 2], na.rm = TRUE)
  thresh <- spacing * thresh_mult

  nn <- get.knn(coords_mat, k = k_candidates)

  adj <- vector("list", n)
  for (i in seq_len(n)) {
    keep <- nn$nn.dist[i, ] <= thresh
    adj[[i]] <- nn$nn.index[i, keep]
  }

  list(adj = adj, spot_names = rownames(coords_mat), spacing = spacing, thresh = thresh)
}

compute_min_ref_graph_distance_by_group <- function(
    dd,
    group_col = "group",
    reference_label = "TrabBone",
    spa_col = "spa",
    out_col = "dist_to_ref_graph"
) {
  meta <- dd@meta.data

  if (!group_col %in% colnames(meta)) {
    stop("Column '", group_col, "' not found in dd@meta.data")
  }

  if (!spa_col %in% colnames(meta)) {
    stop("Column '", spa_col, "' not found in dd@meta.data")
  }

  adj_info <- make_adj_list(dd)
  spot_names <- adj_info$spot_names
  adj <- adj_info$adj

  if (length(spot_names) != length(adj)) {
    stop("Length of 'spot_names' and 'adj' must match.")
  }

  spot_to_idx <- stats::setNames(seq_along(spot_names), spot_names)

  out <- rep(NA_integer_, nrow(meta))
  names(out) <- rownames(meta)

  groups <- unique(meta[[group_col]])
  groups <- groups[!is.na(groups)]

  for (grp in groups) {
    grp_spots <- rownames(meta)[meta[[group_col]] == grp]
    grp_spots <- intersect(grp_spots, spot_names)

    if (length(grp_spots) == 0) {
      next
    }

    ref_spots <- rownames(meta)[
      meta[[group_col]] == grp & meta[[spa_col]] == reference_label
    ]
    ref_spots <- intersect(ref_spots, spot_names)

    if (length(ref_spots) == 0) {
      warning(
        "No reference spots found for group '", grp,
        "' and reference_label = '", reference_label, "'."
      )
      next
    }

    grp_idx <- unname(spot_to_idx[grp_spots])
    ref_idx <- unname(spot_to_idx[ref_spots])

    dist_vec <- rep(NA_integer_, length(spot_names))
    dist_vec[ref_idx] <- 0L

    frontier <- ref_idx
    visited_in_group <- rep(FALSE, length(spot_names))
    visited_in_group[grp_idx] <- TRUE

    current_dist <- 0L

    while (length(frontier) > 0) {
      neigh <- unique(unlist(adj[frontier], use.names = FALSE))

      neigh <- neigh[visited_in_group[neigh]]
      new_nodes <- neigh[is.na(dist_vec[neigh])]

      if (length(new_nodes) == 0) {
        break
      }

      current_dist <- current_dist + 1L
      dist_vec[new_nodes] <- current_dist
      frontier <- new_nodes
    }

    out[grp_spots] <- dist_vec[grp_idx]
  }

  dd[[out_col]] <- out
  dd
}
