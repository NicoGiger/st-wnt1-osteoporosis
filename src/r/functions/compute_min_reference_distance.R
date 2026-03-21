compute_min_ref_graph_distance <- function(
    dd,
    reference_label = "TrabBone",
    spa_col = "spa",
    out_col = "dist_to_ref_graph"
) {

  meta <- dd@meta.data

  if (!spa_col %in% colnames(meta)) {
    stop("Column '", spa_col, "' not found in dd@meta.data")
  }

  ref_spots <- rownames(meta)[meta[[spa_col]] == reference_label]

  if (length(ref_spots) == 0) {
    stop("No spots found for reference_label = ", reference_label)
  }

  # build adjacency graph
  adj_info <- make_adj_list(dd)
  spot_names <- adj_info$spot_names
  adj <- adj_info$adj

  n <- length(spot_names)

  dist_vec <- rep(NA_integer_, n)
  names(dist_vec) <- spot_names

  # initialize BFS
  ref_idx <- match(ref_spots, spot_names)
  ref_idx <- ref_idx[!is.na(ref_idx)]

  if (length(ref_idx) == 0) {
    stop("Reference spots not found in adjacency graph")
  }

  dist_vec[ref_idx] <- 0L
  frontier <- ref_idx

  repeat {

    neigh <- unique(unlist(adj[frontier], use.names = FALSE))

    new_nodes <- neigh[is.na(dist_vec[neigh])]

    if (length(new_nodes) == 0) break

    dist_vec[new_nodes] <- dist_vec[frontier[1]] + 1L

    frontier <- new_nodes
  }

  # write back to metadata
  out <- rep(NA_integer_, nrow(meta))
  names(out) <- rownames(meta)

  common <- intersect(spot_names, rownames(meta))
  out[common] <- dist_vec[common]

  dd[[out_col]] <- out

  dd
}
