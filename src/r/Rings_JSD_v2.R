library(FNN)
library(Matrix)

dd.norm <- seu
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

  list(adj = adj,
       spot_names = rownames(coords_mat),
       spacing = spacing,
       thresh = thresh)
}

build_numeric_rings <- function(dd, nrings = 5,
                                thresh_mult = 1.10,
                                k_candidates = 18){

  adj_info <- make_adj_list(dd,
                            k_candidates = k_candidates,
                            thresh_mult = thresh_mult)

  spot_names <- adj_info$spot_names
  adj <- adj_info$adj

  # initialize numeric ring column
  dd$ring <- NA_integer_

  # TrabBone = 0
  trab_spots <- rownames(dd@meta.data)[dd$spa == "TrabBone"]
  dd$ring[trab_spots] <- 0

  # frontier = TrabBone
  current_frontier <- trab_spots

  for (r in seq_len(nrings)) {

    from_idx <- match(current_frontier, spot_names)
    from_idx <- from_idx[!is.na(from_idx)]

    neigh_idx <- unique(unlist(adj[from_idx], use.names = FALSE))
    neigh_spots <- spot_names[neigh_idx]

    # only unassigned cavity spots
    candidate_spots <- rownames(dd@meta.data)[
      dd$spa == "Cavity" & is.na(dd$ring)
    ]

    new_ring_spots <- intersect(neigh_spots, candidate_spots)

    if (length(new_ring_spots) == 0) break

    dd$ring[new_ring_spots] <- r
    current_frontier <- new_ring_spots
  }

  attr(dd, "ring_spacing") <- adj_info$spacing
  attr(dd, "ring_thresh")  <- adj_info$thresh

  dd
}
dd.norm <- build_numeric_rings(dd.norm, nrings = 100)
subset(dd.norm, nrings)

table(dd.norm$ring, useNA = "ifany")

SpatialFeaturePlot(dd.norm, features = "ring") +
  scale_fill_viridis_c(
    option = "magma",
    direction = -1,
    trans = 'sqrt'
  )



expr <- GetAssayData(dd.norm, assay = "SCT", layer = "data")

rings <- sort(unique(na.omit(dd.norm$ring)))

ring_means <- sapply(rings, function(r) {
  idx <- which(dd.norm$ring == r)
  if (length(idx) == 0) {
    rep(NA_real_, nrow(expr))
  } else if (length(idx) == 1) {
    as.numeric(expr[, idx, drop = TRUE])  # single column
  } else {
    Matrix::rowMeans(expr[, idx, drop = FALSE])
  }
})

rownames(ring_means) <- rownames(expr)
colnames(ring_means) <- paste0("ring_", rings)
ring_counts <- table(dd.norm$ring)
rings_use <- as.integer(names(ring_counts)[ring_counts > 0])
\

SpatialFeaturePlot(dd.norm, features = 'ring', )
