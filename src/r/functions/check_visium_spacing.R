check_visium_spacing <- function(seu) {
  coords <- GetTissueCoordinates(seu)
  xy <- as.matrix(coords[, c("x", "y")])

  nn <- FNN::get.knn(xy, k = 1)
  d1 <- nn$nn.dist[, 1]

  list(
    n_spots = nrow(xy),
    nn_summary = summary(d1),
    nn_sd = sd(d1),
    nn_cv = sd(d1) / mean(d1),
    nn_median = median(d1)
  )
}
