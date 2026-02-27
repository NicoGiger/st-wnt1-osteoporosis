library(Seurat)
library(Matrix)
library(dplyr)
library(SingleCellExperiment)
library(scran)
library(DESeq2)
library(scCustomize)
library(FNN)
library(ggplot2)
library(ggrepel)
library(ggVennDiagram)
library(AnnotationDbi)
library(org.Mm.eg.db)
library(enrichplot)
library(clusterProfiler)
library(plotly)


# Global ####

col13_Giger <-  c("#F0E442", "#009E73",'#99DDFF', "#0072B2",  "#CC79A7",
                  '#BBCC33', "#B66DFF", "#E69F00", "#919492", "#D55E00",
                  "#000000", "#FFFFFF", "#ececec")

# Paths ####
path.in <- paste0("/media/p_drive/Ulm/Wnt/")
path.spa <- paste0("~/Documents/st-wnt1-osteoporosis/metadata/spa.csv")
path.out <- paste0("~/Documents/st-wnt1-osteoporosis/results/")
path_out_figures <- paste0(path.out, "figures/")

# Functions ####
SCT2_norm <- function(data, split.by = "group", assay ="Spatial",
                      slot = "counts"){
  data.list <- SplitObject(data, split.by = split.by)
  data.norm.list <- list()
  dds.list <- list()
  var_feat.list <- list()

  for (sample in data.list){
    sample <- SCTransform(sample, assay = assay, new.assay.name = "SCT",
                          return.only.var.genes = FALSE, variable.features.n = 3000)
    #vars.to.regress = "Wnt1") # If Wnt1 was orders of magnitudes higher
    data.norm.list <- append(data.norm.list, sample)
    var_feat.sample <- list(VariableFeatures(sample, assay = "SCT"))
    var_feat.list <- c(var_feat.list, var_feat.sample)
  }

  # Proper intersection across all groups
  var_feat <- Reduce(intersect, var_feat.list)

  # Merge ALL groups (not only the first two)
  merged <- Reduce(function(x, y) merge(x, y), data.norm.list)

  DefaultAssay(merged) <- "SCT"
  VariableFeatures(merged) <- var_feat
  merged <- PrepSCTFindMarkers(merged)

  # Keep original image (minimal change from your approach)
  merged@images$slice1 <- data@images$slice1
  merged@images$slice1.2 <- NULL

  return(merged)
}
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

build_numeric_rings <- function(dd, nrings = 50, thresh_mult = 1.10, k_candidates = 18,
                                trab_label = "TrabBone", cavity_label = "Cavity") {

  adj_info <- make_adj_list(dd, k_candidates = k_candidates, thresh_mult = thresh_mult)
  spot_names <- adj_info$spot_names
  adj <- adj_info$adj

  # Capture masks ONCE (so it works even if you later change dd$spa)
  trab_spots <- rownames(dd@meta.data)[dd$spa == trab_label]
  cavity_spots <- rownames(dd@meta.data)[dd$spa == cavity_label]

  dd$ring <- NA_integer_
  dd$ring[trab_spots] <- 0L

  frontier <- trab_spots

  for (r in seq_len(nrings)) {
    from_idx <- match(frontier, spot_names)
    from_idx <- from_idx[!is.na(from_idx)]
    if (length(from_idx) == 0) break

    neigh_idx <- unique(unlist(adj[from_idx], use.names = FALSE))
    neigh_spots <- spot_names[neigh_idx]

    # Only cavity spots, and only those not yet assigned
    new_frontier <- neigh_spots[neigh_spots %in% cavity_spots & is.na(dd$ring[neigh_spots])]

    if (length(new_frontier) == 0) break

    dd$ring[new_frontier] <- r
    frontier <- new_frontier
  }

  attr(dd, "ring_spacing") <- adj_info$spacing
  attr(dd, "ring_thresh")  <- adj_info$thresh
  dd
}
##### Load Data #####

dd <- Load10X_Spatial(data.dir = path.in)
dd@images$slice1@scale.factors$spot <- 6/(dd@images$slice1@scale.factors$lowres) #correct image scale
spa <- read.csv(path.spa, stringsAsFactors = F, header = T, row.names = 1) #load metadata
dd <- AddMetaData(dd, spa) # add metadata
dd <- dd %>% subset(!(spa %in% c('out', 'Undefined', 'CortBone', 'GP', 'TrabBone2')) &
                      nCount_Spatial > 250)
SpatialFeaturePlot(dd, features = 'nFeature_Spatial', max.cutoff = 3000)

dd.norm <- SCT2_norm(dd) #Sample-wise SCT normalization

# Rings ####
nrings <- 15
dd.norm <- build_numeric_rings(dd.norm, nrings = nrings)
dd.norm <- subset(dd.norm, ring < nrings)

#SpatialFeaturePlot(dd.norm, features = 'ring')
SpatialFeaturePlot(dd.norm, features = "ring", crop = F) +
  scale_fill_viridis_c(
    option = "magma",
    trans = "sqrt",
    direction = -1
  )
