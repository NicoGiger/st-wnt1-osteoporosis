# Library ####

library(here)
source(here("src","r","libraries.R"))

# Paths ####
path_in <- "/media/p_drive/Ulm/Wnt/"
path_spa <- here("metadata", "spa.csv")


# Functions
source(here("src", "r", "functions", "SCT_splitwise.R"))
source(here("src", "r", "functions", "compute_knn_graph_rings.R"))

# Load Data
obj <- Load10X_Spatial(data.dir=path_in)


obj@images$slice1@scale.factors$spot <- 6/(obj@images$slice1@scale.factors$lowres) #correct image scale
spa <- read.csv(path_spa, stringsAsFactors=F, header=T, row.names=1) #load metadata
obj <- AddMetaData(obj, spa) # aobj metadata
SpatialDimPlot(obj, group.by='spa')
# Add distance metrices

obj <- compute_min_ref_graph_distance_by_group(obj, group_col='group',
                                               reference_label="TrabBone",
                                               spa_col='spa', out_col="graphdist_to_TrabBone")
obj <- compute_min_ref_graph_distance_by_group(obj, group_col='group',
                                               reference_label="Transition",
                                               spa_col='spa', out_col="Prox_Dist")

obj_sub <- subset(obj, spa %in% c("Cavity", "TrabBone", "CortBone"))

# Sample-wise SCT normalization
obj_norm <- SCT_splitwise(obj_sub)

#save normalized seurat object
saveRDS(obj_norm, here("data", "preprocessed_obj.rds"))
