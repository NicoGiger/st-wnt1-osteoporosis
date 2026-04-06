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
obj <- AddMetaData(obj, spa) #  add obj metadata
obj$group <- factor(obj$group, levels = c("ctl OVX", "Wnt1 OVX")) #specicy factor levels

Cavity <- obj %>% subset(spa %in% c("Cavity", "TrabBone"))
CortBone <- obj %>% subset(spa=='CortBone')

# Add distance metrices

obj <- compute_min_ref_graph_distance_by_group(obj, group_col='group',
                                               reference_label="TrabBone",
                                               spa_col='spa', out_col="Trab_gDist")
obj <- compute_min_ref_graph_distance_by_group(obj, group_col='group',
                                               reference_label="Transition",
                                               spa_col='spa', out_col="Distal_gDist")

obj_sub <- subset(obj, spa %in% c("Cavity", "TrabBone", "CortBone"))

# Sample-wise SCT normalization
obj_norm <- SCT_splitwise(obj_sub)

#save normalized seurat object
saveRDS(obj_norm, here("data", "preprocessed_obj.rds"))
