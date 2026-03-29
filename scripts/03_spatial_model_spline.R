# Library
library(here)
source(here("src","r","libraries.R"))

# Filtering parameter
perc_filter <- 0.01

# Paths
path_in <- here("data", "preprocessed_obj.rds")

# Functions
source(here("src","r", "functions", "spline_fit_regression.R"))

# Load data
obj <- readRDS(path_in)
Cavity <- subset(obj, spa %in% c("TrabBone", "Cavity"))
CortBone <- subset(obj, spa=="CortBone")

# Trab graph distance in Cavity
K <- 12 #number of distances
Cavity_K12 <- Cavity %>% subset(Trab_Dist<=K-1) # K-1: trabecular Bone is distance 0
Cavity_K12_split <- SplitObject(Cavity_K12, split.by="group")
Cavity_trab_dist_reg <- spline_limma_abs(Cavity_K12_split$`Control OVX`, Cavity_K12_split$`Wnt1tg OVX`,
                        dist_col="Trab_Dist", assay="SCT",
                        layer="data", perc_filter=perc_filter,
                        filter_assay="Spatial", filter_layer="counts",
                        max_dist=K-1, spline_df=6) # K-1: trabecular Bone is distance 0

saveRDS(Cavity_trab_dist_reg, file=here("data", "Cavity_trab_dist_reg.rds"))


# Prox_dist in Cavity
Cavity_split <- SplitObject(Cavity, split.by="group")
K <- 50
Cavity_prox_dist_reg <- spline_limma_abs(Cavity_split$`Control OVX`, Cavity_split$`Wnt1tg OVX`,
                        dist_col="Prox_Dist", assay="SCT", nbins=9,
                        layer="data", perc_filter=perc_filter,
                        filter_assay="Spatial", filter_layer="counts",
                        max_dist=K, spline_df=4) # K-1: trabecular Bone is distance 0

saveRDS(Cavity_prox_dist_reg, file=here("data", "Cavity_prox_dist_reg.rds"))

# Prox_dist in CortBone
CortBone_split <- SplitObject(CortBone, split.by="group")
K <- 50
CortBone_trab_dist_reg <- spline_limma_abs(CortBone_split$`Control OVX`, CortBone_split$`Wnt1tg OVX`,
                        dist_col="Prox_Dist", assay="SCT", nbins=9,
                        layer="data", perc_filter=perc_filter,
                        filter_assay="Spatial", filter_layer="counts",
                        max_dist=K, spline_df=4)
saveRDS(CortBone_trab_dist_reg, file=here("data", "CortBone_prox_dist_reg.rds"))



