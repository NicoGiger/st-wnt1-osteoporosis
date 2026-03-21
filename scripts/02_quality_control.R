# Library
library(here)
source(here("src","r","libraries.R"))

# Custom color map
source(here("src","r", "color_paletts.R"))

# Paths
path_in <- here("data", "preprocessed_obj.rds")

# Functions
source(here("src","r", "functions", "plot_spots_per_gdist.R"))

# Load data
obj <- readRDS(path_in)
Cavity <- obj %>% subset(spa %in% c("Cavity", "TrabBone"))
CortBone <- obj %>% subset(spa=='CortBone')
# QC plots
plt_groups_annotations <- SpatialDimPlot(obj, group.by='group')
ggsave(plot=plt_groups_annotations,filename="groupe_annotation.png",
       path=here("results", "figures"), dpi=300, width=6, height=3)

plt_tissue_annotations <- SpatialDimPlot(obj, group.by='spa',
                                        cols=c('Cavity'=col13_Giger[1],
                                                 'CortBone'=col13_Giger[2],
                                                 'TrabBone'=col13_Giger[6]))
ggsave(plot=plt_tissue_annotations,filename="tissue_annotation.png",
       path=here("results", "figures"), dpi=300, width=6, height=3)

counts_spatial <- SpatialFeaturePlot(obj, features='nCount_Spatial',
                                     max.cutoff=2000)
ggsave(plot=counts_spatial, filename="counts_spatial.png",
       path=here("results", "figures"), dpi=300, width=6, height=3)

counts_violine <- VlnPlot(obj, features='nCount_Spatial', log=T,
                          group.by='spa', split.by="group", ) +
  theme(axis.text.x=element_text(angle=0, hjust=0.5))
ggsave(plot=counts_violine, filename="counts_violine.png",
       path=here("results", "figures"), dpi=300, width=6, height=3)

features_spatial <- SpatialFeaturePlot(obj, features='nFeature_Spatial', max.cutoff=2000)
ggsave(plot=features_spatial, filename="features_spatial.png",
       path=here("results", "figures"), dpi=300, width=6, height=3)

features_violine <- VlnPlot(obj, features='nFeature_Spatial', log=T,
                            group.by='spa', split.by='group') +
  theme(axis.text.x=element_text(angle=0, hjust=0.5))
ggsave(plot=features_violine, filename="features_violine.png",
       path=here("results", "figures"), dpi=300, width=6, height=3)


# Trabecular distances in Cavity
spatial_trab_dist <- SpatialFeaturePlot(Cavity, features="graphdist_to_TrabBone") +
  scale_fill_viridis_c(option='magma', direction=-1, trans='sqrt')
ggsave(plot=spatial_trab_dist, filename="spatial_trab_dist.png",
       path=here("results", "figures", "Cavity"), dpi=300, width=6, height=3)

spots_per_trab_dist <- plot_spots_per_gdist(Cavity, ring_col="graphdist_to_TrabBone")
ggsave(plot=spots_per_trab_dist, filename="spots_per_trab_dist.png",
       path=here("results", "figures", "Cavity"), dpi=300, width=6, height=3)

nCounts_per_trab_dist <- FeatureScatter(Cavity, feature1="graphdist_to_TrabBone", feature2="nCount_Spatial",
                                   split.by="group", group.by="group", plot.cor=F, log=F) +
  scale_color_manual(values=seurat_cols)
ggsave(plot=nCounts_per_trab_dist, filename="nCounts_per_trab_dist.png",
       path=here("results", "figures", "Cavity"), dpi=300, width=6, height=3)

nFeature_per_trab_dist <- FeatureScatter(Cavity, feature1="graphdist_to_TrabBone", feature2="nFeature_Spatial",
                                   split.by="group", group.by="group", plot.cor=F, log=F) +
  scale_color_manual(values=seurat_cols)
ggsave(plot=nFeature_per_trab_dist, filename="nFeature_per_trab_dist.png",
       path=here("results", "figures", "Cavity"), dpi=300, width=6, height=3)



# Prox distances in Cavity
spatial_prox_dist <- SpatialFeaturePlot(Cavity, features="Prox_Dist") +
  scale_fill_viridis_c(option='magma', direction=-1, trans='sqrt')
ggsave(plot=spatial_prox_dist, filename="spatial_prox_dist.png",
       path=here("results", "figures", "Cavity"), dpi=300, width=6, height=3)

spots_per_prox_dist <- plot_spots_per_gdist(Cavity, ring_col='Prox_Dist')
ggsave(plot=spots_per_prox_dist, filename="spots_per_prox_dist.png",
       path=here("results", "figures", "Cavity"), dpi=300, width=6, height=3)

nCounts_per_prox_dist <- FeatureScatter(Cavity, feature1="Prox_Dist", feature2="nCount_Spatial",
                                   split.by="group", group.by="group", plot.cor=F, log=F) +
  scale_color_manual(values=seurat_cols)
ggsave(plot=nCounts_per_prox_dist, filename="nCounts_per_prox_dist.png",
       path=here("results", "figures", "Cavity"), dpi=300, width=6, height=3)

nFeature_per_prox_dist <- FeatureScatter(Cavity, feature1="Prox_Dist", feature2="nFeature_Spatial",
                                    split.by="group", group.by="group", plot.cor=F, log=F) +
  scale_color_manual(values=seurat_cols)
ggsave(plot=nFeature_per_prox_dist, filename="nFeature_per_prox_dist.png",
       path=here("results", "figures", "Cavity"), dpi=300, width=6, height=3)


# Prox distances in Cortex
spatial_prox_dist <- SpatialFeaturePlot(CortBone, features="Prox_Dist") +
  scale_fill_viridis_c(option='magma', direction=-1, trans='sqrt')
ggsave(plot=spatial_prox_dist, filename="spatial_prox_dist.png",
       path=here("results", "figures", "CortBone"), dpi=300, width=6, height=3)

spots_per_prox_dist <- plot_spots_per_gdist(CortBone, ring_col='Prox_Dist')
ggsave(plot=spots_per_prox_dist, filename="spots_per_prox_dist.png",
       path=here("results", "figures", "CortBone"), dpi=300, width=6, height=3)

nCounts_per_prox_dist <- FeatureScatter(CortBone, feature1="Prox_Dist", feature2="nCount_Spatial",
                                        split.by="group", group.by="group", plot.cor=F, log=F) +
  scale_color_manual(values=seurat_cols)
ggsave(plot=nCounts_per_prox_dist, filename="nCounts_per_prox_dist.png",
       path=here("results", "figures", "CortBone"), dpi=300, width=6, height=3)

nFeature_per_prox_dist <- FeatureScatter(CortBone, feature1="Prox_Dist", feature2="nFeature_Spatial",
                                         split.by="group", group.by="group", plot.cor=F, log=F) +
  scale_color_manual(values=seurat_cols)
ggsave(plot=nFeature_per_prox_dist, filename="nFeature_per_prox_dist.png",
       path=here("results", "figures", "CortBone"), dpi=300, width=6, height=3)

