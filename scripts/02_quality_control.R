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

# QC plots
plt_groups_annotations <- SpatialDimPlot(obj, group.by='group',
                                         cols = c('ctl OVX'=seurat_cols[1],
                                         'Wnt1 OVX'=seurat_cols[2])) +
  theme(
    legend.text = element_text(size = 14),
    legend.title = element_blank()
  )+
  guides(
    fill = guide_legend(override.aes = list(size = 6))
  )
ggsave(plot=plt_groups_annotations,filename="groupe_annotation.png",
       path=here("results", "figures"), dpi=300, width=6, height=3)

plt_tissue_annotations <- SpatialDimPlot(obj, group.by='spa',
                                        cols=c('Cavity'=col13_Giger[1],
                                                 'CortBone'=col13_Giger[2],
                                                 'TrabBone'=col13_Giger[6])) +
  theme(
    legend.text = element_text(size = 14),
    legend.title = element_blank()
  )+
  guides(
    fill = guide_legend(override.aes = list(size = 6))
  )
ggsave(plot=plt_tissue_annotations,filename="tissue_annotation.png",
       path=here("results", "figures"), dpi=300, width=6, height=3)

counts_spatial <- SpatialFeaturePlot(obj, features='nCount_Spatial',
                                     max.cutoff=2000)
ggsave(plot=counts_spatial, filename="counts_spatial.png",
       path=here("results", "figures"), dpi=300, width=6, height=3)

counts_violine <- VlnPlot(obj, features='nCount_Spatial', log=T,
                          group.by='spa', split.by="group", cols = c('Ctll OVX'=seurat_cols[1],
                                                                    'Wnt1 OVX'=seurat_cols[2]))+
  theme(axis.text.x=element_text(angle=0, hjust=0.5))
ggsave(plot=counts_violine, filename="counts_violine.png",
       path=here("results", "figures"), dpi=300, width=6, height=3)

features_spatial <- SpatialFeaturePlot(obj, features='nFeature_Spatial', max.cutoff=2000)
ggsave(plot=features_spatial, filename="features_spatial.png",
       path=here("results", "figures"), dpi=300, width=6, height=3)

features_violine <- VlnPlot(obj, features='nFeature_Spatial', log=T,
                            group.by='spa', split.by='group', cols = c('ctl OVX'=seurat_cols[1],
                                                                       'Wnt1 OVX'=seurat_cols[2])) +
  theme(axis.text.x=element_text(angle=0, hjust=0.5))
ggsave(plot=features_violine, filename="features_violine.png",
       path=here("results", "figures"), dpi=300, width=6, height=3)


# Trabecular distances in Cavity
Cavity <- subset(obj, spa %in% c("Cavity", "TrabBone"))
spatial_Trab_gDist <- SpatialFeaturePlot(Cavity, features="Trab_gDist") +
  scale_fill_viridis_c(option='magma', direction=-1, trans='sqrt')
ggsave(plot=spatial_Trab_gDist, filename="spatial_Trab_gDist.png",
       path=here("results", "figures", "Cavity"), dpi=300, width=6, height=3)

spots_per_Trab_gDist <- plot_spots_per_gdist(Cavity, ring_col="Trab_gDist",
                                            cols = c('ctl OVX'=seurat_cols[1],
                                                     'Wnt1 OVX'=seurat_cols[2]))
ggsave(plot=spots_per_Trab_gDist, filename="spots_per_Trab_gDist.png",
       path=here("results", "figures", "Cavity"), dpi=300, width=6, height=3)


nCounts_per_Trab_gDist <- FeatureScatter(Cavity, feature1="Trab_gDist", feature2="nCount_Spatial",
                                   split.by="group", group.by="group", plot.cor=F, log=F) +
  scale_color_manual(values=seurat_cols)
ggsave(plot=nCounts_per_Trab_gDist, filename="nCounts_per_Trab_gDist.png",
       path=here("results", "figures", "Cavity"), dpi=300, width=6, height=3)

nFeature_per_Trab_gDist <- FeatureScatter(Cavity, feature1="Trab_gDist", feature2="nFeature_Spatial",
                                   split.by="group", group.by="group", plot.cor=F, log=F) +
  scale_color_manual(values=seurat_cols)
ggsave(plot=nFeature_per_Trab_gDist, filename="nFeature_per_Trab_gDist.png",
       path=here("results", "figures", "Cavity"), dpi=300, width=6, height=3)



# Prox distances in Cavity
spatial_Distal_gDist <- SpatialFeaturePlot(Cavity, features="Distal_gDist") +
  scale_fill_viridis_c(option='magma', direction=-1, trans='sqrt')
ggsave(plot=spatial_Distal_gDist, filename="spatial_Distal_gDist.png",
       path=here("results", "figures", "Cavity"), dpi=300, width=6, height=3)

spots_per_Distal_gDist <- plot_spots_per_gdist(Cavity, ring_col='Distal_gDist',
                                            cols = c('ctl OVX'=seurat_cols[1],
                                                     'Wnt1 OVX'=seurat_cols[2]))
ggsave(plot=spots_per_Distal_gDist, filename="spots_per_Distal_gDist.png",
       path=here("results", "figures", "Cavity"), dpi=300, width=6, height=3)

nCounts_per_Distal_gDist <- FeatureScatter(Cavity, feature1="Distal_gDist", feature2="nCount_Spatial",
                                   split.by="group", group.by="group", plot.cor=F, log=F) +
  scale_color_manual(values=seurat_cols)
ggsave(plot=nCounts_per_Distal_gDist, filename="nCounts_per_Distal_gDist.png",
       path=here("results", "figures", "Cavity"), dpi=300, width=6, height=3)

nFeature_per_Distal_gDist <- FeatureScatter(Cavity, feature1="Distal_gDist", feature2="nFeature_Spatial",
                                    split.by="group", group.by="group", plot.cor=F, log=F) +
  scale_color_manual(values=seurat_cols)
ggsave(plot=nFeature_per_Distal_gDist, filename="nFeature_per_Distal_gDist.png",
       path=here("results", "figures", "Cavity"), dpi=300, width=6, height=3)


# Prox distances in Cortex
CortBone <- subset(obj, spa == "CortBone")
spatial_Distal_gDist <- SpatialFeaturePlot(CortBone, features="Distal_gDist") +
  scale_fill_viridis_c(option='magma', direction=-1, trans='sqrt')
ggsave(plot=spatial_Distal_gDist, filename="spatial_Distal_gDist.png",
       path=here("results", "figures", "CortBone"), dpi=300, width=6, height=3)

spots_per_Distal_gDist <- plot_spots_per_gdist(CortBone, ring_col='Distal_gDist',
                                            cols = c('ctl OVX'=seurat_cols[1],
                                                     'Wnt1 OVX'=seurat_cols[2]))
ggsave(plot=spots_per_Distal_gDist, filename="spots_per_Distal_gDist.png",
       path=here("results", "figures", "CortBone"), dpi=300, width=6, height=3)

nCounts_per_Distal_gDist <- FeatureScatter(CortBone, feature1="Distal_gDist", feature2="nCount_Spatial",
                                        split.by="group", group.by="group", plot.cor=F, log=F) +
  scale_color_manual(values=seurat_cols)
ggsave(plot=nCounts_per_Distal_gDist, filename="nCounts_per_Distal_gDist.png",
       path=here("results", "figures", "CortBone"), dpi=300, width=6, height=3)

nFeature_per_Distal_gDist <- FeatureScatter(CortBone, feature1="Distal_gDist", feature2="nFeature_Spatial",
                                         split.by="group", group.by="group", plot.cor=F, log=F) +
  scale_color_manual(values=seurat_cols)
ggsave(plot=nFeature_per_Distal_gDist, filename="nFeature_per_Distal_gDist.png",
       path=here("results", "figures", "CortBone"), dpi=300, width=6, height=3)

