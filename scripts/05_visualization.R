# Library
library(here)
source(here("src","r","libraries.R"))

# Paths
## Seurat object
path_in_obj <- here("data", "preprocessed_obj.rds")
## Cavity trab dist
path_in_cavity_trab_dist <- here("data", "Cavity_trab_dist_reg.rds")
path_in_shift_cavity_trab_dist<- here("data", "gsea_shift_cavity_trab_dist.rds")
path_in_shape_cavity_trab_dist <- here("data", "gsea_shape_cavity_trab_dist.rds")
## Cavity prox dist
path_in_cavity_prox_dist <- here("data", "Cavity_prox_dist_reg.rds")
path_in_shift_cavity_prox_dist <- here("data", "gsea_shift_cavity_prox_dist.rds")
path_in_shape_cavity_prox_dist <- here("data", "gsea_shape_cavity_prox_dist.rds")
## CortBone prox dist
path_in_cortBone_prox_dist <- here("data", "CortBone_prox_dist_reg.rds")
path_in_shift_cortBone_prox_dist <- here("data", "gsea_shift_cortBone_prox_dist.rds")
path_in_shape_cortBone_prox_dist <- here("data", "gsea_shape_cortBone_prox_dist.rds")

# Functions
source(here("src","r", "functions", "spline_fit_plot.R"))
source(here("src","r", "functions", "shift_volcano.R"))
source(here("src","r", "functions", "shape_heatmap.R"))


# Data
## Cavity trab dist

trab_cap12 <- readRDS(path_in_cavity_trab_dist) #spline regression output
obj_trab_cap12 <- readRDS(path_in_obj) %>%  #normalized seurat object
  subset(spa %in% c("Cavity", "TrabBone") & Trab_Dist<max(trab_cap12$info$dist_labels_ok)) #subset to only rings
gsea_shift_trab_cap12 <- readRDS(path_in_shift_cavity_trab_dist)
gsea_shape_trab_cap12 <- readRDS(path_in_shape_cavity_trab_dist)

## Cavity prox dist
trab_bin12 <- readRDS(path_in_cavity_prox_dist) #spline regression output
obj_trab_bin12 <- readRDS(path_in_obj) %>%  #normalized seurat object
  subset(spa %in% c("Cavity", "TrabBone")) #subset to only rings
gsea_shift_trab_bin12 <- readRDS(path_in_shift_cavity_prox_dist)
gsea_shape_trab_bin12 <- readRDS(path_in_shape_cavity_prox_dist)

## CortBone prox dist
prox_cortBone <- readRDS(path_in_cortBone_prox_dist) #spline regression output
obj_cortBone <- readRDS(path_in_obj) %>%  #normalized seurat object
  subset(spa %in% c("CortBone")) #subset to only rings
gsea_shift_prox_cortBone <- readRDS(path_in_shift_cortBone_prox_dist)
gsea_shape_prox_cortBone <- readRDS(path_in_shape_cortBone_prox_dist)

# Main
## Cavity trab dist Plots ####
## Plot individual genes of interest (GOIs)
res <- trab_cap12
obj <- obj_trab_cap12
gsea_shift <- gsea_shift_trab_cap12
gsea_shape <- gsea_shape_trab_cap12
path_out <- here("results", "figures", "cavity_trab_dist")

GOIs <- c("Wnt1","Mmp9", "Col1a2")
dist_col <- "Trab_Dist"

obj_splitted <- obj %>% SplitObject(split.by='group') #split into experimental groups
obj_wt <- obj_splitted$`Control OVX`
obj_oe <- obj_splitted$`Wnt1tg OVX`


## Venn Diagrams
sets <- res$filter$sets
venn_plt <- make_filter_ggvenn(sets)
ggsave(plot=venn_plt, filename= 'venn_plot.png',
       path=path_out, dpi=300, width=6, height=3)

## Spline plots
for (gene in GOIs){
  spline_plt <- spline_fit_plt(res, obj_wt, obj_oe, gene=gene,
                                        dist_col=dist_col)
  ggsave(plot=spline_plt, filename=paste0(gene,"spline-fit.png"),
         path=path_out, dpi=300, width=6, height=3)

  spatial_plt <- SpatialFeaturePlot(obj, features=gene, pt.size.factor=2,
                                    alpha=c(0,3))
  ggsave(plot=spatial_plt, filename=paste0(gene,"spatial_plt.png"),
         path=path_out, dpi=300, width=6, height=3)
}

## Shift volcano
shift_volcano_plt <- shift_volcano(res, GOIs=GOIs, top_pct=0.1, pval_sig=0.1)
ggsave(plot=shift_volcano_plt, filename="shift_volcano_cavity_trab_dist.png",
       path=path_out, dpi=300, width=6, height=3)
## Shape heatmap
shape_heatmap_plt <- shape_heatmap(res, top_n=50, scale_rows=TRUE)
ggsave(plot=shape_heatmap_plt, filename="shape_heatmap_cavity_trab_dist.png",
       path=path_out, dpi=300, width=12, height=6)

## GSEA dotplots
if (gsea_shift@result %>% nrow() != 0){
  shift_dotplot <- dotplot(gsea_shift, split=".sign") + #split=".sign" separates positive vs negative NES.
    facet_grid(. ~ .sign) +
    theme_classic() +
    labs(title="GSEA (SHIFT): overall expression difference")
  ggsave(plot=shift_dotplot, filename="shift_dotplot_cavity_trab_dist.png",
         path=path_out, dpi=300, width=7, height=5)
}

if (gsea_shape@result %>% nrow() != 0){
  shape_dotplot <- dotplot(gsea_shape) +
    theme_classic() +
    labs(title="GSEA (SHAPE): spatial pattern changes")
  ggsave(plot=shape_dotplot, filename="shape_dotplot_cavity_trab_dist.png",
         path=path_out, dpi=300, width=7, height=5)
}


## Cavity prox dist Plots ####

res <- trab_bin12
obj <- obj_trab_bin12
gsea_shift <- gsea_shift_trab_bin12
gsea_shape <- gsea_shape_trab_bin12
dist_col <- "Prox_Dist"
path_out <- here("results", "figures", "cavity_prox_dist")

GOIs <- c("Wnt1","Mmp9", "Col1a2")

## Venn Diagrams
sets <- res$filter$sets
venn_plt <- make_filter_ggvenn(sets)
ggsave(plot=venn_plt, filename= 'venn_plot.png',
       path=path_out, dpi=300, width=6, height=3)

obj_splitted <- obj %>% SplitObject(split.by='group') #split into experimental groups
obj_wt <- obj_splitted$`Control OVX`
obj_oe <- obj_splitted$`Wnt1tg OVX`

for (gene in GOIs){
  spline_plt <- spline_fit_plt(res, obj_wt, obj_oe, gene=gene,
                               dist_col=dist_col)
  ggsave(plot=spline_plt, filename=paste0(gene,"spline-fit.png"),
         path=path_out, dpi=300, width=6, height=3)

  spatial_plt <- SpatialFeaturePlot(obj, features=gene, pt.size.factor=2,
                                    alpha=c(0,3))
  ggsave(plot=spatial_plt, filename=paste0(gene,"spatial_plt.png"),
         path=path_out, dpi=300, width=6, height=3)
}

## Shift volcano
shift_volcano_plt <- shift_volcano(res, GOIs=GOIs, top_pct=0.1, pval_sig=0.1)
ggsave(plot=shift_volcano_plt, filename="shift_volcano_cavity_prox_dist.png",
       path=path_out, dpi=300, width=6, height=3)
## Shape heatmap
shape_heatmap_plt <- shape_heatmap(res, top_n=50, scale_rows=TRUE)
ggsave(plot=shape_heatmap_plt, filename="shape_heatmap_cavity_prox_dist.png",
       path=path_out, dpi=300, width=12, height=6)

## GSEA dotplots
if (gsea_shift@result %>% nrow() != 0){
  shift_dotplot <- dotplot(gsea_shift, split=".sign") + #split=".sign" separates positive vs negative NES.
    facet_grid(. ~ .sign) +
    theme_classic() +
    labs(title="GSEA (SHIFT): overall expression difference")
  ggsave(plot=shift_dotplot, filename="shift_dotplot_cavity_prox_dist.png",
         path=path_out, dpi=300, width=7, height=5)
}

if (gsea_shape@result %>% nrow() != 0){
  shape_dotplot <- dotplot(gsea_shape) +
    theme_classic() +
    labs(title="GSEA (SHAPE): spatial pattern changes")
  ggsave(plot=shape_dotplot, filename="shape_dotplot_cavity_prox_dist.png",
         path=path_out, dpi=300, width=7, height=5)
}

## Cortical Bone prox dist Plots ####
## Plot individual genes of interest (GOIs)

res <- prox_cortBone
obj <- obj_cortBone
gsea_shift <- gsea_shift_prox_cortBone
gsea_shape <- gsea_shape_prox_cortBone
dist_col <- "Prox_Dist"
path_out <- here("results", "figures", "cortBone_prox_dist")

GOIs <- c("Wnt1","Mmp9", "Col1a2")

## Venn Diagrams
sets <- res$filter$sets
venn_plt <- make_filter_ggvenn(sets)
ggsave(plot=venn_plt, filename= 'venn_plot.png',
       path=path_out, dpi=300, width=6, height=3)

obj_splitted <- obj %>% SplitObject(split.by='group') #split into experimental groups
obj_wt <- obj_splitted$`Control OVX`
obj_oe <- obj_splitted$`Wnt1tg OVX`

for (gene in GOIs){
  spline_plt <- spline_fit_plt(res, obj_wt, obj_oe, gene=gene,
                               dist_col=dist_col)
  ggsave(plot=spline_plt, filename=paste0(gene,"_spline-fit.png"),
         path=path_out, dpi=300, width=6, height=3)

  spatial_plt <- SpatialFeaturePlot(obj, features=gene, pt.size.factor=2,
                                    alpha=c(0,3))
  ggsave(plot=spatial_plt, filename=paste0(gene,"_spatial_plt.png"),
         path=path_out, dpi=300, width=6, height=3)
}

## Shift volcano
shift_volcano_plt <- shift_volcano(res, GOIs=GOIs, top_pct=0.1, pval_sig=0.1)
ggsave(plot=shift_volcano_plt, filename="shift_volcano_cortBone_prox_dist.png",
       path=path_out, dpi=300, width=6, height=3)
## Shape heatmap
shape_heatmap_plt <- shape_heatmap(res, top_n=50, scale_rows=TRUE)
ggsave(plot=shape_heatmap_plt, filename="shape_heatmap_cortBone_prox_dist.png",
       path=path_out, dpi=300, width=12, height=6)

## GSEA dotplots
if (gsea_shift@result %>% nrow() != 0){
  shift_dotplot <- dotplot(gsea_shift, split=".sign") + #split=".sign" separates positive vs negative NES.
    facet_grid(. ~ .sign) +
    theme_classic() +
    labs(title="GSEA (SHIFT): overall expression difference")
  ggsave(plot=shift_dotplot, filename="shift_dotplot_cortBone_prox_dist.png",
         path=path_out, dpi=300, width=7, height=5)
}

if (gsea_shape@result %>% nrow() != 0){
  shape_dotplot <- dotplot(gsea_shape) +
    theme_classic() +
    labs(title="GSEA (SHAPE): spatial pattern changes")
  ggsave(plot=shape_dotplot, filename="shape_dotplot_cortBone_prox_dist.png",
         path=path_out, dpi=300, width=7, height=5)
}
