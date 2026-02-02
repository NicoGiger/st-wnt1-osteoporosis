
suppressPackageStartupMessages({
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
})
#global ####
#palette <- DiscretePalette_scCustomize(36, palette = 'polychrome')
col13_Giger <-  c("#F0E442",  '#99DDFF', "#CC79A7", "#009E73", 
                  "#0072B2", "#B66DFF", "#E69F00", "#919492", "#D55E00",'#BBCC33',
                  "#000000", "#FFFFFF", "#ececec")

##### Paths #####
name <- "Wnt"
path.in <- paste0('~/Documents/spatial_transcriptomics/ST_runs/',name,'/')
path.out <- paste0("~/Documents/spatial_transcriptomics/ST_runs/Wnt/analysis_GIN/output/")
path.out.images <- paste0(path.out, "images/")
path.spa <- paste0("~/Documents/spatial_transcriptomics/st-wnt1-osteoporosis/metadata/spa.csv")

##### functions ####

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

k_neighbor_spots <- function(dd, k) {
  target_spots <- dd@meta.data %>% filter(spa == 'TrabBone') %>% rownames()
  Cavity_spots <- dd@meta.data %>% filter(spa == 'Cavity') %>% rownames()
  
  coords <- GetTissueCoordinates(dd)
  coords_mat <- as.matrix(coords[, c("x", "y")])
  spot_names <- rownames(coords_mat)
  
  nn <- get.knn(coords_mat, k = k)
  target_indices <- match(target_spots, spot_names)
  
  # Collect all neighbors for these spots
  neighbor_indices <- unique(as.vector(nn$nn.index[target_indices, ]))
  neighbor_spots <- spot_names[neighbor_indices]
  
  # Remove self-matches and spots defined in other tissues
  neighbor_spots <- setdiff(neighbor_spots, target_spots)
  neighbor_spots <- intersect(neighbor_spots, Cavity_spots)
  
  # Add trab neighbors to meta.data
  dd$spa[neighbor_spots] <- paste0("Trab neighbor ", k)
  
  dd$idents <- paste(dd$spa, dd$group) #add idents
  return(dd)
}
scran_DESeq_DGE <- function(obj){
  DefaultAssay(obj) <- "Spatial"
  cts <- GetAssayData(obj, layer = "counts") 
  md <- obj@meta.data
  md <- md[colnames(cts), , drop = FALSE]   # Ensure rownames(md) match colnames(cts)
  md$condition <- factor(md$group, levels = c("Wt", "WntOverExp"))
  min_prop <- 0.1
  keep_genes <- Matrix::rowSums(cts > 0) >= (ncol(cts) * min_prop) #Keep genes expressed in at least 5% of spots (adjust as desired)
  cts_f <- cts[keep_genes, ]
  message("Kept genes: ", sum(keep_genes), " / ", length(keep_genes))
  sce <- SingleCellExperiment(list(counts = cts_f)) # scran expects a SingleCellExperiment with counts
  #clusters <- quickCluster(sce) # Optional: quick clustering to improve deconvolution stability (recommended by scran)
  # scran deconvolution size factors
  sce <- computeSumFactors(sce)
  sf <- sizeFactors(sce) 
  # DESeq2 with scran size factors
  dds <- DESeqDataSetFromMatrix(
    countData = cts_f,
    colData   = md,
    design    = ~ condition
  )
  # Put scran size factors into DESeq2
  sizeFactors(dds) <- sf
  
  # LRT: test whether condition improves fit vs intercept-only model
  dds <- DESeq(dds, test="LRT", reduced=~1, minmu=1e-6, minRep=Inf)
  res <- results(dds, independentFiltering = FALSE)
  # Tidy results + thresholds
  res_df <- as.data.frame(res) %>%
    tibble::rownames_to_column("gene") %>%
    arrange(padj, pvalue)
  
  # Typical thresholds (adjust to your needs)
  cutP <- 0.05
  cutLFC <- 0.5
  
  res_df <- res_df %>%
    mutate(
      DEG = case_when(
        !is.na(padj) & padj < cutP & log2FoldChange >=  cutLFC ~ "UP",
        !is.na(padj) & padj < cutP & log2FoldChange <= -cutLFC ~ "DOWN",
        !is.na(padj) & padj > cutP & pvalue < 0.01 & log2FoldChange >= cutLFC ~ "SUGGESTIVE UP",
        !is.na(padj) & padj > cutP & pvalue < 0.01 & log2FoldChange <= -cutLFC ~ "SUGGESTIVE DOWN",
        TRUE ~ "NO"
      )
    )
  return(res_df)
}

volcano <- function(scran_deseq_res, filname = NULL, width = 7, height = 14, dpi = 300) {
  plt <- ggplot(data = scran_deseq_res, aes(x = log2FoldChange, y = -log10(pvalue), col = DEG)) +
    geom_point() +
    scale_color_manual(values = c('DOWN' = 'purple', 'NO' = 'grey', 'UP' = 'orange', 
                                  'SUGGESTIVE DOWN' = '#DB7093', 'SUGGESTIVE UP' = '#FFD700')) +
    theme(text = element_text(size = 20)) +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", col = "black") +
    geom_hline(yintercept = -log10(0.01), linetype = "dashed", col = "black") +
    geom_label_repel(
      data = scran_deseq_res %>%
        filter(DEG %in% c("UP", "DOWN")),
      aes(label = gene), label.size = 0,force = 1.5) +
    xlab("Fold Change (log2)") +
    ylab("p-value (-log10)") +
    #xlim(c(-3,3)) +
    #ylim(c(-0,50)) +
    theme_bw()
  ggsave(plot = plt, filename = filname, width = width, height = height, dpi = dpi)
  
}
##### Load Data #####

dd <- Load10X_Spatial(data.dir = path.in)
dd@images$slice1@scale.factors$spot <- 6/(dd@images$slice1@scale.factors$lowres) #correct image scale
spa <- read.csv(path.spa, stringsAsFactors = F, header = T, row.names = 1) #load metadata
dd <- AddMetaData(dd, spa) # add metadata
dd <- dd %>% subset(!(spa %in% c('out', 'Undefined')))

#Sample-wise SCT normalization
dd.norm <- SCT2_norm(dd)

#Determine Trabecular bone neighbouring spots
k_values <- c(6, 18)

dd.norm <- Reduce(
  function(obj, k) k_neighbor_spots(obj, k),
  k_values,
  init = dd.norm
)
#QC plots
SpatialDimPlot_scCustom(dd.norm, colors_use = col13_Giger, group.by = 'spa')
VlnPlot(dd.norm, split.by = 'group', group.by = 'spa', 'nCount_Spatial', 
        log = T, sort = T)
SpatialFeaturePlot(dd.norm, 'nCount_Spatial', max.cutoff = 750)


# Scran Deseq
Trab <- subset(dd.norm, subset = spa %in% c("TrabBone"))
Trabk6 <- subset(dd.norm, subset = spa %in% c("TrabBone", "Trab neighbor 6"))
Trab18 <- subset(dd.norm, subset = spa %in% c("TrabBone", "Trab neighbor 6",  "Trab neighbor 18"))
Cavity <- subset(dd.norm, subset = spa %in% c("TrabBone", "Trab neighbor 6",  "Trab neighbor 18", "Cavity"))

Trab_res <- scran_DESeq_DGE(Trab)
Trabk6_res <- scran_DESeq_DGE(Trabk6)
Trabk18_res <- scran_DESeq_DGE(Trabk18)
Cavity_res <- scran_DESeq_DGE(Cavity)

#Volcano plots
volcano(Cavity_res, filname = '~/Documents/spatial_transcriptomics/st-wnt1-osteoporosis/results/figures/volcano_full_cavity.png', 
        width = 7, height = 10)
volcano(Trab_res, filname = '~/Documents/spatial_transcriptomics/st-wnt1-osteoporosis/results/figures/volcano_trabecular-bone.png',
        width = 7, height = 10)
volcano(Trabk6_res, filname = '~/Documents/spatial_transcriptomics/st-wnt1-osteoporosis/results/figures/volcano_Trabk6.png',
        width = 7, height = 10)
volcano(Trabk18_res, filname = '~/Documents/spatial_transcriptomics/st-wnt1-osteoporosis/results/figures/volcano_Trabk18.png',
        width = 7, height = 10)
























