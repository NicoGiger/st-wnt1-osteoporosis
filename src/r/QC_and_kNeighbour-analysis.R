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

#global ####
#palette <- DiscretePalette_scCustomize(36, palette = 'polychrome')
col13_Giger <-  c("#F0E442", "#009E73",'#99DDFF', "#0072B2",  "#CC79A7",
                  '#BBCC33', "#B66DFF", "#E69F00", "#919492", "#D55E00",
                  "#000000", "#FFFFFF", "#ececec")

##### Paths #####
path.in <- paste0("/media/p_drive/Ulm/Wnt/")
path.spa <- paste0("~/Documents/st-wnt1-osteoporosis/metadata/spa.csv")
path.out <- paste0("~/Documents/st-wnt1-osteoporosis/results/")
path_out_figures <- paste0(path.out, "figures/")


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


scran_DESeq_DGE <- function(obj, min_prop = 0.1){
  DefaultAssay(obj) <- "Spatial"
  cts <- GetAssayData(obj, layer = "counts")
  md <- obj@meta.data
  md <- md[colnames(cts), , drop = FALSE]   # Ensure rownames(md) match colnames(cts)
  md$condition <- factor(md$group, levels = c("Wt", "WntOverExp"))
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

volcano <- function(scran_deseq_res, GOIs, filename = NULL, force = 1.5,
                    max.overlaps = 10, lim = FALSE, ylim = 100,
                    pval_sug = 0.05) {
  genes_up <- scran_deseq_res %>% filter(DEG == 'UP') %>% slice_min(order_by = pvalue, n = 5) %>% pull(gene)
  genes_down <- scran_deseq_res %>% filter(DEG == 'DOWN') %>% slice_min(order_by = pvalue, n = 5) %>% pull(gene)
  genes_updown <- unique(c(genes_up, genes_down))

  plt <- ggplot(data = scran_deseq_res, aes(x = log2FoldChange, y = -log10(pvalue), col = DEG)) +
    geom_point() +
    scale_color_manual(values = c('DOWN' = 'purple', 'NO' = 'grey', 'UP' = 'orange',
                                  'SUGGESTIVE DOWN' = '#DB7093', 'SUGGESTIVE UP' = '#FFD700')) +
    theme(text = element_text(size = 20)) +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", col = "black") +
    geom_hline(yintercept = -log10(pval_sug), linetype = "dashed", col = "black") +
    geom_label_repel(
      data = scran_deseq_res %>%
        filter(DEG != 'NO' & gene %in% unique(c(genes_updown, GOIs))),
      aes(label = gene), label.size = 0,force = force, force_pull = 0.5, max.overlaps = max.overlaps) +
    xlab("Fold Change (log2)") +
    ylab("p-value (-log10)") +

    theme_bw()
 if (lim == T ){
   plt <- plt + xlim(c(-3,3)) +
     ylim(c(-0,ylim))}
  return(plt)
}


##### Load Data #####

dd <- Load10X_Spatial(data.dir = path.in)
dd@images$slice1@scale.factors$spot <- 6/(dd@images$slice1@scale.factors$lowres) #correct image scale
spa <- read.csv(path.spa, stringsAsFactors = F, header = T, row.names = 1) #load metadata
dd <- AddMetaData(dd, spa) # add metadata
dd <- dd %>% subset(!(spa %in% c('out', 'Undefined', 'CortBone', 'GP', 'TrabBone2')) &
                      nCount_Spatial > 250)
SpatialFeaturePlot(dd, features = 'nCount_Spatial', max.cutoff = 3000)
VlnPlot(dd, features = 'nCount_Spatial', group.by = 'spa', log = T)
SpatialDimPlot(dd, group.by = 'spa')

#Sample-wise SCT normalization
dd.norm <- SCT2_norm(dd, )


#Determine Trabecular bone neighbouring spots
k_values <- c(6, 18)

dd.norm <- Reduce(
  function(obj, k) k_neighbor_spots(obj, k),
  k_values,
  init = dd.norm
)
dd.norm$spa <- factor(dd.norm$spa, levels = c('TrabBone', 'Trab neighbor 6', 'Trab neighbor 18', 'Cavity'))


#QC
SpatialFeaturePlot(dd.norm, features = 'nCount_Spatial', min.cutoff = 500, max.cutoff = 2000)
VlnPlot(dd.norm, features = 'nCount_Spatial', log = T, group.by = 'spa', split.by = 'group')
VlnPlot(dd.norm, features = 'nFeature_Spatial', log = T, group.by = 'spa', split.by = 'group')
dd.norm$CountvFeature <- dd.norm$nCount_Spatial/dd.norm$nFeature_Spatial
VlnPlot(dd.norm, features = 'CountvFeature', log = T, group.by = 'spa', split.by = 'group')

#FeatureScatter(dd.norm, feature1 = 'nFeature_Spatial', feature2 = 'nCount_Spatial', log = T, group.by = 'spa')

#clustering
dd.norm <- RunPCA(dd.norm)
dd.norm <- FindNeighbors(dd.norm)
dd.norm <- FindClusters(dd.norm, resolution = 0.7, algorithm = 4)
dd.norm <- RunTSNE(dd.norm, reduction = 'pca')
DimPlot(dd.norm)
SpatialDimPlot(dd.norm)


# Scran Deseq
Trab <- subset(dd.norm, subset = spa %in% c("TrabBone"))
Trabk6 <- subset(dd.norm, subset = spa %in% c("TrabBone", "Trab neighbor 6"))
k6 <- subset(dd.norm, subset = spa %in% c("Trab neighbor 6"))
Trabk18 <- subset(dd.norm, subset = spa %in% c("Trab neighbor 6",  "Trab neighbor 18"))
k18 <- subset(dd.norm, subset = spa %in% c("Trab neighbor 18"))
Trab_Cavity <- subset(dd.norm, subset = spa %in% c("TrabBone", "Trab neighbor 6",  "Trab neighbor 18", "Cavity"))
Cavity <- subset(dd.norm, subset = spa %in% c("Cavity"))


Trab_res <- scran_DESeq_DGE(Trab)

Trabk6_res <- scran_DESeq_DGE(Trabk6)
k6_res <- scran_DESeq_DGE(k6)

Trabk18_res <- scran_DESeq_DGE(Trabk18)
k18_res <- scran_DESeq_DGE(k18)

Trab_Cavity_res <- scran_DESeq_DGE(Trab_Cavity)
Cavity_res <- scran_DESeq_DGE(Cavity)

genes <- rownames(dd.norm)

# Gene of interests
lrp_genes <- grep('Lrp', genes, ignore.case = F, value = T)
fzd_genes <- grep('Fzd', genes, ignore.case = F, value = T)
Wnt_genes <- grep('Wnt', genes, ignore.case = F, value = T)
Col_genes <- grep('Col', genes, ignore.case = F, value = T)
other_genes <- c('Hhip', 'Dusp28')
wnt_feedback <- c(
  "Tle5","Rnf43","Znrf3","Dkk1","Dkk2",
  "Notum","Sfrp1","Sfrp2","Wif1","Nkd1","Nkd2"
)

gois <- unique(c(lrp_genes, fzd_genes, Wnt_genes, Col_genes, other_genes,
                 wnt_feedback))

#wnt_pathway <- unique(c(
#  msigdbr(species = "Mus musculus", category = "H") %>%
#    filter(gs_name == "HALLMARK_WNT_BETA_CATENIN_SIGNALING") %>%
#    pull(gene_symbol),
#  msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP:REACTOME") %>%
#    filter(gs_name == "REACTOME_SIGNALING_BY_WNT") %>%
#    pull(gene_symbol),
#  msigdbr(species = "Mus musculus", category = "C2") %>%
#    filter(gs_name == "KEGG_WNT_SIGNALING_PATHWAY") %>%
#    pull(gene_symbol))
#)

#wnt_pathway <- intersect(wnt_pathway, rownames(dd.norm))




# Volcano plots

volcano(Trab_res, filname = '~/Documents/spatial_transcriptomics/st-wnt1-osteoporosis/results/figures/volcano_trabecular-bone.png',
        width = 7, height = 10, GOIs = gois)

#volcano(Trabk6_res, filname = '~/Documents/spatial_transcriptomics/st-wnt1-osteoporosis/results/figures/volcano_Trabk6.png',
#        width = 7, height = 10, GOIs = gois)
volcano(k6_res, filname = '~/Documents/spatial_transcriptomics/st-wnt1-osteoporosis/results/figures/volcano_k6.png',
        width = 7, height = 10, GOIs = gois)

volcano(Trabk18_res, filname = '~/Documents/spatial_transcriptomics/st-wnt1-osteoporosis/results/figures/volcano_Trabk18.png',
        width = 7, height = 10, GOIs = gois)
volcano(k18_res, filname = '~/Documents/spatial_transcriptomics/st-wnt1-osteoporosis/results/figures/volcano_k18.png',
        width = 7, height = 10, GOIs = gois)

volcano(Trab_Cavity_res, filname = '~/Documents/spatial_transcriptomics/st-wnt1-osteoporosis/results/figures/volcano_full_cavity.png',
        width = 7, height = 10, GOIs = gois, max.overlaps = 50)
volcano(Cavity_res, filname = '~/Documents/spatial_transcriptomics/st-wnt1-osteoporosis/results/figures/volcano_cavity.png',
        width = 7, height = 10, GOIs = gois)


obj <- k6
k6_list <- SplitObject(k6, split.by = 'group')
DefaultAssay(obj) <- "Spatial"
cts <- GetAssayData(obj, layer = "counts")
md <- obj@meta.data
md <- md[colnames(cts), , drop = FALSE]   # Ensure rownames(md) match colnames(cts)
md$condition <- factor(md$group, levels = c("Wt", "WntOverExp"))
min_prop <- 0.1
droped_genes <- Matrix::rowSums(cts > 0) <= (ncol(cts) * min_prop)









############### Test area ############

scran_DESeq_DGE2 <- function(obj, genes_to_keep,
                             cutP=0.05, cutLFC=0.5, cutP_sug=0.05){
  DefaultAssay(obj) <- "Spatial"
  cts <- GetAssayData(obj, layer = "counts")
  md <- obj@meta.data
  md <- md[colnames(cts), , drop = FALSE]   # Ensure rownames(md) match colnames(cts)
  md$condition <- factor(md$group, levels = c("Wt", "WntOverExp"))
  cts_f <- cts[genes_to_keep, ]
  message("Kept genes: ", length(genes_to_keep), " / ", dim(cts)[[1]])
  sce <- SingleCellExperiment(list(counts = cts_f)) # scran expects a SingleCellExperiment with counts

  # scran deconvolution size factors
  clusters <- quickCluster(sce) # Optional: quick clustering to improve deconvolution stability (recommended by scran)
  sce <- scran::computeSumFactors(sce, clusters = clusters)
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

  res_df <- res_df %>%
    mutate(
      DEG = case_when(
        !is.na(padj) & padj < cutP & log2FoldChange >=  cutLFC ~ "UP",
        !is.na(padj) & padj < cutP & log2FoldChange <= -cutLFC ~ "DOWN",
        !is.na(padj) & padj > cutP & pvalue < cutP_sug & log2FoldChange >= cutLFC ~ "SUGGESTIVE UP",
        !is.na(padj) & padj > cutP & pvalue < cutP_sug & log2FoldChange <= -cutLFC ~ "SUGGESTIVE DOWN",
        TRUE ~ "NO"
      )
    )
}

GenesByGroup <- function(obj, group = 'group', min_prop){
  obj_split <- SplitObject(obj, split.by = group)

  keep_genes <- function(x, min_prop) {
    cts <- GetAssayData(x, layer = "counts")
    keep <- Matrix::rowSums(cts > 0) >= ncol(cts) * min_prop
    message("Kept genes: ", sum(keep), " / ", length(keep))
    rownames(cts)[keep]
  }
  genes_by_group <- lapply(obj_split, keep_genes, min_prop = min_prop)
  genes_intersect <- Reduce(intersect, genes_by_group)
  return(list(genes_intersect, genes_by_group))
}


Trab <- subset(dd.norm, subset = spa %in% c("TrabBone"))
Trabk6 <- subset(dd.norm, subset = spa %in% c("TrabBone", "Trab neighbor 6"))
k6 <- subset(dd.norm, subset = spa %in% c("Trab neighbor 6"))
Trabk18 <- subset(dd.norm, subset = spa %in% c("Trab neighbor 6",  "Trab neighbor 18"))
k18 <- subset(dd.norm, subset = spa %in% c("Trab neighbor 18"))
Trab_Cavity <- subset(dd.norm, subset = spa %in% c("TrabBone", "Trab neighbor 6",  "Trab neighbor 18", "Cavity"))
Cavity <- subset(dd.norm, subset = spa %in% c("Cavity"))

#make named list of seurat object to loop through
obj_list <- list(
  Trab = Trab,
  Trabk6 = Trabk6,
  Trabk18 = Trabk18,
  Trab_Cavity = Trab_Cavity,
  k6 = k6,
  k18 = k18
)
obj_list <- list(Trab = Trab,
                 Trabk6 = Trabk6,
                 Trabk18 = Trabk18,
                 Trab_Cavity = Trab_Cavity)

#Run Deseq and
for (obj_name in names(obj_list)){
  obj <- obj_list[[obj_name]]
  genes <- GenesByGroup(obj, min_prop = 0.1)

  assign(paste0(obj_name,'_genes'), genes[[2]] )

  obj_res <- scran_DESeq_DGE2(obj, genes_to_keep = genes[[1]])
  assign(paste0(obj_name,'_res'), obj_res )

  # Venn diagram
  Venn_plt <- ggVennDiagram(genes[[2]], label = "count") +
    ggplot2::theme_void()
  ggsave(plot = Venn_plt,
         filename = paste0(path_out_figures, obj_name, '_Venn.png'),
         width = 10, height = 7, dpi = 300)


  Volcano_plt <- volcano(obj_res, GOIs = gois)
  ggsave(plot = Volcano_plt, filename = paste0(path_out_figures, obj_name,
                                               '_volcano.png'),
         width = 7, height = 10, dpi = 300)
}


### ORA ######
#UP

genes <- Trabk18_genes
only_in_wnt <- setdiff(genes$WntOverExp, genes$Wt)
DEG_SUP <- obj_res %>% filter(DEG %in% c('UP', 'SUGGESTIVE UP')) %>% pull(gene)
DEG_UP <- obj_res %>% filter(DEG %in% c('UP')) %>% pull(gene)
#ONLY_SUP <- union(only_in_wnt, DEG_SUP)

genes_to_test <-DEG_UP
gene_interest_up_df <- bitr(genes_to_test, fromType = "SYMBOL",
                         toType = "ENTREZID", OrgDb = org.Mm.eg.db)
gene_interest_up <- unique(gene_interest_up_df$ENTREZID)
length(gene_interest_up)


uni_df <- bitr(Reduce(union, genes),
              fromType = "SYMBOL",
              toType = "ENTREZID", OrgDb = org.Mm.eg.db)
uni <- unique(uni_df$ENTREZID)
uni %>% length()

test <-mapIds(org.Mm.eg.db, keys = Reduce(union, genes), column = "GENETYPE",
       keytype = "SYMBOL",
       multiVals = "first")
table(test)
###For KEGG pathway, use below for "GO term" use ###

yy <-enrichKEGG(gene_interest_up, organism = 'mmu', keyType = "ncbi-geneid",
                pAdjustMethod = "BH",minGSSize = 0,
                pvalueCutoff  = 1,qvalueCutoff=1)

yy <- setReadable(yy, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")

descriptions <- yy@result %>% filter(pvalue < 0.05) %>% pull(Description)
yy@result %>% filter(p.adjust < 0.05) %>% select(Description)

barplot(yy, showCategory = descriptions, color = 'pvalue')
dotplot(yy, showCategory = descriptions, color = 'pvalue')
cnetplot(yy, showCategory = cats)
heatplot(yy, showCategory = 15)
emapplot(pairwise_termsim(yy), showCategory = 30)
dot_DOWN | dot_UP

#write.xlsx(data.frame(yy), "gene_intrest_KEGGPathway.xlsx")
Idents(Trabk6) <- Trabk6$group
t <- PrepSCTFindMarkers(Trabk6)
markers <- FindMarkers(t, 'Wt', 'WntOverExp')


ego <- enrichGO(gene_intrest,
                OrgDb = org.Mm.eg.db,keyType = "ENTREZID",
                ont = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff= 0.05,qvalueCutoff= 0.05,readable= TRUE)
head(ego)

#### Down ########

not_in_wnt <- setdiff(genes$Wt, genes$WntOverExp)
DEG_SDOWN <- obj_res %>% filter(DEG %in% c('DOWN', 'SUGGESTIVE DOWN')) %>% pull(gene)
NOT_SDOWN <- union(not_in_wnt, DEG_SDOWN)

gene_interest_down_df <- bitr(NOT_SDOWN, fromType = "SYMBOL",
                              toType = "ENTREZID", OrgDb = org.Mm.eg.db)
gene_interest_down <- unique(gene_interest_down_df$ENTREZID)

yy <-enrichKEGG(gene_interest_down, organism = 'mmu', keyType = "ncbi-geneid",
                pAdjustMethod = "BH",minGSSize = 0,
                pvalueCutoff  = 0.05,qvalueCutoff=0.1, universe = uni)
yy@result %>% filter(pvalue < 0.05) %>% pull(Description)
yy@result %>% filter(p.adjust < 0.05) %>% pull(Description)
yy@gene %>% length()

barplot(yy, color = 'p.adjust')
dotplot(yy, )
cnetplot(yy, showCategory = 10)
heatplot(yy, showCategory = 15)
emapplot(pairwise_termsim(yy), showCategory = 30)


# GSEA ####
gsea_df <- Trabk6_res %>% filter(is.finite(log2FoldChange))

gl <- gsea_df$log2FoldChange
names(gl) <- gsea_df$gene

gl <- tapply(gl, names(gl), max)

# drop array class + keep names
gl <- as.vector(gl)
names(gl) <- names(tapply(gsea_df$log2FoldChange, gsea_df$gene, max))
gl <- sort(gl, decreasing = TRUE)

library(msigdbr)
db <- msigdbr(db_species = "MM", species = "mouse")
reactome2 <- db %>%
  filter(gs_subcollection == "CP:REACTOME") %>%
  select(gs_description, gene_symbol) %>%
  distinct()
colnames(reactome2) <- c("term", "gene")
db2 <- db %>% select(gs_description, gene_symbol) %>% distinct()
gsea <- GSEA(
  geneList = gl,
  TERM2GENE = db2,
  pvalueCutoff = 0.05,
  seed = TRUE,
  nPermSimple = 10000
)



#Count the number of significant adjusted p-values.
dim(gsea)[1]

#GSEA plots

## barplot()
gsea.sorted <- gsea@result[order(gsea@result$NES, decreasing = F),]
gsea.sorted$color<-ifelse(gsea.sorted$NES<0, "orange", "purple")
par(mar=c(5,20,3,3))
barplot(gsea.sorted$NES,
        horiz = T, names=gsea.sorted$Description,
        las=2, xlab="NES",
        cex.names = 0.5,
        col=gsea.sorted$color)
abline(v=0)
## gseaplot()
gseaplot(gsea, geneSetID = gsea.sorted$Description[1])


########## Morans I #############
obj.split <- SplitObject(dd.norm, split.by = 'group')
DefaultAssay(obj.split$Wt) <- 'Spatial'
WT <- ScaleData(obj.split$Wt, layer = 'counts' )
WT <- FindSpatiallyVariableFeatures(WT, selection.method = 'moransi')
WT_morans_table <- SVFInfo(
  WT,
  assay = "Spatial",
  method = "moransi"
)
WT_morans_table['Wnt1', ]

DefaultAssay(obj.split$WntOverExp) <- 'Spatial'
WntOverExp <- ScaleData(obj.split$WntOverExp, layer = 'counts' )
WntOverExp <- FindSpatiallyVariableFeatures(WntOverExp, selection.method = 'moransi')
WntOverExp_morans_table <- SVFInfo(
  WntOverExp,
  assay = "Spatial",
  method = "moransi"
)
WntOverExp_morans_table['Wnt1', ]
Wnt1_svg <- WntOverExp_morans_table %>% arrange(MoransI_observed) %>%
  filter(MoransI_p.value < 0.05)
WT_svg <- WT_morans_table %>% arrange(MoransI_observed) %>%
  filter(MoransI_p.value < 0.05)

WntOverExp_morans_table %>% tail()
WntOverExp_morans_table$WTMoransI <- WT_morans_table$MoransI_observed
WntOverExp_morans_table$WT_p.value <- WT_morans_table$MoransI_p.value
WntOverExp_morans_table$deltaMoransI <- WntOverExp_morans_table$MoransI_observed - WntOverExp_morans_table$WTMoransI
morans_genes <- WntOverExp_morans_table %>% arrange(desc(deltaMoransI)) %>% filter(MoransI_p.value < 0.05 | WT_p.value < 0.05) %>% na.omit()
SpatialFeaturePlot(dd.norm, features = rownames(morans_genes)[1:12], ncol = 4)
WntOverExp_morans_table %>% filter(deltaMoransI > 0)
Wnt1spec_svg <- setdiff(rownames(Wnt1_svg), rownames(WT_svg))
WTspec_svg <- setdiff(rownames(WT_svg), rownames(Wnt1_svg))

SpatialFeaturePlot(dd.norm, Wnt1spec_svg[1:12], ncol = 4)
WntOverExp_morans_table[SpatiallyVariableFeatures(WntOverExp), ]


WT_spatial_genes <- WT_morans_table %>% filter(MoransI_p.value < 0.01) %>% rownames()
WntOverExp_spatial_genes <- WntOverExp_morans_table %>% filter(MoransI_p.value < 0.01) %>% rownames()

setdiff(WT_spatial_genes, WntOverExp_spatial_genes)
spatial_Wn_genes <- setdiff(WntOverExp_spatial_genes, WT_spatial_genes)
test <- intersect(only_in_wnt, spatial_Wn_genes)
test_SUP <- union(test, DEG_SUP)
