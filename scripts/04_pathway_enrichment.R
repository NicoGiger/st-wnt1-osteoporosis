# Global
set.seed(123)

# Library
library(here)
source(here("src","r","libraries.R"))

# Paths

##Cavity trab dist
path_in_cavity_Trab_gDist <- here("data", "Cavity_Trab_gDist_reg.rds")
##Cavity prox dist
path_in_cavity_Distal_gDist <- here("data", "Cavity_Distal_gDist_reg.rds")
##CortBone prox dist
path_in_cortBone_Distal_gDist <- here("data", "CortBone_Distal_gDist_reg.rds")


# Load data
cavity_Trab_gDist <- readRDS(path_in_cavity_Trab_gDist)
cavity_Distal_gDist <- readRDS(path_in_cavity_Distal_gDist)
cortBone_Distal_gDist <- readRDS(path_in_cortBone_Distal_gDist)

# Get Hallmark (H) database
db_hallmark <- msigdbr(species="Mus musculus", collection="H") %>%
  transmute(term=gs_name, gene=gene_symbol) %>%
  distinct()


# Main

## Cavity trab dist ######
## SHIFT ranking (signed)
res <- cavity_Trab_gDist
gl_shift <- res$shift$t
names(gl_shift) <- rownames(res$shift)
gl_shift <- gl_shift[is.finite(gl_shift)]
gl_shift <- sort(gl_shift, decreasing=TRUE)

## SHAPE ranking (unsigned)
gl_shape_unsigned <- res$shape$F
names(gl_shape_unsigned) <- rownames(res$shape)
gl_shape_unsigned <- gl_shape_unsigned[is.finite(gl_shape_unsigned)]
gl_shape_unsigned <- sort(gl_shape_unsigned, decreasing=TRUE)

## Run gsea
TERM2GENE_use <- db_hallmark

gsea_shape <- GSEA(
  geneList=gl_shape_unsigned,
  TERM2GENE=TERM2GENE_use,
  pvalueCutoff=0.05,
  seed=TRUE,
  nPermSimple=10000,
  scoreType="pos"
)

gsea_shift <- GSEA(
  geneList=gl_shift,
  TERM2GENE=TERM2GENE_use,
  pvalueCutoff=0.05,
  seed=TRUE,
  nPermSimple=10000)

## save data
saveRDS(gsea_shift, file=here("data", "gsea_shift_cavity_Trab_gDist.rds"))
write.csv2(gsea_shift@result,
           file=here("results", "gsea_shift_cavity_Trab_gDist.csv"), row.names = TRUE)

saveRDS(gsea_shape, file=here("data", "gsea_shape_cavity_Trab_gDist.rds"))
write.csv2(gsea_shape@result,
           file=here("results", "gsea_shape_cavity_Trab_gDist.csv"), row.names = TRUE)


## Cavity prox dist ######

## SHIFT ranking (signed)
res <- cavity_Distal_gDist
gl_shift <- res$shift$t
names(gl_shift) <- rownames(res$shift)
gl_shift <- gl_shift[is.finite(gl_shift)]
gl_shift <- sort(gl_shift, decreasing=TRUE)

## SHAPE ranking (unsigned)
gl_shape_unsigned <- res$shape$F
names(gl_shape_unsigned) <- rownames(res$shape)
gl_shape_unsigned <- gl_shape_unsigned[is.finite(gl_shape_unsigned)]
gl_shape_unsigned <- sort(gl_shape_unsigned, decreasing=TRUE)

## Run gsea
TERM2GENE_use <- db_hallmark

gsea_shape <- GSEA(
  geneList=gl_shape_unsigned,
  TERM2GENE=TERM2GENE_use,
  pvalueCutoff=0.05,
  seed=TRUE,
  nPermSimple=10000,
  scoreType="pos"
)

gsea_shift <- GSEA(
  geneList=gl_shift,
  TERM2GENE=TERM2GENE_use,
  pvalueCutoff=0.05,
  seed=TRUE,
  nPermSimple=10000)

## save data
saveRDS(gsea_shift, file=here("data", "gsea_shift_cavity_Distal_gDist.rds"))
write.csv2(gsea_shift@result,
           file=here("results", "gsea_shift_cavity_Distal_gDist.csv"), row.names = TRUE)

saveRDS(gsea_shape, file=here("data", "gsea_shape_cavity_Distal_gDist.rds"))
write.csv2(gsea_shape@result,
           file=here("results", "gsea_shape_cavity_Distal_gDist.csv"), row.names = TRUE)

## Cortex prox dist ####

## SHIFT ranking (signed)
res <- cortBone_Distal_gDist
gl_shift <- res$shift$t
names(gl_shift) <- rownames(res$shift)
gl_shift <- gl_shift[is.finite(gl_shift)]
gl_shift <- sort(gl_shift, decreasing=TRUE)

## SHAPE ranking (unsigned)
gl_shape_unsigned <- res$shape$F
names(gl_shape_unsigned) <- rownames(res$shape)
gl_shape_unsigned <- gl_shape_unsigned[is.finite(gl_shape_unsigned)]
gl_shape_unsigned <- sort(gl_shape_unsigned, decreasing=TRUE)

## Run gsea
TERM2GENE_use <- db_hallmark

gsea_shape <- GSEA(
  geneList=gl_shape_unsigned,
  TERM2GENE=TERM2GENE_use,
  pvalueCutoff=0.05,
  seed=TRUE,
  nPermSimple=10000,
  scoreType="pos"
)

gsea_shift <- GSEA(
  geneList=gl_shift,
  TERM2GENE=TERM2GENE_use,
  pvalueCutoff=0.05,
  seed=TRUE,
  nPermSimple=10000)

## save data
saveRDS(gsea_shift, file=here("data", "gsea_shift_cortBone_Distal_gDist.rds"))
write.csv2(gsea_shift@result,
           file=here("results", "gsea_shift_cortBone_Distal_gDist.csv"), row.names = TRUE)

saveRDS(gsea_shape, file=here("data", "gsea_shape_cortBone_Distal_gDist.rds"))
write.csv2(gsea_shape@result,
           file=here("results", "gsea_shape_cortBone_Distal_gDist.csv"), row.names = TRUE)
