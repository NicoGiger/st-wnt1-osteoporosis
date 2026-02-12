##SIB enrichment analysis course - Extra excersie for ECTS credits

#path
setwd("~/Documents/SIB_courses/EnrichmentAnalysis/")

#library
library(clusterProfiler)
library(enrichplot)
library(pathview)
library(org.Hs.eg.db)
library(ggplot2)
library(ggrepel)
library(msigdbr)
library(tidyverse) # for bonus code/dplyr/pipe

#body

#import data 
set.seed(1234)


# Import DE table:
NK_vs_Th<-read.csv("/home/ngiger/Documents/SIB_courses/EnrichmentAnalysis/data/NK_vs_Th_diff_gene_exercise_1.csv",
                   header = T)
gl<-#########Pathway_analysis.R##############
library("AnnotationDbi")
library("org.Mm.eg.db")
library("clusterProfiler")
library(gage)
library(org.Mm.eg.db)
library(genefilter)
library(annotate)
library("GO.db")
library("GOstats")
#library(KEGG.db)
library(dplyr)
library(rWikiPathways)
############### use symbol to gene ID #################### here gene name in the vetor ##########
#####Marker gene after extracting from ML model ##########

gene_intrest<-mapIds(org.Mm.eg.db,
       keys = c("Brac2", "Tp53"),
       column = "ENTREZID",
       keytype = "SYMBOL",
       multiVals = "first")

###For KEGG pathway, use below for "GO term" use ###

yy <-enrichKEGG(gene_intrest, 'mmu',keyType = "kegg",
                 pAdjustMethod = "BH",minGSSize=20,
                 pvalueCutoff  = 0.05,qvalueCutoff=0.1)

#barplot(yy1, drop = TRUE)
#dotplot(yy1)
write.xlsx(data.frame(yy), "gene_intrest_KEGGPathway.xlsx")
##########Go Term############################
##########g############################
###change "ont" of "BP", "MF", and "CC" subontologies, or "ALL" for all three.

ego <- enrichGO(gene_intrest,
                OrgDb = org.Mm.eg.db,keyType = "ENTREZID",
                ont = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff= 0.01,qvalueCutoff= 0.05,readable= TRUE) 
head(ego)

#####################plot ####################################
###############################################################

pdf(file = "KEGG_PATH.pdf",width = 6, height = 5)
p<-ggplot(yy[1:4],aes(reorder(x=Description,-log10(qvalue)),y=-log10(qvalue)))+ coord_flip()+geom_bar(fill="grey91", stat="identity",width=.9)+theme_minimal()+ylab(expression(paste('-log10',italic(q),'value')))
p <- p + xlab("Kegg Pathway")+theme(axis.text = element_text(size =13),axis.title=element_text(size=13),axis.ticks = element_line(size = 2), axis.ticks.length = unit(.2, "cm")) 
#p<-p+theme(axis.ticks = element_line(size = 2),axis.line = element_line(size = 1, colour = "black", linetype=2))
p<-p+theme(axis.line = element_line(size = 1, colour = "black"))
p<-p+theme(axis.title.x = element_text(face="bold",size=13))+ theme(aspect.ratio = .4)
#p+geom_text(aes(label = Description), colour="blue", nudge_x = 0.3, hjust = 1) 
p
dev.off()
####################################################################################
####################################################################################
names(gl)<-make.names(NK_vs_Th$symbol, unique = T)
gl<-gl[order(gl, decreasing = T)]

#Import the reactome gene sets using msigdbr() with arguments category="C2" and subcategory="CP:REACTOME", or import from the provided file c2.cp.reactome.v2023.1.Hs.symbols.gmt using read.gmt().
reactome <-clusterProfiler::read.gmt("/home/ngiger/Documents/SIB_courses/EnrichmentAnalysis/data/c2.cp.reactome.v2023.2.Hs.symbols.gmt")


#Use the GSEA() function providing a sorted named vector of t-statistics and the reactome gene sets, and use argument minGSSize=30.
h_NK_vs_Th<-GSEA(gl, TERM2GENE = reactome,
                 minGSSize = 30,
                 pvalueCutoff = 0.05, #filter for only significant gene sets
                 seed=T)
#Count the number of significant adjusted p-values.
dim(h_NK_vs_Th)[1]
#Use barplot() and gseaplot() for the visualization of the results.

## barplot()
h_NK_vs_Th.sorted <- h_NK_vs_Th@result[order(h_NK_vs_Th@result$NES, decreasing = F),]
h_NK_vs_Th.sorted$color<-ifelse(h_NK_vs_Th.sorted$NES<0, "red", "darkblue")
par(mar=c(5,20,3,3))
n = 62
barplot(h_NK_vs_Th.sorted$NES,
        horiz = T, names=h_NK_vs_Th.sorted$Description,
        las=2, xlab="NES",
        cex.names = 0.5,
        col=h_NK_vs_Th.sorted$color) 
abline(v=0)
## gseaplot()
gseaplot(h_NK_vs_Th, geneSetID = h_NK_vs_Th.sorted$Description[1])

