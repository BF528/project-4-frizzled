library(dplyr)
library(Seurat)
library(patchwork)

# read in seurat file from programmer
panc <- readRDS("/projectnb2/bf528/users/frizzled/project_4/Programmer/panc.rds")

# rename cluster with cell type (analysis to do this is lower in script)
cluster_cell_types <- c("Delta_1","Alpha_1","Gamma","Unknown_A","Beta_1","Alpha_2","Unknown_B","Stellate","Beta_2","Acinar","Delta_2","Ductal","Vascular","Macrophage")
names(cluster_cell_types) <- levels(panc)
panc <- RenameIdents(panc, cluster_cell_types)

# identify marker genes for each cluster
# report only the positive ones = enriched marker genes
#marker_genes <- FindAllMarkers(panc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) - threshold gave too many genes
#marker_genes_test <- FindAllMarkers(panc, only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.5) - threshold gave too few genes

marker_genes <- FindAllMarkers(panc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5) # threshold gives just the right amount of genes
# min.pct argument requires a feature to be detected at a minimum percentage in either of the two groups of cells,

# significant marker genes
signficant_marker_genes <- marker_genes[marker_genes$p_val_adj<0.05,]

# write to csv
write.csv(signficant_marker_genes, "sig_marker_genes_cell_types.csv")

# top 2 markers for each cluster
top_2 <- signficant_marker_genes %>% group_by(cluster) %>% top_n(2, avg_log2FC)

# Assigning cell type identity to clusters

# label clusters as a cell type based on marker genes
violin <- VlnPlot(panc, features = c("GCG","INS","SST","PPY","GHRL","KRT19",
                                        "CPA1","PDGFRB","VWF","PECAM1","CD34",
                                        "CD163","CD68","IgG","CD3","CD8","TPSAB1",
                                        "KIT","CPA3"), pt.size=0)
# the following requested variables were not found: IgG, CD3, CD8

# violin without those genes
violin2 <- VlnPlot(panc, features = c("GCG","INS","SST","PPY","GHRL","KRT19",
                                        "CPA1","PDGFRB","VWF","PECAM1","CD34",
                                        "CD163","CD68","TPSAB1",
                                        "KIT","CPA3"), pt.size=0)
# feature plot
feature_plot <- FeaturePlot(panc, features = c("GCG","INS","SST","PPY","GHRL","KRT19",
                               "CPA1","PDGFRB","VWF","PECAM1","CD34",
                               "CD163","CD68","TPSAB1",
                               "KIT","CPA3"))

# Process to label clusters with cell types:
# Use supplemental table in paper and search each gene associated to each cell type in significant marker genes table
# write down clusters associated to that gene
# if a gene can not be assigned to a cluster because it is not found significant, then we assume it was filtered out by thresholds
# if there are multiple genes assigned significantly to clusters, then we compare the different sig expressed genes and perform a literature search to see if they resonate with one cell type more
# if there are no gene assigned to a cluster, than we do literature of other significantly expressed genes in that cluster 

## Alpha **
# 1, 5
alpha.markers <- FindMarkers(panc, ident.1 = 1, ident.2 = 5, min.pct = 0.5, min.diff.pct = 0.5, logfc.threshold = 0.5)
head(alpha.markers, n = 5)

# GCG (Alpha) was the most sig. expressed gene for cluster 1
# RPS6KA5
# CACNA1A - commonly found in pancreas, https://www.genecards.org/cgi-bin/carddisp.pl?gene=CACNA1A&keywords=CACNA1A
##https://discovery.lifemapsc.com/in-vivo-development/pancreas/dorsal-pancreatic-bud/endocrine-progenitor-cells
## "The first wave of endocrine progenitor cells (primary transition) mainly give rise to alpha cells and endocrine cells expressing multiple hormones"
### Final: GCG is alpha for cluster 1

# cluster 5
# CRYBA2 - expressed in pancreas - alpha, beta, and gamma - cluster 1 and 5 heavy
# AC092155.1 - 


## Beta **
# 4, 8
beta.markers <- FindMarkers(panc, ident.1 = 4, ident.2 = 8, min.pct = 0.5, min.diff.pct = 0.5, logfc.threshold = 0.5)
head(beta.markers, n = 5)

# INS (Beta) was the most sig. expressed gene for cluster 4
# different cells did not reveal anything
### Final: INS is beta for cluster 4

# Cluster 8
# MAFA - intestinal endocrine cells
## https://www.proteinatlas.org/ENSG00000182759-MAFA
## pancreatic beta-cell specific


## Delta **
# 0, 10
delta.markers <- FindMarkers(panc, ident.1 = 0, ident.2 = 10, min.pct = 0.5, min.diff.pct = 0.3, logfc.threshold = 0.5)
head(delta.markers, n = 5)
# SST (Delta) was the most sig. expressed gene for cluster 0
### Final: SST is Delta for cluster 0


## Gamma
# 10
# based on the log2fchange threshold, this is in the top 2 for cluster 10
### Final: PPY is Gamma for cluster 10


## Epsilon
# GHRL not found

## Ductal
# 9, 11 - just cluster 11
cluster9.markers <- FindMarkers(panc, ident.1 = 9, ident.2 = 11, min.pct = 0.25)
head(cluster9.markers, n = 5)
## KRT19 is sig in both cluster 11 and 9, but cpa1 (acinar) is only sig in 9
## Final: KRT19 is Ductal for cluster 11

# cluster 2
## SPTBN2 - found highly expressed in the pancreas and ranges in cell types, but Ductal comes up a lot
# https://www.proteinatlas.org/ENSG00000173898-SPTBN2/celltype
# also sig in cluster 5, and 10
## XACT - nothing found
#also sig in cluster 5, and 10
## CLUSTERS 2, 5, AND 20 ARE CLOSE

# cluster 3
## REG3A - pancreatic endocrine cells
# also sig in cluster 9
## REG1B - endocrine glandular cells
# also sig in cluster 9
## CLUSTER 3 AND 9 ARE CLOSE

# cluster 5
## AC092155.1
## CRYBA2
# also sig in cluster 1
## CLUSTER 1 AND 5 ARE CLOSE
cluster_5 <- FindMarkers(panc, ident.1 = 5, ident.2 = 1, min.pct = 0.5, min.diff.pct = 0.3, logfc.threshold = 0.5)
head(cluster_1, n = 5)

# cluster 6 --- Ductal
## TACSTD2 - found in various pancreatic cell types, but ductal was common
# https://www.proteinatlas.org/ENSG00000184292-TACSTD2/celltype/pancreas
# also sig in cluster 9 and 11 (ductal)
## KRT8 - found in various pancreatic cell types, but ductal was common
# https://www.proteinatlas.org/ENSG00000170421-KRT8/celltype/pancreas
# also sig in cluster 9 and 11 (ductal)
## CLUSTER 6, 9, AND 11 ARE CLOSE

cluster_6 <- FindMarkers(panc, ident.1 = 6, ident.2 = c(9,11), min.pct = 0.5, min.diff.pct = 0.3, logfc.threshold = 0.5)
head(cluster_6, n = 5)
# CFB - ductal cells - https://www.proteinatlas.org/ENSG00000243649-CFB/celltype/pancreas
# C3 - pancreatic endocrine cells and ductal cells - https://www.proteinatlas.org/ENSG00000125730-C3/celltype/pancreas

# cluster 8
## MAFA - intestinal endocrine cells
# https://www.proteinatlas.org/ENSG00000182759-MAFA
# pancreatic beta-cell specific
# also sig in cluster 4
## C1QL1 - pancreatic endocrine cells
# also sig in cluster 4
## CLUSTER 4 ANDS 8 ARE CLOSE

cluster_8 <- FindMarkers(panc, ident.1 = 8, ident.2 = 4, min.pct = 0.5, min.diff.pct = 0.3, logfc.threshold = 0.5)
head(cluster_8, n = 5)
# SORL1 - macrophage, ductal cells, and monocytes - https://www.proteinatlas.org/ENSG00000137642-SORL1/celltype/pancreas
#  exocrine cell types: acinar and duct cells.
# SMPD1 - same thing - https://www.proteinatlas.org/ENSG00000166311-SMPD1/celltype/pancreas
# pancreatic exocrine cells and monocyts - https://www.proteinatlas.org/ENSG00000104112-SCG3/celltype
# MAY JUST CALL IT OANCREATIC EXOCRINE CELLS


## Acinar
# 9
## Final: CPA1 is Acinar for cluster 9

## Stellate
# 7

## Vascular
# VWF - 12
# PECAM - 12
# cd34 not found 

## Macrophage
# cd163 not found
# cd 68 - cluster13
# igg not found

## Cytotoxic T
# cd3 not found
# cd8 not found

## Mast
# TPSAB1 - not found and other names for it not found
# KIT - not found and other names for it not found
# CPA3 - not found either

# In the paper "four immune cell types (tissue-resident macrophages, mast cells, B cells, and cytotoxic T cells); and Schwann cells"
# so macrophage, Cytotoxic, Mast, and B cells

# cluster 0 - delta_1
# cluster 1 - alpha_1
# cluster 2 - *gamma
# cluster 3 - Unknown_A
# cluster 4 - beta
# cluster 5 - *alpha_2
# cluster 6 - Unknown_B
# cluster 7 - *stellate
# cluster 8 - *beta, 
# cluster 9 - *acinar
# cluster 10 - *delta_2
# cluster 11 - *ductal
# cluster 12 - *vascular
# cluster 13 - *macrophage

cluster3_violin <- VlnPlot(panc, features = c("REG3A","REG1B"), pt.size=0)
cluster6_violin <- VlnPlot(panc, features = c("TACSTD2","KRT8"), pt.size=0)

# visualize the top marker genes per cluster
top_10 <- marker_genes %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top_5 <- marker_genes %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)

# heat map
heatmap_10 <- DoHeatmap(panc, features = top_10$gene, group.bar = TRUE, angle=70, size = 2)
heatmap_5 <- DoHeatmap(panc, features = top_5$gene, group.bar = TRUE, angle=70, size = 4)
png(filename = "heatmap.png", width = 1200, height = 800)
heatmap_5
dev.off()

# visualize the clustered cells using a projection method - umap
umap <- RunUMAP(panc, dims = 1:10)
umapp <- DimPlot(umap, reduction = "umap", label=TRUE)
png(filename = "umap.png", width = 1200, height = 800)
umapp
dev.off()

# find novel marker genes

# find markers for each cell type
#extract by clusters
cluster0.markers <- FindMarkers(panc, ident.1 = 'Delta_1', min.pct = 0.25, only.pos = TRUE, logfc.threshold = 0.5)
cluster1.markers <- FindMarkers(panc, ident.1 = 'Alpha_1', min.pct = 0.25, only.pos = TRUE, logfc.threshold = 0.5)
cluster2.markers <- FindMarkers(panc, ident.1 = 'Gamma', min.pct = 0.25, only.pos = TRUE, logfc.threshold = 0.5)
cluster3.markers <- FindMarkers(panc, ident.1 = 'Unknown_A', min.pct = 0.25, only.pos = TRUE, logfc.threshold = 0.5) 
cluster4.markers <- FindMarkers(panc, ident.1 = 'Beta_1', min.pct = 0.25, only.pos = TRUE, logfc.threshold = 0.5) 
cluster5.markers <- FindMarkers(panc, ident.1 = 'Alpha_2', min.pct = 0.25, only.pos = TRUE, logfc.threshold = 0.5) 
cluster6.markers <- FindMarkers(panc, ident.1 = 'Unknown_B', min.pct = 0.25, only.pos = TRUE, logfc.threshold = 0.5)
cluster7.markers <- FindMarkers(panc, ident.1 = 'Beta_2', min.pct = 0.25, only.pos = TRUE, logfc.threshold = 0.5)
cluster8.markers <- FindMarkers(panc, ident.1 = 'Acinar', min.pct = 0.25, only.pos = TRUE, logfc.threshold = 0.5)
cluster9.markers <- FindMarkers(panc, ident.1 = 'Delta_2', min.pct = 0.25, only.pos = TRUE, logfc.threshold = 0.5)
cluster10.markers <- FindMarkers(panc, ident.1 = 'Ductal', min.pct = 0.25, only.pos = TRUE, logfc.threshold = 0.5)
cluster11.markers <- FindMarkers(panc, ident.1 = 'Vascular', min.pct = 0.25, only.pos = TRUE, logfc.threshold = 0.5)
cluster12.markers <- FindMarkers(panc, ident.1 = 'Macrophage', min.pct = 0.25, only.pos = TRUE, logfc.threshold = 0.5)

# top 6 differentially expressed genes for each cell type
top_6_marker <- signficant_marker_genes %>% group_by(cluster) %>% top_n(6, avg_log2FC)

# stricter min percent
novel_genes <- FindAllMarkers(panc, only.pos = TRUE, min.pct = 0.8, logfc.threshold = 1)
# stricter adjusted p value
signficant_novel_genes <- novel_genes[novel_genes$p_val_adj<0.005,]
write.csv(signficant_novel_genes, "novel_genes.csv")
