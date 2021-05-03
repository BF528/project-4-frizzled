#---------------------------------
#Title: BF528 Project 4 Programmer
#Author: Janvee Patel
#Date: 04/12/2021
#---------------------------------

#-------------------------------------------------------------------------------------------------------------------------
#Additional Information about the data in the Seurat Object taken from Satijalab
#Citation: Tools for single cell genomics. (n.d.). https://satijalab.org/seurat/

#nCount_RNA = number of molecules detected within a cell
#nFeature_RNA = number of genes detected in each cell
#low nFeature_RNA for a cell = dying cell/empty droplet
#high nCount_RNA and/or n_Feature_RNA = cell might be doublet/multiplet
#-------------------------------------------------------------------------------------------------------------------------

#load libraries
library(Seurat)
library(tximport)
library(fishpond)
library(plyr)
library(dplyr)
library(Matrix)
library(RColorBrewer)

#Part 1
#Load the UMI counts matrix - Salmon Quant Alevin file
files <- file.path("/projectnb/bf528/users/frizzled/project_4/Data_Curator/Alevin_Output/Whitelist_mean_counts/alevin/quants_mat.gz")
txi <- tximport(files, type="alevin")
txicts <- txi$counts

#convert the Ensembl ids to Gene Symbols
ens_ids <- rownames(txicts)

#check what txi$counts rownames look like
txicts[1:10, 1:2]

#Downloaded Ensembl Mart CSV file for Human Genome with Ensembl Ids, Ensembl Version Ids, Gene Names, and HGNC Symbols
#Remove the 3 duplicates of Ensembl IDs mapping to the same gene symbol
#in the CSV file: ENSG00000230417.12, ENSG00000254876.5, ENSG00000276085.1
#duplicated gene symbols were made unique with "."
id_map <- read.csv(file="mart_export.txt")
id_map <- id_map[, c(2,3)]
id_map <- id_map[!duplicated(id_map$Gene.stable.ID.version), ]
#id_map$Gene.name <- make.unique(id_map$Gene.name, sep=".")

#Map the Ensembl ids to the Gene Symbols using the Ensembl Mart CSV mapping file
#join() from Plyr does not reorder the rows
ens_ids <- as.data.frame(ens_ids)
colnames(ens_ids) <- c("Gene.stable.ID.version")
mapped_ids <- join(ens_ids, id_map, by="Gene.stable.ID.version")

#set the rownames as the gene symbols
txicts <- txicts[rownames(txicts) %in% mapped_ids$Gene.stable.ID.version, ]
rownames(txicts) <- mapped_ids$Gene.name

#double check what the rownames of txicts look like
txicts[1:10, 1:2]

#check initial values of the number of genes and cells
dim(txicts)

#---------------------------------------------------------------------------------------------------------
#consider minimum # of count of genes/cell (complexity of cells)
#consider minimum # of count of cells/gene (rarity of genes)
#Applied an approach that is used in an MIT tutorial
#Citation: Single-cell RNA-seq Demo (10x non-small cell lung cancer). (n.d.). http://barc.wi.mit.edu/education/hot_topics/scRNAseq_2020/SingleCell_Seurat_2020.html
genes_per_cell <- Matrix::colSums(txicts > 0)
cells_per_gene <- Matrix::rowSums(txicts > 0)

hist(log10(genes_per_cell+1),main='Genes per Cell')
hist(log10(cells_per_gene+1),main='Cells per Gene')
plot(sort(genes_per_cell), xlab='Cells', log='y', ylab="Genes per Cell (sorted)", main='Complexity of Cells (Genes per Cell Ordered)')

#create Seurat Object
#studying dataset of pancreatic cells
panc <- CreateSeuratObject(counts = txicts, project = "panc_sc", min.cells = 3, min.features = 50)
panc

#see mitochondrial QC metrics
panc[["percent.mt"]] <- PercentageFeatureSet(panc, pattern = "^MT-")
x <- panc[["percent.mt"]]

#see QC metrics for the first 10 cells
head(panc@meta.data, 10)

#visualize QC metrics as violin plots
VlnPlot(panc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#visualize feature-feature relationships
plot1 <- FeatureScatter(panc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(panc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- FeatureScatter(panc, feature1 = "nFeature_RNA", feature2 = "percent.mt")
p4 <- plot1 + plot2 + plot3

#filter out low-quality cells
panc <- subset(panc, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 48)

#determine number of cells after filtering out low-quality cells
ncol(x=panc)

#Check plots after filtering out low-quality cells
VlnPlot(panc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(panc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(panc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- FeatureScatter(panc, feature1 = "nFeature_RNA", feature2 = "percent.mt")
plot1 + plot2 + plot3


#----------------------------------------------------------------------------------------------------------------------------
#Part 2
#Normalize counts matrix - ensure expression values across cells are on comparable scale (Satijalab)
#divide counts for each gene by total counts in cell -> multiply this value for each gene by scale factor and log transform
panc <- NormalizeData(panc, normalization.method = "LogNormalize", scale.factor = 10000)

#Feature Selection
#Filter to only include highly variable features - filter out low variance genes
#features with high cell-to-cell variation (highly expressed in certain cells, lowly expressed in other cells) (Satijalab)
panc <- FindVariableFeatures(panc, selection.method = "vst", nfeatures = 2000)

#Top 10 most highly variable genes
top10 <- head(VariableFeatures(panc), 10)
top10

#Plot variable features
plot4 <- VariableFeaturePlot(panc) #without labels
plot5 <- LabelPoints(plot = plot4, points=top10, repel=TRUE, xnudge = 0, ynudge = 0) #with labels
plot4
plot5

#Scaling (Satijalab)
#linear transformation prior to dimensional reduction
#adjust expression of each gene -> mean expression across cells is 0
#scales expression of each gene -> variance across cells is 1
all_genes <- rownames(panc)
panc <- ScaleData(panc, features=all_genes)

#Linear Dimensional Reduction
#Perform PCA
panc <- RunPCA(panc, features=VariableFeatures(object=panc))

#Method 1 - Look at PCA results
print(panc[["pca"]], dims = 1:5, nfeatures = 5)
DimPlot(panc, reduction=pca)

#Method 2 - Look at PCA results
VizDimLoadings(panc, dims = 1:2, reduction = "pca")

#Method 3 - Look at PCA results
DimPlot(panc, reduction = "pca")

#Method 4 - DimHeatmap
DimHeatmap(panc, dims = 1, cells = 500, balanced = TRUE)

#Determine Dimensionality of Dataset
panc <- JackStraw(panc, num.replicate = 100)
panc <- ScoreJackStraw(panc, dims = 1:20)
JackStrawPlot(panc, dims=1:20)

#Plot elbow plot
ElbowPlot(panc)

#---------------------------------------------------------------------------------------------------------------------------
#Part 3
#Determine clusters
panc <- FindNeighbors(panc, dims=1:10)
panc <- FindClusters(panc, resolution = 0.5)

#Cluster IDs for first 5 cells
head(Idents(panc), 5)

#Non-Linear Dimensional Reduction
panc <- RunUMAP(panc, dims=1:10)
DimPlot(panc, reduction = "umap")

#Save to RDS
saveRDS(panc, "panc.rds")

#Create pie chart to display relative proportions of count in each cluster
#stores the cluster number and count of cell number
cluster_counts <- as.data.frame(table(Idents(panc)))
cluster_counts <- paste(cluster_counts$Var1, cluster_counts$Freq, sep=" - ")

#stores the cluster number and the relative proportions
cluster_proportions <- as.data.frame(prop.table(table(Idents(panc))))
cluster_proportions$Freq <- cluster_proportions$Freq %>% `*`(100) %>% round(3)
cluster_proportions$Percents <- cluster_proportions$Freq
cluster_proportions$Percents <- paste(cluster_proportions$Percents, "%", sep="")

#plot the piechart
cols = colorRampPalette(brewer.pal(8, "Dark2"))(14)
pie(cluster_proportions$Freq, labels=cluster_proportions$Percents, col=cols, main = "Relative Proportions of Cell Numbers\nFor the Identified Clusters")
legend(1.45, 1, legend=cluster_counts, cex = 0.8, fill = cols, bty="n", title="Clusters")
