library(Seurat)
library(tximport)
library(biomaRt)
files <- file.path("combinedSRR_output/alevin/quants_mat.gz")
file.exists(files)
txi <- tximport(files, type="alevin")
#replace Ensembl IDs with gene symbols (does not work due to some symbols being blank or "")
genes<-txi[["counts"]]@Dimnames[[1]]
genes <- sub("[.][0-9]*","",genes)
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
gene_IDs <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values = genes, mart= mart)
genedf<-as.data.frame(genes,col.names = "gene_id")
left_join(genedf, gene_IDs, by = c("genes"="ensembl_gene_id"))
genedf<-genedf %>% distinct(genes, .keep_all = TRUE)
txi[["counts"]]@Dimnames[[1]]<-genedf$hgnc_symbol
#run Seurat
pbmc <- CreateSeuratObject(counts = txi$counts , min.cells = 3, min.features = 200, project = "10X_PBMC")
#plot violin plot - no outstanding points
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
#plot features (molecules detected) vs Count (genes detected)
FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
pbmc <- subset(pbmc)
#subset based on genes per cell, and cells per gene
pbmc <- CreateSeuratObject(counts = txi$counts , min.cells = 3, min.features = 200)
#subset based on genes compared to molecules per cell
pbmc <- subset(pbmc, subset = nFeature_RNA<nCount_RNA | nFeature_RNA<2000)
#normalize data, scaling by 10000 first
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
#find top 2000 features based on different expression
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 500)
#plot variance
plot1 <- VariableFeaturePlot(pbmc)
#scale data
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
#perform and plot PCA
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
DimPlot(pbmc, reduction = "pca")
pbmc <- FindNeighbors(pbmc, dims = 1:5)
pbmc <- FindClusters(pbmc, resolution = 0.5)
#trying the same with umap
pbmc <- RunUMAP(pbmc, dims = 1:5)
#make pca results for pie chart 
pbmc@reductions[["pca"]]
#save object to transfer to analyst
saveRDS(pbmc, file = "pbmc_output.rds")

