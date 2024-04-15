setwd("~/temporary/")
list.files()


library(dplyr)
library(Seurat)
library(patchwork)
library(harmony)


adata <- Read10X(data.dir = "count_matrix/")
adata <- CreateSeuratObject(counts = adata, 
                            project = "endometrium_epithelial", min.cells = 3, min.features = 200)
adata

adata[["percent.mt"]] <- PercentageFeatureSet(adata, pattern = "^MT-")
p1 <- VlnPlot(adata, features =c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


adata <- subset(adata, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 15) 
adata <- NormalizeData(adata, normalization.method = "LogNormalize", scale.factor = 10000)
adata <- FindVariableFeatures(adata, selection.method = "vst", nfeatures = 3000)
top10 <- head(VariableFeatures(adata), 10)

p2 <- VariableFeaturePlot(adata)
p3 <- LabelPoints(plot = p2, points = top10, repel = TRUE)

adata <- ScaleData(adata, vars.to.regress = c("nCount_RNA",
                                          "percent.mt"), do.par = TRUE, num.cores = 3)

regev_lab_cell_cycle_genes <- read.delim("~/scRNA/files/regev_lab_cell_cycle_genes.txt",header=F)
s.genes <- regev_lab_cell_cycle_genes[1:43,1]
g2m.genes <- regev_lab_cell_cycle_genes[44:97,1]
adata <- CellCycleScoring(object = adata,
                          s.features = s.genes,
                          g2m.features = g2m.genes,
                          set.ident = TRUE)

adata <- RunPCA(adata, features = c(s.genes, g2m.genes))

p4 <- VizDimLoadings(adata, dims = 1:2, reduction = "pca")  ## visualise features 
p5 <- DimPlot(adata) ## visualise cells 
p6 <- DimHeatmap(adata, dims = 1:9, cells = 500, balanced = TRUE) ## explore if feature is included in downstream analysis 
p7 <- ElbowPlot(adata) ## determine dimensionality of the dataset

adata <- ScaleData(adata,vars.to.regress = c("nCount_RNA", "percent.mt", "S.Score", "G2M.Score"),
                   do.par = TRUE, num.cores = 3)

p8 <- VizDimLoadings(adata, dims = 1:2, reduction = "pca")  ## visualise features 
p9 <- DimPlot(adata) ## visualise cells 
p10 <- DimHeatmap(adata, dims = 1:15, cells = 500, balanced = TRUE) ## explore if feature is included in downstream analysis 
p11 <- ElbowPlot(adata) ## determine dimensionality of the dataset

adata <- RunPCA(adata, features = VariableFeatures(object = adata))

p12 <- VizDimLoadings(adata, dims = 1:2, reduction = "pca")  ## visualise features 
p13 <- DimPlot(adata) ## visualise cells 
p14 <- DimHeatmap(adata, dims = 1:15, cells = 500, balanced = TRUE) ## explore if feature is included in downstream analysis 
p15 <- ElbowPlot(adata) ## determine dimensionality of the dataset

adata <- FindNeighbors(adata, dims = 1:15)
adata <- FindClusters(adata, resolution = 0.5)

adata <- RunUMAP(adata, dims = 1:15)
p16 <- DimPlot(adata, reduction = "umap", group.by = "old.ident")

metadata <- read.csv("~/temporary/metadata.csv", row.names = 1)
metadata_old <- adata@meta.data
metadata <- merge(metadata, metadata_old, by="row.names")
rownames(metadata) <- metadata$Row.names
metadata$Row.names <- NULL
metadata <- metadata[rownames(adata@meta.data), ]
adata@meta.data <- metadata

var_to_regress <- c("old.ident", "BiopsyType", "DonorID", "Location")
adata <- RunHarmony(adata, 
                  group.by.vars = var_to_regress, 
                  reduction = "pca", 
                  assay.use = "RNA", 
                  reduction.save = "harmony", 
                  lambda = length(var_to_regress))
harmony_stdev <- adata@reductions$harmony@stdev
harmony_above_2 <- sum(harmony_stdev > 2)

adata <- FindNeighbors(adata, dims = 1:harmony_above_2, reduction="harmony")
adata <- FindClusters(adata, resolution = 0.5, reduction="harmony")
adata <- RunUMAP(adata, dims = 1:15, reduction="harmony")
table(adata$orig.ident, adata@meta.data$seurat_clusters)
Idents(object = adata) <- adata@meta.data$seurat_clusters

p17 <- DimPlot(adata, reduction = "umap", group.by = "old.ident", label = T, pt.size = 0.8) + NoLegend()
p18 <- DimPlot(adata, reduction = "umap", label = T, pt.size = 0.8) + NoLegend()
p19 <- DimPlot(adata, reduction = "umap", label = T, pt.size = 0.5, group.by = "Epithelial.celltype") + NoLegend()


saveRDS(adata, file = "endometrium_epithelial.rds")








