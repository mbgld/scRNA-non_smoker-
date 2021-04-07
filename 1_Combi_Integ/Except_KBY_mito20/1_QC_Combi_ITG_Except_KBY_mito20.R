# Get previous dataset
#load("/media/yschang/W/Merged/1_Combi_Integ/Except_KBY_mito20/Results/1_QC_Combi_ITG_Except_KBY_mito20.RData")
load("F:/Merged/1_Combi_Integ/Except_KBY_mito20/Results/1_QC_Combi_ITG_Except_KBY_mito20.RData")

# set up
library(Seurat); library(dplyr)
# YS_ITG.outdir = "/media/yschang/W/Merged/1_Combi_Integ/Except_KBY_mito20/Results/"
YS_ITG.outdir = "F:/Merged/1_Combi_Integ/Except_KBY_mito20/Results/"
### Part A. Merging YS cohorts

## 1. Call each dataset
# KJA dataset
# NL
KJA_NL_raw.data <- Read10X(data.dir = "/media/yschang/T/KJA/A_data_NL/filtered_feature_bc_matrix/"); KJA_NL_original_seurat.obj <- CreateSeuratObject(counts = KJA_NL_raw.data, project = "KJA_NL", min.cells = 3, min.features = 200); rm(KJA_NL_raw.data)
# Tu
KJA_Tu_raw.data <- Read10X(data.dir = "/media/yschang/T/KJA/A_data_Tu/filtered_feature_bc_matrix/"); KJA_Tu_original_seurat.obj <- CreateSeuratObject(counts = KJA_Tu_raw.data, project = "KJA_Tu", min.cells = 3, min.features = 200); rm(KJA_Tu_raw.data)

# LYJ dataset
# NL
LYJ_NL_raw.data <- Read10X(data.dir = "/media/yschang/T/LYJ/A_data_NL/filtered_feature_bc_matrix/"); LYJ_NL_original_seurat.obj <- CreateSeuratObject(counts = LYJ_NL_raw.data, project = "LYJ_NL", min.cells = 3, min.features = 200); rm(LYJ_NL_raw.data)
# Tu
LYJ_Tu_raw.data <- Read10X(data.dir = "/media/yschang/T/LYJ/A_data_Tu/filtered_feature_bc_matrix/"); LYJ_Tu_original_seurat.obj <- CreateSeuratObject(counts = LYJ_Tu_raw.data, project = "LYJ_Tu", min.cells = 3, min.features = 200); rm(LYJ_Tu_raw.data)

# KBY dataset
# NL
# KBY_NL_raw.data <- Read10X(data.dir = "/media/yschang/T/KBY/A_data_NL/filtered_feature_bc_matrix/"); # KBY_NL_original_seurat.obj <- CreateSeuratObject(counts = KBY_NL_raw.data, project = "KBY_NL", min.cells = 3, min.features = 200); rm(KBY_NL_raw.data)
# Tu
# KBY_Tu_raw.data <- Read10X(data.dir = "/media/yschang/T/KBY/A_data_Tu/filtered_feature_bc_matrix/"); KBY_Tu_original_seurat.obj <- CreateSeuratObject(counts = KBY_Tu_raw.data, project = "KBY_Tu", min.cells = 3, min.features = 200); rm(KBY_Tu_raw.data)

# SSS dataset
# Tu
SSS_Tu_raw.data <- Read10X(data.dir = "/media/yschang/S/SSS/A_data_Tu/filtered_feature_bc_matrix/"); SSS_Tu_original_seurat.obj <- CreateSeuratObject(counts = SSS_Tu_raw.data, project = "SSS_Tu", min.cells = 3, min.features = 200); rm(SSS_Tu_raw.data)

# YJO dataset
# NL
YJO_NL_raw.data <- Read10X(data.dir = "/media/yschang/S/YJO/A_data_NL/filtered_feature_bc_matrix/"); YJO_NL_original_seurat.obj <- CreateSeuratObject(counts = YJO_NL_raw.data, project = "YJO_NL", min.cells = 3, min.features = 200); rm(YJO_NL_raw.data)
# Tu
YJO_Tu_raw.data <- Read10X(data.dir = "/media/yschang/S/YJO/A_data_Tu/filtered_feature_bc_matrix/"); YJO_Tu_original_seurat.obj <- CreateSeuratObject(counts = YJO_Tu_raw.data, project = "YJO_Tu", min.cells = 3, min.features = 200); rm(YJO_Tu_raw.data)

# NHK dataset
#	NL
NKH_NL_raw.data <- Read10X(data.dir = "/media/yschang/S/NKH/A_data_NL/filtered_feature_bc_matrix/"); NKH_NL_original_seurat.obj <- CreateSeuratObject(counts = NKH_NL_raw.data, project = "NKH_NL", min.cells = 3, min.features = 200); rm(NKH_NL_raw.data)
# Tu
NKH_Tu_raw.data <- Read10X(data.dir = "/media/yschang/S/NKH/A_data_Tu/filtered_feature_bc_matrix/"); NKH_Tu_original_seurat.obj <- CreateSeuratObject(counts = NKH_Tu_raw.data, project = "NKH_Tu", min.cells = 3, min.features = 200); rm(NKH_Tu_raw.data)

# LDS dataset
# NL
LDS_NL_raw.data <- Read10X(data.dir = "/media/yschang/S/LDS/A_data_NL/filtered_feature_bc_matrix/"); LDS_NL_original_seurat.obj <- CreateSeuratObject(counts = LDS_NL_raw.data, project = "LDS_NL", min.cells = 3, min.features = 200); rm(LDS_NL_raw.data)
# TU
LDS_Tu_raw.data <- Read10X(data.dir = "/media/yschang/S/LDS/A_data_Tu/filtered_feature_bc_matrix/"); LDS_Tu_original_seurat.obj <- CreateSeuratObject(counts = LDS_Tu_raw.data, project = "LDS_Tu", min.cells = 3, min.features = 200); rm(LDS_Tu_raw.data)


## 2. Merge
YS_Combi_seurat.obj <- merge(KJA_NL_original_seurat.obj, y = c( KJA_Tu_original_seurat.obj, LYJ_NL_original_seurat.obj, LYJ_Tu_original_seurat.obj, SSS_Tu_original_seurat.obj, YJO_NL_original_seurat.obj, YJO_Tu_original_seurat.obj, NKH_NL_original_seurat.obj, NKH_Tu_original_seurat.obj, LDS_NL_original_seurat.obj, LDS_Tu_original_seurat.obj), add.cell.ids = c('KJA_NL', 'KJA_Tu', 'LYJ_NL', 'LYJ_Tu', 'SSS_Tu', 'YJO_NL', 'YJO_Tu', 'NKH_NL', 'NKH_Tu', 'LDS_NL', 'LDS_Tu')); rm(KJA_NL_original_seurat.obj, KJA_Tu_original_seurat.obj, LYJ_NL_original_seurat.obj, LYJ_Tu_original_seurat.obj, SSS_Tu_original_seurat.obj, YJO_NL_original_seurat.obj, YJO_Tu_original_seurat.obj, NKH_NL_original_seurat.obj, NKH_Tu_original_seurat.obj, LDS_NL_original_seurat.obj, LDS_Tu_original_seurat.obj)


## 3. QC (1. percent.mito:  20 <; 2. UMIs count: 100 ~ 150,000; 3. gene count: 200 ~ 10,000)
# Get mito percentage and QC
YS_Combi_seurat.obj[["percent.mt"]] <- PercentageFeatureSet(YS_Combi_seurat.obj, pattern = "^MT-")
YS_Combi.sqj <- subset(x = YS_Combi_seurat.obj, subset = percent.mt < 20 & nCount_RNA >100 & nCount_RNA < 150000 & nFeature_RNA > 200 & nFeature_RNA < 10000)
# Visualize
VlnPlot(object = YS_Combi.sqj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) 
plot1 <- FeatureScatter(YS_Combi.sqj, feature1 = "nCount_RNA", feature2 = "percent.mt"); plot2 <- FeatureScatter(YS_Combi.sqj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA"); CombinePlots(plots=list(plot1, plot2))


## 4. LogNormalization
YS_Combi.snj = NormalizeData(object = YS_Combi.sqj, normalization.method = "LogNormalize", scale.factor = 10000); rm(YS_Combi.sqj) 


## 5. Scaling
# Variable feature selection; default nfeatures=2000. 
YS_Combi.snj <- FindVariableFeatures(YS_Combi.snj, selection.method = "vst", nfeatures = 2000); v.genes <- VariableFeatures(YS_Combi.snj) # Visualize
LabelPoints(plot = VariableFeaturePlot(YS_Combi.snj), points = v.genes[1:15],repel = T)
#	Scaling
YS_Combi.ssj <- ScaleData(object = YS_Combi.snj, features = rownames(YS_Combi.snj)); rm(YS_Combi.snj)


## 6. Dimension Reduction
# Linear Dimensional Reduction by PCA. 
YS_Combi.spj <- RunPCA(object = YS_Combi.ssj, features = v.genes); rm(YS_Combi.ssj)
# Print PCA components
print(YS_Combi.spj[["pca"]], dims = 1:5, nfeatures=5) 
# Visualize data
VizDimLoadings(YS_Combi.spj, dims = 1:2, reduction ="pca")
DimPlot(YS_Combi.spj, reduction = 'pca')
DimHeatmap(YS_Combi.spj, dims = 1, cells = 500, balanced =TRUE)
DimHeatmap(YS_Combi.spj, dims = 1:15, cells = 500, balanced =TRUE)

## NON-LINEAR DIMENSION REDUCTION 
ElbowPlot(object = YS_Combi.spj)
YS_Combi.spj <- JackStraw(object = YS_Combi.spj, num.replicate = 100)
YS_Combi.spj <- ScoreJackStraw(object = YS_Combi.spj, dims = 1:20)
JackStrawPlot(object = YS_Combi.spj, dims = 1:20)


## 7. Clustering analysis
### 7.1 Determine the number of cluster
YS_Combi.spj <- FindNeighbors(object = YS_Combi.spj, dims = 1:10) # GRAPH based clustering. Defaults for dims is 10. 
YS_Combi.spj <- FindClusters(object = YS_Combi.spj, resolution = 0.2) # resolution is critical for cluster number: use 0.2 ~ 1.2. small number increase cluster number whereas large number decreases. 
### 7.2 Determine distance between clusters using UMAP #################
YS_Combi.suj <- RunUMAP(object = YS_Combi.spj, dims = 1:10)  # Determine the distance between cluster using dims. Default value is 10.

## 8. Visualize
DimPlot(object = YS_Combi.suj, reduction = "umap", label = T)
DimPlot(object = YS_Combi.suj, group.by= 'orig.ident', label = T)

##### END of merging Every cases ###############################################
################################################################################

##### Part B. Integration of Dataset
### 1. Dataset Preprocessing 
YS_dataset.list <- SplitObject(YS_Combi.suj, split.by = 'orig.ident')
YS_dataset.list <- lapply(YS_dataset.list, FUN = function(x){
  x<-NormalizeData(x)
  x<-FindVariableFeatures(x, selection.method = 'vst', nfeatures=2000)
})
YS.anchors <- FindIntegrationAnchors(object.list = YS_dataset.list, dims=1:30)

### 2. Integration
YS_Itg.sbj <- IntegrateData(anchorset=YS.anchors, dims=1:30)
DefaultAssay(YS_Itg.sbj) <- 'integrated'; rm(YS_dataset.list, YS.anchors)

## 3. Scaling: .sbj starting point
YS_Itg.ssj <- ScaleData(object = YS_Itg.sbj, features = rownames(YS_Itg.sbj))

## 4. Dimensional Reduction
# 4.1 PCA. 
YS_Itg.spj <- RunPCA(object = YS_Itg.ssj, npcs = 30, verbose = FALSE); rm(object = YS_Itg.ssj); print(YS_Itg.spj[["pca"]], dims = 1:5, nfeatures=5)
# Visualize data
VizDimLoadings(YS_Itg.spj, dims = 1:2, reduction ="pca"); DimPlot(YS_Itg.spj, reduction = 'pca'); DimHeatmap(YS_Itg.spj, dims = 1, cells = 500, balanced =TRUE); DimHeatmap(YS_Itg.spj, dims = 1:15, cells = 500, balanced =TRUE)
# 4.2  NON-LINEAR DIMENSION REDUCTION
ElbowPlot(object = YS_Itg.spj)
YS_Itg.spj <- JackStraw(object = YS_Itg.spj, num.replicate = 100)
YS_Itg.spj <- ScoreJackStraw(object = YS_Itg.spj, dims = 1:20)
JackStrawPlot(object = YS_Itg.spj, dims = 1:20)

### 5. Determine the number of cluster
YS_Itg.spj <- FindNeighbors(object = YS_Itg.spj, dims = 1:10) 
YS_Itg.spj <- FindClusters(object = YS_Itg.spj, resolution = 0.2)

### 6. Determine distance between clusters using UMAP
YS_Itg.suj <- RunUMAP(object = YS_Itg.spj, reduction = 'pca',  dims = 1:10)

### 7. Visualize
DimPlot(object = YS_Itg.suj, reduction = "umap", label = T)
DimPlot(object = YS_Itg.suj, group.by= 'orig.ident', label = T)


##################################################
save.image(paste0(YS_ITG.outdir, "1_QC_Combi_ITG_Except_KBY_mito20.RData"))

rm(list=ls())
gc()
q()
