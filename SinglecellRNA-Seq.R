 # Analysis of Single cell RNA-Seq data

# Load libraries

library(Seurat)
library(tidyverse)

# Load the NSCLC dataset.

nslc.sparse.m <- Read10X_h5(filename =
                            "20K_NSLC_DTC_3p_nextgem_Multiplex_count_raw_feature_bc_matrix.h5")
str(nslc.sparse.m)
cts <- nslc.sparse.m$'Gene Expression'


# Initialize the Seurat object with the raw (non-normalized data).

nslc.sparse.obj <- CreateSeuratObeject(counts = cts, project = 'NSCLC', min.cells = 3,
                                      min.features = 200)
str(nslc.sparse.obj)


# 1. QC --------------------------------
# % MT reads

nslc.sparse.obj[["persent.mt"]] <- PercentageFeatureSet(nslc.sparse.obj, pattern = "^MT-")
vieew(nslc.sparse.obj@meet.data)

VlnPlot(nslc.sparse.obj, feature1 = c("nFeature_RNA","nCount_RNA","percent.mt"), ncol = 3)
FeatureScatter(nslc.sparse.obj, feature1 ="nCount_RNA", feature2 = "nFeature_RNA") +
     geom_smooth(method = 'lm')
     
# 2. Filtering ------------------------ 

nslc.sparse.obj <- subset(nslc.sparse.obj, subset = nFeature_RNA >200 & nFeature_RNA < 2500 &
                           percent.mt < 5)
                           
# 3. Normalize data ---------------------------

nslc.sparse.obj <- NormalizedData(nslc.sparse.obj)

# 4. Identify highly variable features --------------

nslc.sparse.obj <- FindVariableFeature(nslc.sparse.obj, selection.method = 'vst', 
                                         nfeature = 2000) 

# identify the 10 most highly variable genes

top10 <-  head(VariableFeatures(nslc.sparse.obj))

# plot variable features with and without labels 

plot1 <- VariableFeaturePlot(nslc.sparse.obj)
LabelPoints(plot =plot1, points = top10, repel = TRUE)

# 5. Scaling -----------------------------------------

all.genes <- rownames(nslc.sparse.obj)
nslc.sparse.obj <- ScaleData(nslc.sparse.obj, features = all.genes)

str(nslc.sparse.obj)

# 6. Perform Linear dimensionality reduction ------------------------

nslc.sparse.obj <- RunPCA(nslc.sparse.obj, features = VariableFeatures(object = nslc.sparse.obj))

# visualize PCA results

print(nslc.sparse.obj[["pca"]], dims = 1.5, nfeatures = 5)
DimHeatmap(nslc.sparse.obj, dims = 1, cells = 500, balanced = TRUE)

# determine dimensionality of the data

ElbowPlot(nslc.sparse.obj)

# 7. Clusting -------------------------------------------

nslc.sparse.obj <- FindNeighbors(nslc.sparse.obj, dims=1:15)

# understanding resolution

nslc.sparse.obj <- FindClusters(nslc.sparse.obj, resolution = c(0.1,0.3,0.5,0.7, 1))
view(nslc.sparse.obj@meta.data)

DimPlot(nslc.sparse.obj, grouped.by = "RNA_snn_res.0.5", label = TRUE)

# setting identity of clusters

Idents(nslc.sparse.obj)
Idents(nslc.sparse.obj) <- "RNA_snn_res.0.1"
Idents(nslc.sparse.obj)

# non-linear dimensionally reduction-----------------------------
# If you haven't installed UMAP, you can do so via reticulate::py_install(packages = "umap-learn")

RNA_snn_res <- RunUMAP(RNA_snn_res, dims = 1:15)

# note that you can set 'label = TRUE' or use the labelClusters function to help label individaul clusters

DimPlot(RNA_snn_res, reduction = "umap")


