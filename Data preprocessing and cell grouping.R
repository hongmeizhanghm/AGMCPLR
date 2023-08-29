library(Seurat)

Seurat_tmp <- CreateSeuratObject(counts = sc_dataset,project="sc",min.cells=0.1*ncol(sc_dataset), min.features = 0)
Seurat_tmp <- NormalizeData(object = Seurat_tmp, normalization.method = "LogNormalize", scale.factor = 10000, verbose = F)
Seurat_tmp <- FindVariableFeatures(Seurat_tmp, selection.method = "vst", verbose = F)
Seurat_tmp <- ScaleData(Seurat_tmp, verbose = F)
Seurat_tmp <- RunPCA(Seurat_tmp, features = VariableFeatures(Seurat_tmp), verbose = F)
Seurat_tmp <- FindNeighbors(Seurat_tmp, dims = 1:10, verbose = F)
Seurat_tmp <- FindClusters(object = Seurat_tmp, resolution = 0.6, verbose = F)
cell_group <- Idents(Seurat_tmp)

