# ti-scMR

Trajectory-inference-based dynamic single-cell Mendelian randomization (ti-scMR)

## Tutorial
Here we used a simulated toy dataset to illustrate the workflow of ti-scMR. We first process the sc-RNA count matrix using Seurat v5.
```R
library(Seurat)
seurat_obj <- readRDS("example_sc.rds")

# Quality Control
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-",assay = 'RNA')
seurat_obj[["percent.rb"]] <- PercentageFeatureSet(seurat_obj, pattern = "^RP[SL]",assay = 'RNA')
seurat_obj$log10GenesPerUMI <- log10(seurat_obj$nFeature_RNA) / log10(seurat_obj$nCount_RNA)
seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & percent.mt < 5 & percent.rb < 5)  #nUMI = nCount_RNA,nGene = nFeature_RNA

# Log normalization
seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize")

# High variable features
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 3000)

# Batch correction
library(harmony)
seurat_obj <- RunHarmony(seurat_obj,"id")

# Dimension reduction
seurat_obj <- seurat_obj %>% RunUMAP(reduction="harmony",dims=1:20) %>%
  FindNeighbors(reduction="harmony",dims=1:20)%>%
  FindClusters(resolution = 0.5)

# Scaling
seurat_obj <- ScaleData(seurat_obj)

# Clustering

# Trajectory inference
library(SingleCellExperiment)
library(slingshot)
library(scater)
library(scran)
library(bluster)
sce <- as.SingleCellExperiment(seurat_A)
sce <- slingshot(sce,reducedDim="UMAP")
```
