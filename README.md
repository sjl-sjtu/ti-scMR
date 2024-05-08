# ti-scMR

Trajectory-inference-based dynamic single-cell Mendelian randomization (ti-scMR)

## Installation
The package can be installed by

## Tutorial
Here we used a simulated toy dataset to illustrate the workflow of ti-scMR. We first process the sc-RNA count matrix using `Seurat` v5.
```R
library(data.table)
library(tidyverse)
library(Seurat)

seurat_obj <- readRDS("example_sc.rds")

# Quality Control
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-",assay = 'RNA')
seurat_obj[["percent.rb"]] <- PercentageFeatureSet(seurat_obj, pattern = "^RP[SL]",assay = 'RNA')
seurat_obj$log10GenesPerUMI <- log10(seurat_obj$nFeature_RNA) / log10(seurat_obj$nCount_RNA)
seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & percent.mt < 5 & percent.rb < 5)

# Log normalization
seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)

# High variable features
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)

# Scaling
seurat_obj <- ScaleData(seurat_obj)

# PCA
seurat_obj <- RunPCA(seurat_obj)

# Batch correction
library(harmony)
seurat_obj <- RunHarmony(seurat_obj,"id")

# Dimension reduction
seurat_obj <-  RunUMAP(seurat_obj, reduction="harmony",dims=1:10)

# Clustering
seurat_obj <- FindNeighbors(seurat_obj , dims = 1:10)
seurat_obj <- FindClusters(seurat_obj , resolution = 0.5)

# Subset
SetIdent(seurat_obj,value="seurat_clusters")
UMAPPlot(seurat_obj)
seurat_A <- subset(seurat_obj,idents = 0)
```

We chose cells of type A, and conducted trajectory inference via `slingshot`
```R
# Trajectory inference
library(SingleCellExperiment)
library(slingshot)
sce <- as.SingleCellExperiment(seurat_A)
sce <- slingshot(sce,reducedDim="UMAP")
df <- sce@colData
df$slingshot <- NULL
df <- as.data.table(df,keep.rownames = T) 
df %>% fwrite("ms_ti.csv")
```

Get candidate genes via differential expression analysis
```R
# DE
seurat_A <- SetIdent(seurat_A,value="disease")
DefaultAssay(seurat_A) <- "RNA" 
markers <- FindMarkers(seurat_A, ident.1 = 0, ident.2 = 1)
markers %>% rownames_to_column(var = "gene") %>% as_tibble() %>% write_csv("MS_markers.csv")
```

Caculate the cumulative expression effects via PACE
```R
# first
```

Prepare data for eQTL. For this simulated dataset, we do not specify the locations of genes and SNPs.
```R
library(data.table)
library(tidyverse)
cum_mat <- fread("pace_cum.csv")
```

Conduct eQTL mapping
