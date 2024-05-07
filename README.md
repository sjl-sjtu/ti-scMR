# ti-scMR

Trajectory-inference-based dynamic single-cell Mendelian randomization (ti-scMR)

## Tutorial
Here we used a simulated toy dataset to illustrate the workflow of ti-scMR. We first process the sc-RNA count matrix using `Seurat` v5.
```R
library(Seurat)
seurat_obj <- readRDS("example_sc.rds")

# Quality Control
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-",assay = 'RNA')
seurat_obj[["percent.rb"]] <- PercentageFeatureSet(seurat_obj, pattern = "^RP[SL]",assay = 'RNA')
seurat_obj$log10GenesPerUMI <- log10(seurat_obj$nFeature_RNA) / log10(seurat_obj$nCount_RNA)
seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & percent.mt < 5 & percent.rb < 5)

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

```

We chose cells of type A, and conducted trajectory inference via `slingshot`
```R
# Trajectory inference
library(SingleCellExperiment)
library(slingshot)
library(scater)
library(scran)
library(bluster)
sce <- as.SingleCellExperiment(seurat_A)
sce <- slingshot(sce,reducedDim="UMAP")
```

Get candidate genes via differential expression analysis
```R
# DE
seurat_Oligo <- SetIdent(seurat_Oligo,value="outcome")
DefaultAssay(seurat_Oligo) <- "RNA" 
markers <- FindMarkers(seurat_Oligo, ident.1 = "Control", ident.2 = "MS")
markers %>% rownames_to_column(var = "gene") %>% as_tibble() %>% write_csv("Oligo_Ctr_MS_markers.csv")
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
