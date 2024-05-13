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
dfti <- sce@colData
dfti$slingshot <- NULL
dfti <- as.data.table(dfti,keep.rownames = T) 
dfti %>% fwrite("example_ti.csv")
```

Get candidate genes via differential expression analysis
```R
# DE
seurat_A <- SetIdent(seurat_A,value="disease")
DefaultAssay(seurat_A) <- "RNA" 
markers <- FindMarkers(seurat_A, ident.1 = 0, ident.2 = 1)
markers <- markers %>% rownames_to_column(var = "gene") %>% as_tibble()
genelist <- markers%>%filter(p_val_adj<0.05)%>%pull(gene)
markers %>% write_csv("example_markers.csv")
```

Caculate the cumulative expression effects via PACE
```R
# first extract the processed RNA abundance levels
df_rna <- seurat_A@assays$RNA$scale.data %>% t() %>% as_tibble()
df_rna$cellid <- colnames(seurat_A)
df_rna <- df_rna %>% left_join(seurat_A@meta.data,by="cellid") %>%
  dplyr::select(id,cellid,any_of(genelist),disease) %>% 
  left_join(dfti[,.(cellid,slingPseudotime_1)],by="cellid")

# transfer pseudotime points to time period (by rounding the pseudotime)
df_rna <- df_rna %>% mutate(pseudotime=round(slingPseudotime_1))
df_rna %>% write_csv("express.csv")

# restore the cumulative effects
df_cum <- cum_expression(df_rna,genelist,id_var="id",pseudotime_var="pseudotime")
df_cum %>% write_csv("pace_cum.csv")
```

Prepare data for eQTL. For this simulated dataset, we do not specify the locations of genes and SNPs. We include genomic PCs as covariates. Please check requirments of data format for package `MatrixEQTL` if you want to add these information in real data analysis.
```R
# Gene expression data
df <- fread("pace_cum.csv")
gene <- t(df[,genelist,with=F])
colnames(gene) <- df$id
gene <- as.data.table(gene)
gene$geneid <- genelist
gene <- gene[,c("geneid",df$id),with=FALSE]
gene %>% write_delim("GE_pace.txt")

# Genotype data
dat <- fread("genotypes.csv")
snplist <- setdiff(colnames(snps),"id")
snp <- t(dat[,snplist,with=FALSE])
colnames(snp) <- df$id
snp <- as.data.table(snp)
snp[,snpid:= snplist]
snp <- snp[,c("snpid",df$id),with=F]
snp %>% fwrite("SNP.txt",row.names=FALSE,sep=" ")

# Genomic PCs as covariates
fit_prcomp <- prcomp(dat[,paste0("snp",1:p),with=FALSE],center = T,scale. = T)
m <- 10
covs <- fit_prcomp$x[,paste0("PC",1:m)]
covs <- covs |> t() |> as_tibble()
colnames(covs) <- df$id
covs$id <- paste0("PC",1:m)
covs <- covs%>%select(id,df$id)
covs %>% fwrite("Covariates.txt",row.names=FALSE,sep=" ")

# We do not consider locations of SNPs and genes in this toy example.
```

Conduct eQTL mapping
```R
# using MatrixEQTL tools

eqtl %>% write_csv("eqtl_pace.csv")
```

We finally conducted MR analysis
```R
# IV

```
