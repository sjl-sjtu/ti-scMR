# ti-scMR

Trajectory-inference-based dynamic single-cell Mendelian randomization (ti-scMR)

![](https://github.com/sjl-sjtu/ti-scMR/blob/main/Figs/ti-scMR.jpg)

![](https://github.com/sjl-sjtu/ti-scMR/blob/main/Figs/flowchart.jpg)

## Tutorial
Here we used a simulated toy dataset to illustrate the workflow of ti-scMR. Before starting the workflow, you need to prepare at least a single cell object (e.g. Seurat object `seurat_sc.rds` or h5ad file for Scanpy) and a file containing all genotypes `genotypes.csv` (It can be processed from VCF or other formats via PLINK, refer `application/onek1k/genotype.sh`). We provided a toy dataset on `example` folder. For real analysis, you should prepare other files like the locations of each SNP and gene, as well as covariates information.

We first process the sc-RNA count matrix using `Seurat` v5 (<https://satijalab.org/seurat/>).
```R
library(data.table)
library(tidyverse)
library(Seurat)

seurat_obj <- readRDS("seurat_sc.rds")

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

# Batch correction (if necessary)
# library(harmony)
# seurat_obj <- RunHarmony(seurat_obj,"batch")

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

We chose cells of type A, and conducted trajectory inference via `slingshot` (<https://nbisweden.github.io/workshop-archive/workshop-scRNAseq/2020-01-27/labs/compiled/slingshot/slingshot.html>)
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

Caculate the cumulative expression effects via PACE. We provide a function `cum_expression` (<https://github.com/sjl-sjtu/ti-scMR/blob/main/R/cum_expression.R>) to calculate cumulative effects of all candidate genes, which is implemented based on R package `fdapace` (<https://cran.r-project.org/web/packages/fdapace/vignettes/fdapaceVig.html>).
```R
# first extract the processed RNA abundance levels
df_rna <- seurat_A@assays$RNA$scale.data %>% t() %>% as_tibble()
df_rna$cellid <- colnames(seurat_A)

# combine with cell information
df_rna <- df_rna %>% left_join(seurat_A@meta.data,by="cellid") %>%
  dplyr::select(id,cellid,any_of(genelist),disease) %>% 
  left_join(dfti[,.(cellid,slingPseudotime_1)],by="cellid")

# transfer pseudotime points to time period (by rounding the pseudotime)
df_rna <- df_rna %>% mutate(pseudotime=round(slingPseudotime_1))
df_rna %>% write_csv("express.csv")

# restore the cumulative effects
source("cum_expression.R")
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

Conduct eQTL mapping. If you have location information for genes and SNPs, it is recommended to conduct *cis*- and *trans*-eQTL mapping separately. Please refer `MatrixEQTL` (<https://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/runit.html>) for more tutorials.
```R
# using MatrixEQTL tools
library(MatrixEQTL)

base.dir = getwd()
useModel = modelLINEAR 
SNP_file_name = paste(base.dir, "/SNP.txt", sep="") 
expression_file_name = paste(base.dir, "/GE_pace.txt", sep="") 
covariates_file_name = paste(base.dir, "/Covariates.txt", sep="")
output_file_name = tempfile() 
pvOutputThreshold = 1e-4 
errorCovariance = numeric() 

snps = SlicedData$new() 
snps$fileDelimiter = " "       
snps$fileOmitCharacters = "NA" 
snps$fileSkipRows = 1        
snps$fileSkipColumns = 1     
snps$fileSliceSize = 2000      
snps$LoadFile( SNP_file_name ) 

gene = SlicedData$new()
gene$fileDelimiter = " " 
gene$fileOmitCharacters = "NA"
gene$fileSkipRows = 1
gene$fileSkipColumns = 1
gene$fileSliceSize = 2000
gene$LoadFile(expression_file_name)

cvrt = SlicedData$new();
cvrt$fileDelimiter = " "; 
cvrt$fileOmitCharacters = "NA";
cvrt$fileSkipRows = 1;  
cvrt$fileSkipColumns = 1; 
if(length(covariates_file_name)>0) {
  cvrt$LoadFile(covariates_file_name);
}


maf.list = vector('list', length(snps))
for(sl in 1:length(snps)) {
  slice = snps[[sl]];
  maf.list[[sl]] = rowMeans(slice,na.rm=TRUE)/2;
  maf.list[[sl]] = pmin(maf.list[[sl]],1-maf.list[[sl]]);
}
maf = unlist(maf.list)
snps$RowReorder(maf>=0.01);

me = Matrix_eQTL_engine(
  snps = snps,
  gene = gene,
  cvrt = cvrt,
  output_file_name = output_file_name, 
  pvOutputThreshold = pvOutputThreshold,
  useModel = useModel,
  errorCovariance = errorCovariance,
  verbose = TRUE,
  min.pv.by.genesnp = FALSE,  pvalue.hist = "qqplot",
  noFDRsaveMemory = FALSE)

eqtl <- me$all$eqtls
eqtl$gene |> unique() %>% length()
eqtl %>% write_csv("eqtl_pace.csv")
```

We finally conducted MR analysis. We provide an R function `sc_mr` (<https://github.com/sjl-sjtu/ti-scMR/blob/main/R/sc_mr.R>) to conduct single-cell Mendelian randomization.
```R
library(doParallel)
library(foreach)
registerDoParallel(40)
getDoParWorkers()

library(fdapace)
library(tidyverse)
library(data.table)
library(pracma)
library(sandwich)
library(aod)
library(ivreg)
library(glmnet)
library(pROC)
library(boot)
library(boot.pval)

source("sc_mr.R")

loaded_packages <- search()[grepl("package:", search())]
loaded_packages <- sub("package:", "", loaded_packages)

snps <- fread("genotypes.csv")
cum_mat <- fread("pace_cum.csv")
eqtl <- fread("eqtl_pace.csv")[FDR<0.05,]

# LD pruning of IVs
stepPrune <- function(df,dfgwas,cutoff){
  snplist <- dfgwas %>% pull(snps)
  candidate <- c()
  while(length(snplist)>0){
    dfgwas <- dfgwas%>%filter(snps %in% snplist)
    j <- dfgwas[which.min(dfgwas$FDR),"snps"] %>% as.character()
    candidate <- c(candidate,j)
    snplist <- setdiff(snplist,j)
    calcLD <- function(i,j,df){
      return(cor(df%>%pull(i),df%>%pull(j))^2)  # Rogers and Huff (2008)
    }
    LDs <- sapply(snplist,calcLD,j,df) 
    snplist <- snplist[which(LDs<cutoff)]
  }
  return(candidate)
}

res <- foreach(geneName = genelist,.packages = loaded_packages) %dopar% {
      outcome <- "disease"
      IV <- eqtl %>% filter(gene==geneName) %>% filter(FDR<0.05) %>% pull(snps)
      dfex <- cum_mat %>% left_join(snps %>% select_at(vars(any_of(IV))), by="id")
      if(length(IV)>0){
        dfgwas <- eqtl %>% filter(gene==geneName) %>% filter(snps %in% IV) %>% select(snps,FDR)
        IV <- stepPrune(dfex,dfgwas,0.1)
      }
      return(sc_mr(dfex,geneName,outcome,IV,method="logit_lasso",id_var="id",bootstraps=5000))
    }
res <- reduce(res,rbind)
colnames(res) <- c("IVnum","F","beta","ci_l","ci_u","se","p")
res <- as_tibble(res)
res$gene <- genelist
res$padj <- p.adjust(res$p,method = "BH")
res %>% write_csv("results.csv")
```

## Citation
Sun, J., Dong, Q., Wei, J., Gao, Y., Yu, Z., Hu, X., & Zhang, Y. (2025). ti-scMR: trajectory-inference-based dynamic single-cell Mendelian randomization identifies causal genes underlying phenotypic differences. NAR Genomics and Bioinformatics, 7(3), lqaf082

## Contact
If you have any questions, contact me via <jianles@andrew.cmu.edu> or <sjl-2017@sjtu.alumni.edu.cn>.
