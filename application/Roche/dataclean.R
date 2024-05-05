setwd("/lustre/home/acct-clsyzs/clsyzs/SunJianle/singleCellMR/Roche")
library(Seurat)
library(data.table)
library(tidyverse)
library(gridExtra)

count <- Matrix::readMM("ms_lesson_raw/ms_lesions_snRNAseq_cleaned_counts_matrix_2023-09-12.mtx")
genes <- fread("ms_lesson_raw/genes.txt")
barcodes <- read.csv("ms_lesson_raw/barcodes.txt")
colnames(count) <- barcodes$cell_id
rownames(count) <- genes$ensembl
rownames(barcodes) <- barcodes$cell_id
seurat_obj <- CreateSeuratObject(counts=count, meta.data=barcodes)
barcodes

# info <- barcodes[,.(individual_id_anon,sex)]
# setnames(info, "individual_id_anon", "#IID")
# setnames(info, "sex", "SEX")
# info <- unique(info)
# info |> fwrite("EGAF00006717330/sexinfo.txt",sep=" ")

seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-",assay = 'RNA') #线粒体信息
seurat_obj[["percent.rb"]] <- PercentageFeatureSet(seurat_obj, pattern = "^RP[SL]",assay = 'RNA') #核糖体信息
# VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), ncol = 4)
# plot1 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
# plot2 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# plot3 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "percent.rb")
# plot4 <- FeatureScatter(seurat_obj, feature1 = "percent.rb", feature2 = "percent.mt")
# plot1 + plot2 + plot3 + plot4
seurat_obj$log10GenesPerUMI <- log10(seurat_obj$nFeature_RNA) / log10(seurat_obj$nCount_RNA)
# seurat_obj@meta.data %>%
#   ggplot(aes(x = log10GenesPerUMI, color = sample_id_anon, fill = sample_id_anon)) +
#   geom_density(alpha = 0.2) +
#   theme_classic() +
#   geom_vline(xintercept = 0.8) +
#   guides(fill=FALSE,color=FALSE)
#quality control
seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & percent.mt < 5 & percent.rb < 5)  #nUMI = nCount_RNA,nGene = nFeature_RNA

saveRDS(seurat_obj,file="seurat_original.rds")
