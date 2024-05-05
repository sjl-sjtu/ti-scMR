setwd("/lustre/home/acct-clsyzs/clsyzs/SunJianle/singleCellMR/Roche")
library(Seurat)
library(data.table)
library(tidyverse)
library(gridExtra)

count <- Matrix::readMM("ms_lesson_raw/ms_lesions_snRNAseq_cleaned_counts_matrix_2023-09-12.mtx")
genes <- fread("ms_lesson_raw/genes.txt")
barcodes <- read.csv("ms_lesson_raw/barcodes.txt")
rownames(barcodes) <- barcodes$cell_id
colnames(count) <- barcodes$cell_id
rownames(count) <- genes$ensembl

seurat_obj <- CreateSeuratObject(counts=count,meta.data=barcodes)
barcodes

info <- barcodes[,.(individual_id_anon,sex)]
setnames(info, "individual_id_anon", "#IID")
setnames(info, "sex", "SEX")
info <- unique(info)
info |> fwrite("EGAF00006717330/sexinfo.txt",sep=" ")

seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-",assay = 'RNA') #线粒体信息
seurat_obj[["percent.rb"]] <- PercentageFeatureSet(seurat_obj, pattern = "^RP[SL]",assay = 'RNA') #核糖体信息
VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), ncol = 4)
plot1 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "percent.rb")
plot4 <- FeatureScatter(seurat_obj, feature1 = "percent.rb", feature2 = "percent.mt")
plot1 + plot2 + plot3 + plot4
seurat_obj$log10GenesPerUMI <- log10(seurat_obj$nFeature_RNA) / log10(seurat_obj$nCount_RNA)
seurat_obj@meta.data %>%
  ggplot(aes(x = log10GenesPerUMI, color = sample_id_anon, fill = sample_id_anon)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8) +
  guides(fill=FALSE,color=FALSE)
#quality control
seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & percent.mt < 5 & percent.rb < 5)  #nUMI = nCount_RNA,nGene = nFeature_RNA

#normalization
seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize")
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 3000)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(seurat_obj), 10)
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(seurat_obj)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1/plot2
plot2

seurat_obj <- ScaleData(seurat_obj)

seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(seurat_obj))
ElbowPlot(seurat_obj)
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:17, reduction = "pca")
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)

seurat_obj <- RunUMAP(seurat_obj, dims = 1:17, n.components = 2)

seurat_obj <- SetIdent(seurat_obj,value = "type_broad")
p1 <- DimPlot(seurat_obj, reduction = "umap", pt.size = 1)
seurat_obj <- SetIdent(seurat_obj,value = "type_fine")
p2 <- DimPlot(seurat_obj, reduction = "umap", pt.size = 1)+NoLegend()
p1 / p2
seurat_obj <- SetIdent(seurat_obj,value = "diagnosis")
p3 <- DimPlot(seurat_obj, reduction = "umap", pt.size = 1)
p1+p3

seurat_obj <- SetIdent(seurat_obj,value = "diagnosis")
umap_plots <- list()
types <- unique(seurat_obj$type_broad)
disease <- unique(seurat_obj$diagnosis)
for (i in disease) {
  subset_cells <- subset(seurat_obj,idents = i)
  subset_cells <- SetIdent(subset_cells,value = "type_broad")
  if(i==disease[4]){
    umap_plot <- DimPlot(subset_cells, reduction = "umap", pt.size = 0.1,order=types)+labs(title=i)
  }else{
    umap_plot <- DimPlot(subset_cells, reduction = "umap", pt.size = 1,order=types)+labs(title=i)
  }
  umap_plots[[i]] <- umap_plot
}
combined_plots <- grid.arrange(grobs = umap_plots, nrow = 2, ncol = 2)
print(combined_plots)

saveRDS(seurat_obj,"Roche_seurat.rds")



#####
seurat_obj <- readRDS("Roche_seurat.rds")
unique(seurat_obj$type_broad)
seurat_ctr_Oligo <- subset(seurat_obj, subset = type_broad == "Oligodendrocytes" & diagnosis == "CTR")
seurat_ms_Oligo <- subset(seurat_obj, subset = type_broad == "Oligodendrocytes" & diagnosis %in% c("SPMS", "PPMS", "RRMS"))
saveRDS(seurat_ctr_Oligo,"seurat_ctr_Oligo.rds")
saveRDS(seurat_ms_Oligo,"seurat_ms_Oligo.rds")


#DE
seurat_Oligo <- subset(seurat_obj, subset = type_broad == "Oligodendrocytes")
seurat_Oligo@meta.data <- seurat_Oligo@meta.data %>% mutate(outcome = dplyr::recode(seurat_Oligo$diagnosis,  "CTR" = "Control", "SPMS"="MS", "PPMS"="MS", "RRMS"="MS")) 
seurat_Oligo <- SetIdent(seurat_Oligo,value="outcome")
DefaultAssay(seurat_Oligo) <- "RNA" 
markers <- FindMarkers(seurat_Oligo, ident.1 = "Control", ident.2 = "MS")
markers %>% rownames_to_column(var = "gene") %>% as_tibble() %>% write_csv("Oligo_Ctr_MS_markers2.csv")

unique(seurat_obj$type_broad)
seurat_ctr_Oligo <- subset(seurat_obj, subset = type_broad == "Oligodendrocytes" & diagnosis == "CTR")
seurat_ms_Oligo <- subset(seurat_obj, subset = type_broad == "Oligodendrocytes" & diagnosis %in% c("SPMS", "PPMS", "RRMS"))
saveRDS(seurat_ctr_Oligo,"seurat_ctr_Oligo.rds")
saveRDS(seurat_ms_Oligo,"seurat_ms_Oligo.rds")

####
library(harmony)
seurat_obj <- RunHarmony(seurat_obj,"sample_id_anon")
seurat_obj <- seurat_obj %>% RunUMAP(reduction="harmony",dims=1:20) %>%
  FindNeighbors(reduction="harmony",dims=1:20)%>%
  FindClusters(resolution = 0.5)

seurat_obj <- SetIdent(seurat_obj,value = "type_broad")
p1 <- DimPlot(seurat_obj, reduction = "umap", pt.size = 1)
seurat_obj <- SetIdent(seurat_obj,value = "type_fine")
p2 <- DimPlot(seurat_obj, reduction = "umap", pt.size = 1)+NoLegend()
p1 / p2
seurat_obj <- SetIdent(seurat_obj,value = "diagnosis")
p3 <- DimPlot(seurat_obj, reduction = "umap", pt.size = 1)
p1+p3

seurat_obj <- SetIdent(seurat_obj,value = "diagnosis")
umap_plots <- list()
types <- unique(seurat_obj$type_broad)
disease <- unique(seurat_obj$diagnosis)
for (i in disease) {
  subset_cells <- subset(seurat_obj,idents = i)
  subset_cells <- SetIdent(subset_cells,value = "type_broad")
  if(i==disease[4]){
    umap_plot <- DimPlot(subset_cells, reduction = "umap", pt.size = 0.1,order=types)+labs(title=i)
  }else{
    umap_plot <- DimPlot(subset_cells, reduction = "umap", pt.size = 1,order=types)+labs(title=i)
  }
  umap_plots[[i]] <- umap_plot
}
combined_plots <- grid.arrange(grobs = umap_plots, nrow = 2, ncol = 2)
print(combined_plots)
saveRDS(seurat_obj,"Roche_seurat_batch_corr.rds")
seurat_ctr_Oligo <- subset(seurat_obj, subset = type_broad == "Oligodendrocytes" & diagnosis == "CTR")
seurat_ms_Oligo <- subset(seurat_obj, subset = type_broad == "Oligodendrocytes" & diagnosis %in% c("SPMS", "PPMS", "RRMS"))
saveRDS(seurat_ctr_Oligo,"seurat_ctr_Oligo_batch.rds")
saveRDS(seurat_ms_Oligo,"seurat_ms_Oligo_batch.rds")

###########
## use after hamony
seurat_ms_Oligo <- readRDS("seurat_ms_Oligo_batch.rds")
seurat_ctr_Oligo <- readRDS("seurat_ctr_Oligo_batch.rds")

library(SingleCellExperiment)
library(slingshot)
library(scater)
library(scran)
library(bluster)

# reconst_data <- readRDS("scvi_Oligo_reconst.rds")
# seurat_oligo < readRDS("seurat_Oligo_scvi.rds")

sce <- as.SingleCellExperiment(seurat_ms_Oligo)
SetIdent(seurat_ms_Oligo,value = "type_fine")
DimPlot(seurat_ms_Oligo, reduction = "umap", pt.size = 1,group.by="type_fine")


plotUMAP(sce, colour_by="type_fine")
colnames(sce)|>length()

outlier <- which((sce@int_colData@listData$reducedDims@listData$UMAP[,2] < 1) | 
                   (sce@int_colData@listData$reducedDims@listData$UMAP[,1] < -4) #| 
                   #((sce@int_colData@listData$reducedDims@listData$UMAP[,1] > -3) & (sce@int_colData@listData$reducedDims@listData$UMAP[,2] < 2)) 
) |>
  names()
sce <- sce[,!colnames(sce) %in% outlier]
plotUMAP(sce, colour_by="type_fine")

# library(scater)
# sce <- logNormCounts(sce)
library(scran)
dec <- modelGeneVar(sce)
hvg <- getTopHVGs(dec, prop=0.1)
# library(scater)
# sce <- runPCA(sce, ncomponents=15, subset_row=hvg)
# library(bluster)
# colLabels(sce) <- clusterCells(sce, use.dimred='PCA',
#                                BLUSPARAM=NNGraphParam(cluster.fun="louvain")) 
# sce <- runUMAP(sce, dimred = 'PCA')
p1 <- plotUMAP(sce, colour_by="label")
p2 <- plotUMAP(sce, colour_by="type_fine")
plotUMAP(sce, colour_by="diagnosis")
plotUMAP(sce, colour_by="individual_id_anon")
p1+p2

# sce <- slingshot(sce,sce$type_fine,reducedDim="UMAP")
# sce <- slingshot(sce,reducedDim="UMAP")
sce <- slingshot(sce,reducedDim="UMAP")
p1 <- function(){par(mar = c(1,1,1,1));plot.new();
  plot(sce@int_colData@listData$reducedDims@listData$UMAP[,1:2],#seurat_obj1@reductions$umap@cell.embeddings[, 1:2], 
       col= as.factor(sce$type_fine),#seurat_obj1$predicted.celltype.l2,
       pch = 16,cex=0.5,xlab="UMAP_1",ylab="UMAP_2",main="MS");
  lines(SlingshotDataSet(sce)@curves$Lineage1, lwd=2, col='black');
}
p2 <- function(){par(mar = c(1,1,1,1));plot.new();
  legend("top", legend = unique(as.factor(sce$type_fine)), col = unique(as.factor(sce$type_fine)),
         pch = 16, title = "cell type",border = NA)}
pp2 <- cowplot::plot_grid(p1, p2, align="v", ncol=2,rel_widths = c(3, 1))
pp2
cowplot::save_plot("ms_Oligo_batch.png",pp2,dpi=600,base_height=6,base_width=14)
saveRDS(sce,"ms_Oligo_batch.rds")

df <- sce@colData
df$slingshot <- NULL
aggregate(slingPseudotime_1 ~ type_fine, data = df, FUN = mean)
# aggregate(slingPseudotime_2 ~ cell_type, data = df, FUN = mean)
df <- as.data.table(df,keep.rownames = T) 
df %>% fwrite("ms_Oligo_batch.csv")


#
sce <- as.SingleCellExperiment(seurat_ctr_Oligo)
plotUMAP(sce, colour_by="type_fine")

outlier <- which((sce@int_colData@listData$reducedDims@listData$UMAP[,2] < 1) | 
                   (sce@int_colData@listData$reducedDims@listData$UMAP[,1] < -4) #| 
                 #((sce@int_colData@listData$reducedDims@listData$UMAP[,1] > -3) & (sce@int_colData@listData$reducedDims@listData$UMAP[,2] < 2)) 
) |>
  names()
sce <- sce[,!colnames(sce) %in% outlier]
plotUMAP(sce, colour_by="type_fine")

# sce <- logNormCounts(sce)
dec <- modelGeneVar(sce)
hvg <- getTopHVGs(dec, prop=0.1)
# sce <- runPCA(sce, ncomponents=15, subset_row=hvg)
# colLabels(sce) <- clusterCells(sce, use.dimred='PCA',
#                                BLUSPARAM=NNGraphParam(cluster.fun="louvain")) 
# sce <- runUMAP(sce, dimred = 'PCA')
plotUMAP(sce, colour_by="label")

plotUMAP(sce, colour_by="sample_id_anon")
plotUMAP(sce, colour_by="individual_id_anon")

sce <- slingshot(sce,reducedDim="UMAP")
p1 <- function(){par(mar = c(1,1,1,1));plot.new();
  plot(sce@int_colData@listData$reducedDims@listData$UMAP[,1:2],#seurat_obj1@reductions$umap@cell.embeddings[, 1:2], 
       col= as.factor(sce$type_fine),#seurat_obj1$predicted.celltype.l2,
       pch = 16,cex=0.5,xlab="UMAP_1",ylab="UMAP_2",main="Control");
  lines(SlingshotDataSet(sce)@curves$Lineage1, lwd=2, col='black');
}
p2 <- function(){par(mar = c(1,1,1,1));plot.new();
  legend("top", legend = unique(as.factor(sce$type_fine)), col = unique(as.factor(sce$type_fine)),
         pch = 16, title = "cell type",border = NA)}
pp2 <- cowplot::plot_grid(p1, p2, align="v", ncol=2,rel_widths = c(3, 1))
pp2
cowplot::save_plot("ctr_Oligo_batch.png",pp2,dpi=600,base_height=6,base_width=14)
saveRDS(sce,"ctr_Oligo_batch.rds")

df <- sce@colData
df$slingshot <- NULL
aggregate(slingPseudotime_1 ~ type_fine, data = df, FUN = mean)
# aggregate(slingPseudotime_2 ~ cell_type, data = df, FUN = mean)
df <- as.data.table(df,keep.rownames = T) 
df %>% fwrite("ctr_Oligo_batch.csv")

length(unique(sce$sample_id_anon))
length(unique(sce$individual_id_anon))
