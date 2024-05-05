setwd("/lustre/home/acct-clsyzs/clsyzs/SunJianle/singleCellMR/GSE196829/expression_matrix")

library(Seurat)
library(data.table)
library(tidyverse)

seurat_obj <- readRDS("../onek1k.rds")
seurat_obj
# seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200  & nFeature_RNA < 2500)
# seurat_obj

seurat_obj <- SetIdent(seurat_obj,value = "predicted.celltype.l2")
DimPlot(seurat_obj, reduction = "umap", pt.size = 1)

seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
seurat_obj <- ScaleData(seurat_obj)

library(SingleCellExperiment)
library(slingshot)
library(scater)
library(scran)
library(bluster)
seurat_obj <- SetIdent(seurat_obj,value = "predicted.celltype.l2")
seurat_obj2 <- subset(seurat_obj, idents = c("B naive","B intermediate","Plasmablast"))

# seurat_obj2 <- SetIdent(seurat_obj2,value = "donor_id")
# DimPlot(seurat_obj2, reduction = "umap", pt.size = 1)+NoLegend()

sce <- as.SingleCellExperiment(seurat_obj2)

plotUMAP(sce, colour_by="predicted.celltype.l2")
# plotUMAP(sce, colour_by="donor_id")

outlier <- which((sce@int_colData@listData$reducedDims@listData$UMAP[,2]<5) | (sce@int_colData@listData$reducedDims@listData$UMAP[,1] < -4)) |> names()
sce <- sce[,!colnames(sce) %in% outlier]
plotUMAP(sce, colour_by="predicted.celltype.l2")
# sce <- slingshot(sce)
sce <- slingshot(sce, clusterLabels = 'predicted.celltype.l2', reducedDim = 'UMAP',
                 start.clus = 'B naive', end.clus = 'Plasmablast')
sce$slingshot

p1 <- function(){par(mar = c(1,1,1,1));plot.new();
  plot(sce@int_colData@listData$reducedDims@listData$UMAP[,1:2],#seurat_obj1@reductions$umap@cell.embeddings[, 1:2], 
       col= sce$predicted.celltype.l2,#seurat_obj1$predicted.celltype.l2,
       pch = 16,cex=0.1);
  lines(SlingshotDataSet(sce)@curves$Lineage1, lwd=2, col='black');
}
p2 <- function(){par(mar = c(1,1,1,1));plot.new();
  legend("top", legend = unique(sce$predicted.celltype.l2), col = unique(sce$predicted.celltype.l2), 
         pch = 16, title = "cell type",border = NA)}
pp2 <- cowplot::plot_grid(p1, p2, align="v", ncol=2,rel_widths = c(5, 2))
pp2
cowplot::save_plot("plasa_batchcorr.png",pp2,dpi=600,base_height=6,base_width=10)
saveRDS(sce,"plasa_batchcorr.rds")

plot.new()
plot(sce@int_colData@listData$reducedDims@listData$UMAP[,1:2],#seurat_obj1@reductions$umap@cell.embeddings[, 1:2], 
     col= sce$predicted.celltype.l2,#seurat_obj1$predicted.celltype.l2,
     pch = 16,cex=0.1);
lines(SlingshotDataSet(sce)@curves$Lineage1, lwd=2, col='black')
legend("top",legend = unique(sce$predicted.celltype.l2), col = unique(sce$predicted.celltype.l2), 
       pch = 16, title = "cell type",border = NA,xpd=T,ncol=3,inset = c(0,-0.3))

df <- sce@colData
df$slingshot <- NULL
aggregate(slingPseudotime_1 ~ cell_type, data = df, FUN = mean)
df <- as.data.table(df,keep.rownames = T) 
df
df %>% fwrite("plasa_batchcorr.csv")

seurat_obj2 <- subset(seurat_obj, idents = c("B naive","B intermediate","B memory"))
sce <- as.SingleCellExperiment(seurat_obj2)
plotUMAP(sce, colour_by="predicted.celltype.l2")
outlier <- which((sce@int_colData@listData$reducedDims@listData$UMAP[,2]<5) | (sce@int_colData@listData$reducedDims@listData$UMAP[,1] < -4) | (sce@int_colData@listData$reducedDims@listData$UMAP[,1] > 4)) |> names()
sce <- sce[,!colnames(sce) %in% outlier]
plotUMAP(sce, colour_by="predicted.celltype.l2")

sce <- slingshot(sce, clusterLabels = 'predicted.celltype.l2', reducedDim = 'UMAP',
                 start.clus = 'B naive', end.clus = 'B memory')
sce$slingshot
p1 <- function(){par(mar = c(1,1,1,1));plot.new();
  plot(sce@int_colData@listData$reducedDims@listData$UMAP[,1:2],#seurat_obj1@reductions$umap@cell.embeddings[, 1:2], 
       col= sce$predicted.celltype.l2,#seurat_obj1$predicted.celltype.l2,
       pch = 16,cex=0.1);
  lines(SlingshotDataSet(sce)@curves$Lineage1, lwd=2, col='black');
}
p2 <- function(){par(mar = c(1,1,1,1));plot.new();
  legend("top", legend = unique(sce$predicted.celltype.l2), col = unique(sce$predicted.celltype.l2), 
         pch = 16, title = "cell type",border = NA)}
pp2 <- cowplot::plot_grid(p1, p2, align="v", ncol=2,rel_widths = c(5, 2))
pp2

par(mar = c(5,5,5,5), xpd=TRUE)
plot.new()
plot(sce@int_colData@listData$reducedDims@listData$UMAP[,1:2],#seurat_obj1@reductions$umap@cell.embeddings[, 1:2], 
     col= sce$predicted.celltype.l2,#seurat_obj1$predicted.celltype.l2,
     pch = 16,cex=0.1);
lines(SlingshotDataSet(sce)@curves$Lineage1, lwd=2, col='black')
legend("top",legend = unique(sce$predicted.celltype.l2), col = unique(sce$predicted.celltype.l2), 
       pch = 16, title = "cell type",border = NA,xpd=T,ncol=3,inset = c(0,-0.3))

cowplot::save_plot("Bmemory_batchcorr.png",pp2,dpi=600,base_height=6,base_width=10)
saveRDS(sce,"Bmemory_batchcorr.rds")

df <- sce@colData
df$slingshot <- NULL
aggregate(slingPseudotime_1 ~ cell_type, data = df, FUN = mean)
df <- as.data.table(df,keep.rownames = T) 
df %>% fwrite("B_memory_batchcorr.csv")

seurat_obj1 <- subset(seurat_obj, idents = c("B naive","B intermediate",'B memory',"Plasmablast"))
saveRDS(seurat_obj1,"seurat_B_bc.rds")

sce <- readRDS("plasa_batchcorr.rds")
