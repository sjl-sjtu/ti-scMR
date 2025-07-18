library(Seurat)
library(data.table)
library(tidyverse)
library(SingleCellExperiment)
library(slingshot)
library(tradeSeq)
library(scater)
library(scran)
library(bluster)

seurat_obj <- readRDS("seurat_B_bc.rds")
seurat_obj <- SetIdent(seurat_obj,value = "predicted.celltype.l2")
# seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- subset(seurat_obj,idents=c("B naive","B intermediate",
                                        "B memory","Plasmablast"))
sce <- as.SingleCellExperiment(seurat_obj)
rm(seurat_obj)
# plotUMAP(sce, colour_by="predicted.celltype.l2")
outlier <- which((sce@int_colData@listData$reducedDims@listData$UMAP[,2]<5) | (sce@int_colData@listData$reducedDims@listData$UMAP[,1] < -4)) |> names()
sce <- sce[,!colnames(sce) %in% outlier]
# plotUMAP(sce, colour_by="predicted.celltype.l2")

counts(sce) <- logcounts(sce)
sce <- logNormCounts(sce)
dec <- modelGeneVar(sce)
top_hvgs <- getTopHVGs(dec, n = 2000)

selected_cells <- colData(sce) %>% 
  as_tibble(rownames="cell_id") %>%
  group_by(predicted.celltype.l2) %>%
  sample_frac(size = 0.1) %>%
  pull(cell_id)

sce <- slingshot(sce, clusterLabels = 'predicted.celltype.l2', 
                 reducedDim = 'UMAP',
                 start.clus = 'B naive', 
                 end.clus = c("B memory",'Plasmablast'))
p1 <- function(){par(mar = c(1,1,1,1));plot.new();
  plot(sce@int_colData@listData$reducedDims@listData$UMAP[,1:2],#seurat_obj1@reductions$umap@cell.embeddings[, 1:2],
       col= sce$predicted.celltype.l2,#seurat_obj1$predicted.celltype.l2,
       pch = 16,cex=0.1);
  lines(SlingshotDataSet(sce)@curves$Lineage1, lwd=2, col='black');
  lines(SlingshotDataSet(sce)@curves$Lineage2, lwd=2, col='yellow')
}
p2 <- function(){par(mar = c(1,1,1,1));plot.new();
  legend("top", legend = unique(sce$predicted.celltype.l2), col = unique(sce$predicted.celltype.l2),
         pch = 16, title = "cell type",border = NA)}
pp2 <- cowplot::plot_grid(p1, p2, align="v", ncol=2,rel_widths = c(5, 2))
pp2

sce_hvg <- sce[top_hvgs, selected_cells]
counts <- counts(sce_hvg)
pseudotime <- slingPseudotime(sce, na = FALSE)[selected_cells,]
cellWeights <- slingCurveWeights(sce_hvg)


sce_hvg <- fitGAM(counts = counts, pseudotime = pseudotime, 
                  cellWeights = cellWeights,
                  nknots = 6, verbose = FALSE)

endRes <- diffEndTest(sce_hvg)
endRes |> as_tibble(rownames="genes") |> arrange(pvalue) |> write_csv("trajectory_difference_B.csv")

# look at significant genes
o <- order(endRes$waldStat, decreasing = TRUE)
sigGene <- names(sce_hvg)[o[1]]
plotSmoothers(sce_hvg, counts, sigGene)

