# ti-scMR

Trajectory-inference-based dynamic single-cell Mendelian randomization (ti-scMR)

## Tutorial
Here we used a simulated toy dataset to illustrate the workflow of ti-scMR. We first process the sc-RNA count matrix using Seurat v5.
```R
library(Seurat)
seurat_obj <- readRDS("example_sc.rds")
```
