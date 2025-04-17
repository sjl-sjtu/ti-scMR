library(tidyverse)
library(data.table)
library(Seurat)
library(GIFT)
library(ctwas)

ge <- fread("sim_sc.csv")
out <- fread("../phenotypes_logit.csv")
outcomes <- paste0("Y",1:4)
label <- read_csv("../label.csv")
out <- ge %>% dplyr::select(pseudotime,state,id,cellid) %>%
  right_join(out %>% dplyr::select(id,all_of(outcomes)), by="id") %>% as.data.frame()
rownames(out) <- out$cellid
seurat_obj <- CreateSeuratObject(counts=ge%>%dplyr::select(starts_with("gene"))%>%t(),meta.data = out)
seurat_obj <- NormalizeData(seurat_obj)

res1 <- FindMarkers(seurat_obj,ident.1 = 0,
            ident.2 = 1,
            group.by = "Y1", logfc.threshold = 0)
gene1 <- res1 %>% filter(p_val_adj<0.05) %>% rownames()
res1$label <- label$Y1
res2 <- FindMarkers(seurat_obj,ident.1 = 0,
                    ident.2 = 1,
                    group.by = "Y2", logfc.threshold = 0)
res2$label <- label$Y2
gene2 <- res2 %>% filter(p_val_adj<0.05) %>% rownames()
res3 <- FindMarkers(seurat_obj,ident.1 = 0,
                    ident.2 = 1,
                    group.by = "Y3", logfc.threshold = 0)
gene3 <- res3 %>% filter(p_val_adj<0.05) %>% rownames()
res3$label <- label$Y3
res4 <- FindMarkers(seurat_obj,ident.1 = 0,
                    ident.2 = 1,
                    group.by = "Y4", logfc.threshold = 0)
gene4 <- res4 %>% filter(p_val_adj<0.05) %>% rownames()
res4$label <- label$Y4

FDR1 <- nrow(res1%>%filter(p_val_adj<0.05)%>%filter(label==0))/nrow(res1%>%filter(p_val_adj<0.05))
power1 <- nrow(res1%>%filter(label==1)%>%filter(p_val_adj<0.05))/nrow(res1%>%filter(label==1))
FDR2 <- nrow(res2%>%filter(p_val_adj<0.05)%>%filter(label==0))/nrow(res2%>%filter(p_val_adj<0.05))
power2 <- nrow(res2%>%filter(label==1)%>%filter(p_val_adj<0.05))/nrow(res2%>%filter(label==1))
FDR3 <- nrow(res3%>%filter(p_val_adj<0.05)%>%filter(label==0))/nrow(res3%>%filter(p_val_adj<0.05))
power3 <- nrow(res3%>%filter(label==1)%>%filter(p_val_adj<0.05))/nrow(res3%>%filter(label==1))
FDR4 <- nrow(res4%>%filter(p_val_adj<0.05)%>%filter(label==0))/nrow(res4%>%filter(p_val_adj<0.05))
power4 <- nrow(res4%>%filter(label==1)%>%filter(p_val_adj<0.05))/nrow(res4%>%filter(label==1))

tibble(FDR=c(FDR1,FDR2,FDR3,FDR4),power=c(power1,power2,power3,power4)) %>% write_csv("DE_results_binary.csv")

###################
ge <- fread("sim_sc_scale.csv")
snps <- fread("../geno.csv")
X <- ge %>% group_by(id) %>% 
  summarise(across(starts_with("gene"), mean, na.rm = TRUE)) %>%
  dplyr::select(starts_with("gene")) %>%
  as.data.frame() %>%
  as.matrix()
gene <- colnames(X)

out <- fread("../phenotypes_quant.csv")
# dfiv <- fread("../dfiv.csv")
label <- read_csv("../label.csv")

outcomes <- paste0("Y",1:4)
# Y <- out %>% dplyr::select(all_of(outcomes)) %>% as.matrix()

p <- snps %>% dplyr::select(starts_with("snp")) %>% ncol()

eqtl <- fread("eqtl_avg.csv")[pvalue<0.01,]
snplist <- sapply(gene,function(g){
  eqtl %>% filter(gene==g)%>%pull(snps)
})
# pindex <- rep(p,length(gene))
pindex <- sapply(snplist,length)

Z <- snps %>% dplyr::select(starts_with("snp")) %>%
  mutate(across(everything(), ~ (.-mean(., na.rm = TRUE)) / sd(., na.rm = TRUE))) %>%
  as.data.frame() %>% as.matrix() #%>% as("dgCMatrix")


Z <- Z[, unlist(snplist)]


Y <- out$Y1
result <- GIFT_individual(X, Y, Z, Z, gene, pindex, 
                          maxiter=100, tol=1e-3, pleio=0, 
                          ncores=4, filter=T)
result

Y <- out$Y1
res1 <- GIFT_individual(X, Y, Z, Z, gene, pindex, 
                          maxiter=100, tol=1e-3, pleio=0, 
                          ncores=4, filter=T)
res1$label <- label$Y1
Y <- out$Y2
res2 <- GIFT_individual(X, Y, Z, Z, gene, pindex, 
                        maxiter=100, tol=1e-3, pleio=0, 
                        ncores=4, filter=T)
res2$label <- label$Y2
Y <- out$Y3
res3 <- GIFT_individual(X, Y, Z, Z, gene, pindex, 
                        maxiter=100, tol=1e-3, pleio=0, 
                        ncores=4, filter=T)
res3$label <- label$Y3
Y <- out$Y4
res4 <- GIFT_individual(X, Y, Z, Z, gene, pindex, 
                        maxiter=100, tol=1e-3, pleio=0, 
                        ncores=4, filter=T)
res4$label <- label$Y4

FDR1 <- nrow(res1%>%filter(p<0.05)%>%filter(label==0))/nrow(res1%>%filter(p<0.05))
power1 <- nrow(res1%>%filter(label==1)%>%filter(p<0.05))/nrow(res1%>%filter(label==1))
FDR2 <- nrow(res2%>%filter(p<0.05)%>%filter(label==0))/nrow(res2%>%filter(p<0.05))
power2 <- nrow(res2%>%filter(label==1)%>%filter(p<0.05))/nrow(res2%>%filter(label==1))
FDR3 <- nrow(res3%>%filter(p<0.05)%>%filter(label==0))/nrow(res3%>%filter(p<0.05))
power3 <- nrow(res3%>%filter(label==1)%>%filter(p<0.05))/nrow(res3%>%filter(label==1))
FDR4 <- nrow(res4%>%filter(p<0.05)%>%filter(label==0))/nrow(res4%>%filter(p<0.05))
power4 <- nrow(res4%>%filter(label==1)%>%filter(p<0.05))/nrow(res4%>%filter(label==1))

tibble(FDR=c(FDR1,FDR2,FDR3,FDR4),power=c(power1,power2,power3,power4)) %>% write_csv("GIFT_results_binary.csv")


# eqtl <- fread("eqtl_pace.csv")[FDR<0.01,]
# timeLength <- 20
# g <- 500
# genelist <- paste0("gene",1:g)
# outcomes <- paste0("Y",1:4)
# exposureDat <- read_csv("exposureDat.csv")
# cum_mat <- read_csv("pace_cum.csv")
# exposureDat <- exposureDat %>% left_join(out%>%dplyr::select(id,any_of(outcomes)),by="id")
# cum_mat <- cum_mat %>% left_join(out%>%dplyr::select(id,any_of(outcomes)),by="id")

