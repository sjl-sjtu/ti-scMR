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

out <- fread("../phenotypes_logit.csv")
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

library(tidyverse)
res1 <- read_csv("gift_y1.csv")
res1$p_adj <- p.adjust(res1$p, method = "BH")
FDR1 <- nrow(res1%>%filter(p_adj<0.05)%>%filter(label==0))/nrow(res1%>%filter(p_adj<0.05))
power1 <- nrow(res1%>%filter(label==1)%>%filter(p_adj<0.05))/nrow(res1%>%filter(label==1))

res2 <- read_csv("gift_y2.csv")
res2$p_adj <- p.adjust(res2$p, method = "BH")
FDR2 <- nrow(res2%>%filter(p_adj<0.05)%>%filter(label==0))/nrow(res2%>%filter(p_adj<0.05))
power2 <- nrow(res2%>%filter(label==1)%>%filter(p_adj<0.05))/nrow(res2%>%filter(label==1))

res3 <- read_csv("gift_y3.csv")
res3$p_adj <- p.adjust(res3$p, method = "BH")
FDR3 <- nrow(res3%>%filter(p_adj<0.05)%>%filter(label==0))/nrow(res3%>%filter(p_adj<0.05))
power3 <- nrow(res3%>%filter(label==1)%>%filter(p_adj<0.05))/nrow(res3%>%filter(label==1))

res4 <- read_csv("gift_y4.csv")
res4$p_adj <- p.adjust(res4$p, method = "BH")
FDR4 <- nrow(res4%>%filter(p_adj<0.05)%>%filter(label==0))/nrow(res4%>%filter(p_adj<0.05))
power4 <- nrow(res4%>%filter(label==1)%>%filter(p_adj<0.05))/nrow(res4%>%filter(label==1))

gift <- tibble(FDR=c(FDR1,FDR2,FDR3,FDR4),power=c(power1,power2,power3,power4))
de <- read_csv("DE_results_binary.csv")
tisc1 <- read_csv("binary_rep_mixed_FDR.csv")
tisc1 <- tisc1 %>% group_by(outcome) %>% summarise_all(mean,na.rm=T) %>% 
  dplyr::select(outcome,avg_linear,pace_linear_lasso)
tisc2 <- read_csv("binary_rep_mixed_power.csv")
tisc2 <- tisc2 %>% group_by(outcome) %>% summarise_all(mean,na.rm=T) %>% 
  dplyr::select(outcome,avg_linear,pace_linear_lasso)
colnames(tisc1) <- c("outcome","vanilla scMR","ti-scMR")
colnames(tisc2) <- c("outcome","vanilla scMR","ti-scMR")
tisc1$DE <- de$FDR
tisc2$DE <- de$FDR
tisc1$GIFT <- gift$FDR
tisc2$GIFT <- gift$power
tisc1 <- tisc1 %>% pivot_longer(cols=2:5, names_to = "Method", values_to = "FDR")
tisc2 <- tisc2 %>% pivot_longer(cols=2:5, names_to = "Method", values_to = "Power")

tisc1$Method <- factor(tisc1$Method, levels = c("DE","GIFT","vanilla scMR","ti-scMR"))
tisc2$Method <- factor(tisc2$Method, levels = c("DE","GIFT","vanilla scMR","ti-scMR"))
library(ggthemes)
library(ggsci)
g1 <- ggplot(tisc1)+
  geom_bar(aes(outcome,FDR,fill = Method),stat="identity",position='dodge',colour="black")+
  theme_classic()+
  ggtitle("False Discovery Rate")+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_fill_aaas(alpha=0.6)
g2 <- ggplot(tisc2)+
  geom_bar(aes(outcome,Power,fill = Method),stat="identity",position='dodge',colour="black")+
  theme_classic()+
  ggtitle("Power")+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_fill_aaas(alpha=0.6)
library(patchwork)
g1+g2


