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

snps <- fread("../geno.csv")
ge <- fread("sim_sc_scale.csv")
eqtl <- fread("eqtl_pace.csv")[FDR<0.05,]
exposureDat0 <- read_csv("exposureDat.csv")
cum_mat0 <- read_csv("pace_cum.csv")
timeLength <- 20
p <- 10000 #genome size
n <- 500 #sample size
g <- 100 #gene num
t <- 20 #total time
genelist <- paste0("gene",1:g)

source("scMR.R")

loaded_packages <- search()[grepl("package:", search())]
loaded_packages <- sub("package:", "", loaded_packages)

df <- fread("../true_sc.csv")
df <- df %>% 
  mutate(t=pseudotime/t) %>%
  arrange(id,t)%>%
  group_by(id) %>% 
  group_modify(~{.x %>% arrange(t) %>% 
      mutate(diff_t = t - lag(t,n = 1))}) %>%
  ungroup()
df[is.na(df$diff_t),"diff_t"] <- 0

k <- 4
gt <- 20
genenum <- rep(gt,k)
outcomes <- paste0("Y",1:k)

tr <- 50 # repeats
results <- list()
FDRs <- list()
powers <- list()

set.seed(12345)

for(r in 1:tr){
  print(r)

  # simulate phenotypes
  truegene <- replicate(k,paste0("gene",sample(1:g,gt)),simplify = F)
  
  out <- tibble(id = df %>% group_by(id) %>% group_keys() %>% pull(id))
  beta1 <- rnorm(genenum[1],0.1,0.05)
  out$Y1 <- df %>% group_by(id) %>% 
    group_modify(~{.x%>%arrange(t)}) %>% 
    group_map(~{
      sum((t(.x[,truegene[[1]]])*beta1)%*%as.matrix(.x[,"diff_t"]))
    }) %>% 
    unlist()
  beta2 <- beta1 #rnorm(genenum[2],0.05,0.01)
  out$Y2 <- df %>% group_by(id) %>% 
    group_modify(~{.x%>%arrange(t)}) %>% 
    group_map(~{
      sum(t(.x[,truegene[[2]]]*(cos(.x$t)%o%beta2))%*%as.matrix(.x[,"diff_t"]))
    }) %>% 
    unlist()
  beta3 <- beta1 #rnorm(genenum[2],0.05,0.01)
  out$Y3 <- df %>% group_by(id) %>% 
    group_modify(~{.x%>%arrange(t)}) %>% 
    group_modify(~{.x%>%mutate(indi = as.numeric(t<=0.2 | t>=0.8))}) %>% 
    group_map(~{
      sum(t(.x[,truegene[[3]]]*(.x$indi%o%beta3))%*%as.matrix(.x[,"diff_t"]))
    }) %>% 
    unlist()
  beta4 <- beta1 #rnorm(genenum[2],0.05,0.01)
  out$Y4 <- df %>% group_by(id) %>% 
    group_modify(~{.x%>%arrange(t)}) %>% 
    group_modify(~{.x%>%mutate(indi = as.numeric(t<=0.2 | t>=0.8))}) %>% 
    group_map(~{
      sum(t(.x[,truegene[[4]]]*((.x$indi*cos(.x$t))%o%beta4))%*%as.matrix(.x[,"diff_t"]))
    }) %>% 
    unlist()
  
  pl <- 100
  pleiosnps <- paste0("snp",sample(1:p,pl))
  gamma <- rnorm(pl,0.1,0.05)
  pe <- map_dbl(1:n,function(i) unlist(snps[i,..pleiosnps])%*%gamma)
  out[,2:(k+1)] <- apply(out[,2:(k+1)],2,function(x) x+pe)
  
  u <- rnorm(n,0,1)
  out[,2:(k+1)] <- apply(out[,2:(k+1)],2,function(x) x+u)
  colnames(out) <- c("id",paste0("Y",1:k))
  
  label <- tibble(gene=paste0("gene",1:g))
  labels <- map_dfc(1:k,function(i) as.numeric(paste0("gene",1:g) %in% truegene[[i]]))
  colnames(labels) <- paste0("Y",1:k)
  label <- as_tibble(cbind(label,labels))
  
  exposureDat <- exposureDat0 %>% left_join(out%>%dplyr::select(id,any_of(outcomes)),by="id")
  cum_mat <- cum_mat0 %>% left_join(out%>%dplyr::select(id,any_of(outcomes)),by="id")
  
  for(i in 1:k){
    outcome <- paste0("Y",i)

    res <- foreach(geneName = genelist,
                          .combine = 'rbind',.packages = loaded_packages) %dopar% {
                            MR_pace_linear(geneName, pca=F, m=0,  k=0)
                          }
    colnames(res) <- c("IVnum","F","beta","se","p")
    res <- as_tibble(res)
    res$gene <- genelist
    res$label <- label%>%pull(any_of(outcome))
    res$padj <- p.adjust(res$p,method = "BH")
    res$p[res$p == 0] <- 1e-16
    res$padj[res$padj == 0] <- 1e-16
    res %>% write_csv(paste0("quant_rep/pace_linear_",outcome,"_",r,".csv"))
    auc1 <- auc(res$label,-log10(res$p))
    auc2 <- auc(res$label,-log10(res$padj))
    auc3 <- auc(res$label,as.numeric(res$padj<0.05))
    FDR1 <- nrow(res%>%filter(padj<0.05)%>%filter(label==0))/nrow(res%>%filter(padj<0.05))
    power1 <- nrow(res%>%filter(label==1)%>%filter(padj<0.05))/gt
    
    res <- foreach(geneName = genelist,.packages = loaded_packages) %dopar% {
      MR_pace_lasso_linear(geneName, repeats = 5000)
    }
    res <- reduce(res,rbind)
    colnames(res) <- c("IVnum","F","beta","ci_l","ci_u","se","p")
    res <- as_tibble(res)
    res$gene <- genelist
    res$label <- label%>%pull(any_of(outcome))
    res$padj <- p.adjust(res$p,method = "BH")
    res$p[res$p == 0] <- 1e-16
    res$padj[res$padj == 0] <- 1e-16
    res %>% write_csv(paste0("quant_rep/pace_linear_lasso_",outcome,"_",r,".csv"))
    auc4 <- auc(res$label,-log10(res$p))
    auc5 <- auc(res$label,-log10(res$padj))
    auc6 <- auc(res$label,as.numeric(res$padj<0.05))
    FDR2 <- nrow(res%>%filter(padj<0.05)%>%filter(label==0))/nrow(res%>%filter(padj<0.05))
    power2 <- nrow(res%>%filter(label==1)%>%filter(padj<0.05))/gt
    
    res <- foreach(geneName = genelist,
                  .combine = 'rbind',.packages = loaded_packages) %dopar% {
                    mr_linear(geneName, pca=F, m=0, k=0)
                  }
    colnames(res) <- c("IVnum","F","beta","se","p")
    res <- as_tibble(res)
    res$gene <- genelist
    res$label <- label%>%pull(any_of(outcome))
    res$padj <- p.adjust(res$p,method = "BH")
    res$p[res$p == 0] <- 1e-16
    res$padj[res$padj == 0] <- 1e-16
    res %>% write_csv(paste0("quant_rep/avg_linear_",outcome,"_",r,".csv"))
    auc7 <- auc(res$label,-log10(res$p))
    auc8 <- auc(res$label,-log10(res$padj))
    auc9 <- auc(res$label,as.numeric(res$padj<0.05))
    FDR3 <- nrow(res%>%filter(padj<0.05)%>%filter(label==0))/nrow(res%>%filter(padj<0.05))
    power3 <- nrow(res%>%filter(label==1)%>%filter(padj<0.05))/gt

    res <- foreach(geneName = genelist,.packages = loaded_packages) %dopar% {
      MR_lasso_linear(geneName, repeats = 5000)
    }
    res <- reduce(res,rbind)
    colnames(res) <- c("IVnum","F","beta","ci_l","ci_u","se","p")
    res <- as_tibble(res)
    res$gene <- genelist
    res$label <- label%>%pull(any_of(outcome))
    res$padj <- p.adjust(res$p,method = "BH")
    res$p[res$p == 0] <- 1e-16
    res$padj[res$padj == 0] <- 1e-16
    res %>% write_csv(paste0("quant_rep/avg_linear_lasso_",outcome,"_",r,".csv"))
    auc10 <- auc(res$label,-log10(res$p))
    auc11 <- auc(res$label,-log10(res$padj))
    auc12 <- auc(res$label,as.numeric(res$padj<0.05))
    FDR4 <- nrow(res%>%filter(padj<0.05)%>%filter(label==0))/nrow(res%>%filter(padj<0.05))
    power4 <- nrow(res%>%filter(label==1)%>%filter(padj<0.05))/gt
    
    results <- append(results,list(c(r,outcome,auc1,auc2,auc3,
                                     auc4,auc5,auc6,auc7,auc8,auc9,
                                     auc10,auc11,auc12)))
    FDRs <- append(FDRs,list(c(r,outcome,FDR1,FDR2,FDR3,FDR4)))
    powers <- append(powers,list(c(r,outcome,power1,power2,power3,power4)))
  }
}

results <- do.call(rbind,results)
results <- results %>% as_tibble() 
colnames(results) <- c("repeats","outcome","pace_linear","pace_linear2","pace_linear3",
                       "pace_linear_lasso","pace_linear_lasso2","pace_linear_lasso3",
                       "avg_linear","avg_linear2","avg_linear3","avg_linear_lasso",
                       "avg_linear_lasso2","avg_linear_lasso3")
results %>% write_csv("quant_rep_mixed_auc.csv",append=T)

FDRs <- do.call(rbind,FDRs)
FDRs <- FDRs %>% as_tibble() 
colnames(FDRs) <- c("repeats","outcome","pace_linear",
                    "pace_linear_lasso","avg_linear","avg_linear_lasso")
FDRs %>% write_csv("quant_rep_mixed_FDR.csv",append=T)

powers <- do.call(rbind,powers)
powers <- powers %>% as_tibble() 
colnames(powers) <- c("repeats","outcome","pace_linear","pace_linear_lasso",
                      "avg_linear","avg_linear_lasso")
powers %>% write_csv("quant_rep_mixed_power.csv",append=T)
