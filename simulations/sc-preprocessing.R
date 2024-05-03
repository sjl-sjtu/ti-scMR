library(Seurat)
library(data.table)
library(tidyverse)

p <- 10000 #genome size
n <- 500 #sample size
g <- 100 #gene num
t <- 20 #total time

# pre-processing
setwd("mixed_sample")
df <- fread("sim_sc.csv")
genelist <- paste0("gene",1:g)
mtx <- round(t(df[,genelist,with=F]))
rownames(mtx) <- colnames(df[,genelist,with=F])
colnames(mtx) <- df$cellid
cell_meta <- as.data.frame(df[,c("id","pseudotime","state","cellid")])
rownames(cell_meta) <- df$cellid
seurat_obj <- CreateSeuratObject(counts = mtx, meta.data = cell_meta)

seurat_obj <- NormalizeData(seurat_obj,normalization.method = "LogNormalize",
                            scale.factor = 10000)
df <- seurat_obj@assays$RNA@data %>% as.matrix() %>% t() %>% as_tibble()
df$id <- seurat_obj@meta.data$id
df$cellid <- colnames(seurat_obj)
df$pseudotime <- seurat_obj@meta.data$pseudotime
df$state <- seurat_obj@meta.data$state
df %>% write_csv("sim_sc_norm.csv")

all.genes <- rownames(seurat_obj)
seurat_obj <- ScaleData(seurat_obj, features = all.genes)
df <- seurat_obj@assays$RNA@scale.data %>% t() %>% as_tibble()
df$id <- seurat_obj@meta.data$id
df$cellid <- colnames(seurat_obj)
df$pseudotime <- seurat_obj@meta.data$pseudotime
df$state <- seurat_obj@meta.data$state
df %>% write_csv("sim_sc_scale.csv")

#######################################
#PACE
ge <- fread("sim_sc_scale.csv")
genelist <- paste0("gene",1:g)
exposureDat <- ge%>%
  group_by(id,pseudotime) %>%
  summarise(across(c(paste0("gene",1:g)),mean)) %>%
  ungroup() %>%
  drop_na()%>%
  mutate(sub=id)
exposureDat %>% write_csv("exposureDat.csv")

fullPACE <- function(geneName){
  s <- exposureDat %>% 
    # dplyr::select(id,any_of(outcomes),sub,pseudotime,any_of(geneName)) %>%
    dplyr::select(id,sub,pseudotime,any_of(geneName)) %>%
    dplyr::rename(expression=any_of(geneName)) %>%
    group_by(sub)  %>% 
    arrange(pseudotime) %>%
    group_map(~.x)
  N = length(s)
  lens <- map_dbl(s,nrow)
  if(sum(lens<=1)!=0){
    s <- s[-which(lens<=1)]
  }
  
  ly <- map(s,function(x) x$expression)
  lt <- map(s,function(x) x$pseudotime)
  res <- FPCA(ly,lt)
  Ti <- max(res$obsGrid)
  Tib <- min(res$obsGrid)
  Ti_est <- res$workGrid
  x_it <- t(matrix(replicate(length(s),res$mu),ncol=length(s)))+res$xiEst%*%t(res$phi)  #first term: mu(t) ->> n*t; second term: sum_1^K ksi_ik phi_k(t) ->> n*t
  
  getCurrEffe <- function(i){
    # ranges <- Ti_est<=Ti&Ti_est>=Tib
    # cums <- trapz(x=Ti_est[ranges],y=x_it[i,][ranges])
    # # cummeans <- cums/(Ti_est[ranges][length(Ti_est[ranges])]-Ti_est[ranges][1])
    # cums <- integrate(function(x) approx(Ti_est,x_it[i,],xout=x)$y,Tib,Ti)$value
    cums <- trapz(x=Ti_est,y=x_it[i,])
    return(cums)
  }
  xcum <- map_dbl(1:length(s),getCurrEffe)
  return(xcum)
}

library(doParallel)
library(foreach)
library(fdapace)
library(pracma)
registerDoParallel()
getDoParWorkers()
loaded_packages <- search()[grepl("package:", search())]
loaded_packages <- sub("package:", "", loaded_packages)

start_time <- Sys.time()
cum_mat <- foreach(geneName = genelist, .combine = 'cbind',.packages = loaded_packages) %dopar% {
  fullPACE(geneName)
}
end_time <- Sys.time()
run_time <- end_time - start_time
print(run_time)

colnames(cum_mat) <- genelist
cum_mat <- as_tibble(cum_mat)
cum_mat$sub <- exposureDat %>% group_by(sub) %>% group_keys() %>% pull(sub)
outcomes <- paste0("Y",1:7)
cum_mat <- cum_mat %>% left_join(exposureDat  %>% dplyr::select(id,any_of(outcomes),sub) %>% 
                                   dplyr::distinct(sub,.keep_all = T),by="sub")
cum_mat
cum_mat %>% write_csv("pace_cum.csv")


#################################
#preparation for eqtl

# df <- fread("sim_sc_norm.csv")
# df <- fread("true_sc.csv")
df <- fread("sim_sc_scale.csv")
# df <- fread("sim_sc_recons.csv")
# df <- fread("sim_sc_recons_log.csv")

df <- df %>%
  group_by(id) %>%
  summarise(across(c(paste0("gene",1:g)),function(x) mean(x[x!=0]))) %>%
  ungroup()

gene <- t(df[,c(paste0("gene",1:g)),with=F])
colnames(gene) <- df$id
gene <- as.data.table(gene)
gene$geneid <- paste0("gene",1:g)
gene <- gene[,c("geneid",df$id),with=FALSE]
gene %>% write_delim("GE.txt")

df <- fread("pace_cum.csv")
gene <- t(df[,c(paste0("gene",1:g)),with=F])
colnames(gene) <- df$id
gene <- as.data.table(gene)
gene$geneid <- paste0("gene",1:g)
gene <- gene[,c("geneid",df$id),with=FALSE]
gene %>% write_delim("GE_pace.txt")

dat <- fread("../geno.csv")
snp <- t(dat[,paste0("snp",1:p),with=FALSE])
colnames(snp) <- df$id
snp <- as.data.table(snp)
snp[,snpid:= paste0("snp",1:p)]
snp <- snp[,c("snpid",df$id),with=F]
snp %>% fwrite("SNP.txt",row.names=FALSE,sep=" ")


fit_prcomp <- prcomp(dat[,paste0("snp",1:p),with=FALSE],center = T,scale. = T)

m <- 10
covs <- fit_prcomp$x[,paste0("PC",1:m)]
covs <- covs|>t()|>as_tibble()
colnames(covs) <- df$id
covs$id <- paste0("PC",1:m)
covs <- covs%>%select(id,df$id)
covs %>% fwrite("Covariates.txt",row.names=FALSE,sep=" ")


#########################
# prepare for eQTL
df <- fread("sim_sc_scale.csv")

df <- df %>%
  mutate(grid=paste(id,pseudotime,sep="_")) %>%
  group_by(id,pseudotime,grid) %>%
  summarise(across(c(paste0("gene",1:g)),function(x) mean(x[x!=0]))) %>%
  ungroup()
gene <- t(df[,c(paste0("gene",1:g))])
colnames(gene) <- df$grid
gene <- as.data.table(gene)
gene$geneid <- paste0("gene",1:g)
gene <- gene[,c("geneid",df$grid),with=FALSE]
gene %>% write_delim("GE_sc.txt")

dat <- fread("../geno.csv")
dat <- as.data.table(df)[,.(id,grid)][dat,on=.(id)]
snp <- t(dat[,paste0("snp",1:p),with=FALSE])
colnames(snp) <- df$grid
snp <- as.data.table(snp)
snp[,snpid:= paste0("snp",1:p)]
snp <- snp[,c("snpid",df$grid),with=F]
snp %>% fwrite("SNP_sc.txt",row.names=FALSE,sep=" ")


covs <- as_tibble(t(df$pseudotime))
colnames(covs) <- df$grid
covs$id <- "pseudotime"
covs <- covs%>%select(id,df$grid)
covs %>% fwrite("Covariates_sc_t.txt",row.names=FALSE,sep=" ")

