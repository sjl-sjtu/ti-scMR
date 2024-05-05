setwd("~/SunJianle/singleCellMR/Roche/eQTL")
library(data.table)
library(tidyverse)

###PACE
library(Seurat)
library(data.table)
library(tidyverse)
library(fdapace)
library(pracma)
library(doParallel)
library(foreach)
cl <- makeCluster(10)
registerDoParallel(10)
getDoParWorkers()
loaded_packages <- search()[grepl("package:", search())]
loaded_packages <- sub("package:", "", loaded_packages)

markers <- fread("../Oligo_Ctr_MS_markers3.csv")
genes <- markers[p_val_adj<0.05&abs(avg_log2FC)>0.8,gene]

seurat_ms_Oligo <- readRDS("../seurat_ms_Oligo_batch.rds")
seurat_ctr_Oligo <- readRDS("../seurat_ctr_Oligo_batch.rds")
seurat_ms_Oligo <- ScaleData(seurat_ms_Oligo,genes)
seurat_ctr_Oligo <- ScaleData(seurat_ctr_Oligo,genes)

df_ms_ti <- fread("../ms_Oligo_batch.csv")

df_ms_scale <- seurat_ms_Oligo@assays$RNA$scale.data |> t()
df_ms_scale <- as.data.table(df_ms_scale,keep.rownames = T)

df_ms <- df_ms_scale[,c("rn",genes),with=F][df_ms_ti[,.(rn,individual_id_anon,sex,age_at_death,ident,slingPseudotime_1)],on="rn"]

df_ms %>% fwrite("ms_scale_bc.csv")

df_ms$slingPseudotime_1 <- round(df_ms$slingPseudotime_1)
df_ms <- df_ms%>%dplyr::rename(pseudotime=slingPseudotime_1,id=individual_id_anon)

exposureDat <- df_ms%>%
  group_by(id,pseudotime) %>%
  summarise(across(genes,mean)) %>%
  ungroup() %>%
  drop_na()%>%
  mutate(sub=id)
exposureDat %>% write_csv("exposureDat_ms_bc.csv")

fullPACE <- function(geneName){
  s <- exposureDat %>% 
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
    cums <- trapz(x=Ti_est,y=x_it[i,])#/(Ti-Tib)
    return(cums)
  }
  xcum <- map_dbl(1:length(s),getCurrEffe)
  return(xcum)
}


start_time <- Sys.time()
cum_mat <- foreach(geneName = genes,.combine = 'cbind', .packages = loaded_packages) %dopar% {
  fullPACE(geneName)
}
end_time <- Sys.time()
run_time <- end_time - start_time
print(run_time)

#

cum_mat |> dim()

colnames(cum_mat) <- genes
cum_mat <- as_tibble(cum_mat)
cum_mat$sub <- exposureDat %>% group_by(sub) %>% group_keys() %>% pull(sub)
cum_mat$outcome <- 1
cum_mat %>% write_csv("pace_cum_ms_bc.csv")


df_ctr_ti <- fread("../ctr_Oligo_batch.csv")

df_ctr_scale <- seurat_ctr_Oligo@assays$RNA$scale.data |> t()
df_ctr_scale <- as.data.table(df_ctr_scale,keep.rownames = T)

df_ctr <- df_ctr_scale[,c("rn",genes),with=F][df_ctr_ti[,.(rn,individual_id_anon,sex,age_at_death,ident,slingPseudotime_1)],on="rn"]

df_ctr %>% fwrite("ctr_scale_bc.csv")

df_ctr$slingPseudotime_1 <- round(df_ctr$slingPseudotime_1)
df_ctr <- df_ctr%>%dplyr::rename(pseudotime=slingPseudotime_1,id=individual_id_anon)

exposureDat <- df_ctr%>%
  group_by(id,pseudotime) %>%
  summarise(across(genes,mean)) %>%
  ungroup() %>%
  drop_na()%>%
  mutate(sub=id)
exposureDat %>% write_csv("exposureDat_ctr_bc.csv")

fullPACE <- function(geneName){
  s <- exposureDat %>% 
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
    # cuctr <- trapz(x=Ti_est[ranges],y=x_it[i,][ranges])
    # # cummeans <- cuctr/(Ti_est[ranges][length(Ti_est[ranges])]-Ti_est[ranges][1])
    # cuctr <- integrate(function(x) approx(Ti_est,x_it[i,],xout=x)$y,Tib,Ti)$value
    cuctr <- trapz(x=Ti_est,y=x_it[i,])#/(Ti-Tib)
    return(cuctr)
  }
  xcum <- map_dbl(1:length(s),getCurrEffe)
  return(xcum)
}


start_time <- Sys.time()
cum_mat <- foreach(geneName = genes,.combine = 'cbind', .packages = loaded_packages) %dopar% {
  fullPACE(geneName)
}
end_time <- Sys.time()
run_time <- end_time - start_time
print(run_time)

#

cum_mat |> dim()

colnames(cum_mat) <- genes
cum_mat <- as_tibble(cum_mat)
cum_mat$sub <- exposureDat %>% group_by(sub) %>% group_keys() %>% pull(sub)
cum_mat$outcome <- 0
cum_mat %>% write_csv("pace_cum_ctr_bc.csv")


####
cum_mat1 <- fread("pace_cum_ms_bc.csv")
cum_mat0 <- fread("pace_cum_ctr_bc.csv")
cum_mat<-rbind(cum_mat0,cum_mat1)
genes <- grep("^ENSG",colnames(cum_mat),value=T)
cum_mat <- cum_mat[,id:=sub]
cum2 <- cum_mat[,c("id",genes),with=F]
cum3 <- t(cum2[,..genes])
cum3 <- as.data.table(cum3)
colnames(cum3) <- cum2$id
cum3$geneid <- colnames(cum2)[-1]
cum3 <- cum3[,c("geneid",cum2$id),with=F]
cum2$id[1:5]

idinfo <- fread("../EGAF00006717330/imputed_hg38/id_list.txt",header = F)
idinfo
ids <- idinfo$V1
ids2 <- intersect(colnames(cum3),ids)
fwrite(data.table(ID=ids2),"id_slt.txt")

ids2 <- fread("id_slt.txt")[,ID]
cum3 <- cum3[,c("geneid",ids2),with=F]
cum3 %>% write_delim("GE.txt")

###
library(biomaRt)
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
gene_ids <- rownames(seurat_ms_Oligo)
gene_info <- getBM(attributes = c("ensembl_gene_id","external_gene_name",
                                  "chromosome_name","start_position", 
                                  "end_position","strand"), 
                   filters = "ensembl_gene_id", values = gene_ids, mart = ensembl)
gene_info <- as.data.table(gene_info)
gene_info |> as.data.table() |> fwrite("gene_list.txt",sep = "\t")

gene_info <- fread("gene_list.txt",sep = "\t")
genes <- cum3$geneid
genepos <- gene_info[ensembl_gene_id%in%genes][order(match(ensembl_gene_id, genes))]
genepos[,.(ensembl_gene_id,chromosome_name,start_position,end_position)] %>% write_delim("gene_pos.txt")

cum3 <- cum3[geneid%in%genepos$ensembl_gene_id][order(match(geneid,genepos$ensembl_gene_id))]
cum3 %>% write_delim("GE.txt")

####
covs <- fread("../EGAF00006717330/imputed_hg38/plink2.eigenvec")
info <- fread("../ms_lesson_raw/barcodes.txt")
info <- info[,.(individual_id_anon,sex,age_scale)]|>unique()
info <- info[sex == 'M', sex := 1][sex == 'F', sex := 0]
info$sex <- as.numeric(info$sex)
info <- covs[info[individual_id_anon%in%ids2,.(individual_id_anon,sex,age_scale)],on=c("#IID"="individual_id_anon")]
info
covname <- colnames(info)
idss <- info$`#IID`
info <- t(info[,-1])
info <- as.data.table(info)
colnames(info) <- idss
info$id <- covname[-1]
info<- info[,c("id",ids2),with=F]
info %>% fwrite("Covariates.txt",sep=" ")


###
snps <- fread("../EGAF00006717330/imputed_hg38/merged_slt_transpose.txt")
snps[,1:10]
snps_mat <- snps[-seq(1,6),]
colnames(snps_mat) <- c("id",snps[2,-1]|>unname()|>as_vector())
snps_mat[,1:10]
snps_mat <-snps_mat[,c("id",ids2),with=F]
snps_mat |> fwrite("SNP.txt",sep=" ")

snpsinfo <- str_split(snps_mat$id,"[:_]",simplify = T)
snpsinfo <- as_tibble(snpsinfo)
colnames(snpsinfo) <- c("chr","pos","alt","ref","as")
snpsinfo$chr <- gsub("[^0-9]", "", snpsinfo$chr)
snpsinfo <- snpsinfo[,1:2]
snpsinfo <- snpsinfo%>%mutate_all(as.numeric)
snpsinfo$snpid <- snps_mat$id
snpsinfo <- snpsinfo%>%dplyr::select(snpid,chr,pos)
snpsinfo |> fwrite("snp_pos.txt",sep=" ")
