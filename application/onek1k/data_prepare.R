library(Seurat)
library(data.table)
library(tidyverse)
library(fdapace)
library(pracma)
library(doParallel)
library(foreach)
cl <- makeCluster(20)
registerDoParallel(20)
getDoParWorkers()
loaded_packages <- search()[grepl("package:", search())]
loaded_packages <- sub("package:", "", loaded_packages)


markers <- fread("B_markers.csv")
genes <- markers[p_val_adj<0.01&abs(avg_log2FC)>2,gene]

seurat_obj <- readRDS("seurat_B_bc.rds")
seurat_obj <- SetIdent(seurat_obj,value = "predicted.celltype.l2")
seurat_obj <- ScaleData(seurat_obj,genes)

seurat_obj_bm <- subset(seurat_obj,idents=c("B naive","B intermediate",
                                          "B memory"))

df_bm_ti <- fread("B_memory_batchcorr.csv")

df_bm_scale <- seurat_obj_bm@assays$RNA@scale.data |> t()
df_bm_scale <- as.data.table(df_bm_scale,keep.rownames = T)

df_bm <- df_bm_scale[,c("rn",genes),with=F][df_bm_ti[,.(rn,donor_id,cell_type,sex,slingPseudotime_1)],on="rn"]

df_bm %>% fwrite("B_memory_scale_bc.csv")

df_bm$slingPseudotime_1 <- round(df_bm$slingPseudotime_1)
df_bm <- df_bm%>%dplyr::rename(pseudotime=slingPseudotime_1,id=donor_id)

exposureDat <- df_bm%>%
  group_by(id,pseudotime) %>%
  summarise(across(genes,mean)) %>%
  ungroup() %>%
  drop_na()%>%
  mutate(sub=id)
exposureDat %>% write_csv("exposureDat_bm_bc.csv")

# exposureDat <- read_csv("exposureDat_bm_bc.csv")
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

cum_mat |> dim()
colnames(cum_mat) <- genes
cum_mat <- as_tibble(cum_mat)
cum_mat$sub <- exposureDat %>% group_by(sub) %>% group_keys() %>% pull(sub)
cum_mat$outcome <- 0
cum_mat %>% write_csv("pace_cum_bm_bc.csv")

#######
seurat_obj_pl <- subset(seurat_obj,idents=c("B naive","B intermediate",
                                            "Plasmablast"))
df_pl_ti <- fread("plasa_batchcorr.csv")
df_pl_scale <- seurat_obj_pl@assays$RNA@scale.data |> t()
df_pl_scale <- as.data.table(df_pl_scale,keep.rownames = T)

df_pl <- df_pl_scale[,c("rn",genes),with=F][df_pl_ti[,.(rn,donor_id,cell_type,sex,slingPseudotime_1)],on="rn"]

df_pl %>% fwrite("Plasma_scale_bc.csv")

df_pl$slingPseudotime_1 <- round(df_pl$slingPseudotime_1)
df_pl <- df_pl%>%dplyr::rename(pseudotime=slingPseudotime_1,id=donor_id)

exposureDat <- df_pl%>%
  group_by(id,pseudotime) %>%
  summarise(across(genes,mean)) %>%
  ungroup() %>%
  drop_na()%>%
  mutate(sub=id)
exposureDat %>% write_csv("exposureDat_pl_bc.csv")

fullPACE <- function(geneName){
  s <- exposureDat %>% 
    dplyr::select(id,sub,pseudotime,any_of(geneName)) %>%
    dplyr::rename(expression=any_of(geneName)) %>%
    group_by(sub)  %>% 
    arrange(pseudotime) %>%
    group_map(~.x)
  N = length(s)
  
  ly <- map(s,function(x) x$expression)
  lt <- map(s,function(x) x$pseudotime)
  res <- FPCA(ly,lt)
  Ti <- max(res$obsGrid)
  Tib <- min(res$obsGrid)
  Ti_est <- res$workGrid
  x_it <- t(matrix(replicate(length(s),res$mu),ncol=length(s)))+res$xiEst%*%t(res$phi)  #first term: mu(t) ->> n*t; second term: sum_1^K ksi_ik phi_k(t) ->> n*t
  
  getCurrctfe <- function(i){
    # ranges <- Ti_est<=Ti&Ti_est>=Tib
    # cums <- trapz(x=Ti_est[ranges],y=x_it[i,][ranges])
    # # cummeans <- cums/(Ti_est[ranges][length(Ti_est[ranges])]-Ti_est[ranges][1])
    # cums <- integrate(function(x) approx(Ti_est,x_it[i,],xout=x)$y,Tib,Ti)$value
    cums <- trapz(x=Ti_est,y=x_it[i,])
    return(cums)
  }
  xcum <- map_dbl(1:length(s),getCurrctfe)
  return(xcum)
}


start_time <- Sys.time()
cum_mat <- foreach(geneName = genes, .combine = 'cbind',.packages = loaded_packages) %dopar% {
  fullPACE(geneName)
}
end_time <- Sys.time()
run_time <- end_time - start_time
print(run_time)

cum_mat |> dim()

colnames(cum_mat) <- genes
cum_mat <- as_tibble(cum_mat)
cum_mat$sub <- exposureDat %>% group_by(sub) %>% group_keys() %>% pull(sub)
cum_mat$outcome <- 1
cum_mat %>% write_csv("pace_cum_pl_bc.csv")


#################
library(data.table)
library(tidyverse)
setwd("/lustre/home/acct-clsyzs/clsyzs/SunJianle/singleCellMR/GSE196829/expression_matrix")

cum_mat <- fread("pace_cum_bm_bc.csv")
# cum_mat <- fread("pace_cum_pl_bc.csv")
cum_mat <- cum_mat[,id:=sub]
colnames(cum_mat)[(ncol(cum_mat)-5):ncol(cum_mat)]
genes <- grep("^ENSG",colnames(cum_mat),value=T)
cum2 <- cum_mat[,c("id",genes),with=F]
cum2[,1:5]
cum3 <- t(cum2[,..genes])
cum3 <- as.data.table(cum3)
colnames(cum3) <- cum2$id
cum3$geneid <- colnames(cum2)[-1]
cum3 <- cum3[,c("geneid",cum2$id),with=F]
cum2$id[1:5]

idinfo <- map_dfc(unique(cum2$id), function(x) str_split(x,"_")[[1]])
idinfo <- t(idinfo) |> as_tibble()
colnames(idinfo) <- c("FID","IID")
idinfo %>% write_delim("../genotype_chr/imputed/id_slt.txt",delim="\t")

id_order <- read_delim("/lustre/home/acct-clsyzs/clsyzs/SunJianle/singleCellMR/GSE196829/genotype_chr/imputed/slt_id_order.txt")
# id_order <- id_order %>% rowwise()%>% mutate(id=paste(FID,IID,sep="_"))
ids <- id_order$IID
cum3 <- cum3[,c("geneid",ids),with=F]

setwd("../eqtl/")
cum3 %>% write_delim("GE_bm.txt")
# cum3 %>% write_delim("GE_pl.txt")


#####
geneinfo <- fread("../expression_matrix/gene_list.txt",sep = "\t")
genes <- cum3$geneid
genepos <- geneinfo[ensembl_gene_id%in%genes][order(match(ensembl_gene_id, genes))]
genepos[,.(ensembl_gene_id,chromosome_name,start_position,end_position)] %>% write_delim("gene_pos.txt")

covs <- fread("../genotype_chr/imputed/plink2.eigenvec")

info <- fread("../GSE196829_series_matrix.txt",header = F,skip = 30, fill=T,sep="\t")
info <- info[-c(37,39),]
ns <- info$V1
info <- t(info[,-1])
info <- as.data.table(info)
colnames(info) <- ns
library(janitor)
info <- janitor::clean_names(info)
fwrite(info,"../id_info.txt",sep="\t")

info <- fread("../id_info.txt",sep="\t")
info[,id:=sub(".*: ", "", sample_characteristics_ch1)]
info <- info[id%in%ids][order(match(id, ids))]
info[,sex:=sub(".*: ", "", sample_characteristics_ch1_2)]
info[,.(id,sex)]
info[,sex2:=recode(sex,Male=0,Female=1,Unknown=0.5)]
info[,sex2]
info1 = info[,.(id,sex2)]
setnames(info1,"sex2","sex")
info1 <- covs[info1,on=c(`#IID`="id")]
info2 <- t(info1[,-1])
colnames(info2) <- info1$`#IID`
info2 <- as.data.table(rownames_to_column(as.data.frame(info2), var = "var"))
info2[,1:5]
info2 %>% fwrite("Covariates.txt",sep=" ")

snpinfo <- fread("../genotype_chr/imputed/imputed_snp_info.txt")
snplist <- fread("../genotype_chr/imputed/common_snp_slt.txt",sep=" ",header = F)
snplist <- snplist[-seq(1,6),]
colnames(snpinfo) <- c("id","chr","pos","ref","alt")
snpinfo[,snpid:=paste(id,ref,sep="_")]
snpinfo2 <- snpinfo[snpid%in%snplist$V1,.(snpid,chr,pos)][order(match(snpid, snplist$V1))]
snpinfo2 |> fwrite("snp_pos.txt",sep=" ")

snps <- fread("/lustre/home/acct-clsyzs/clsyzs/SunJianle/singleCellMR/GSE196829/genotype_chr/imputed/merged_slt_transpose.txt")
snps[,1:10]
snps_mat <- snps[-seq(1,6),]
colnames(snps_mat) <- snps[2,]|>unname()|>as_vector()
snps_mat[,1:10]
snps_mat |> fwrite("SNP.txt",sep=" ")


################
setwd("/lustre/home/acct-clsyzs/clsyzs/SunJianle/singleCellMR/GSE196829/eqtl/")
cum1 <- fread("GE_bm.txt")
cum2 <- fread("GE_pl.txt")
colnames(cum1)[-1] <- paste(colnames(cum1)[-1],"m",sep="_")
colnames(cum2)[-1] <- paste(colnames(cum2)[-1],"p",sep="_")
cum3 <- cum1[cum2,on="geneid"]
fwrite(cum3,"GE_all.txt",sep=" ")

snps <- fread("SNP.txt")
snps1 <- snps
colnames(snps)[-1] <- paste(colnames(snps)[-1],"m",sep="_")
colnames(snps1)[-1] <- paste(colnames(snps1)[-1],"p",sep="_")
snps <- snps[snps1,on="IID"]
snps <- snps[!duplicated(IID),]
fwrite(snps,"SNP_all.txt",sep=" ")

covs <- fread("Covariates.txt")
covs1 <- covs
covs <- covs%>%rbind(data.frame(matrix(c("outcome",rep(0,ncol(covs)-1)),nrow=1)),use.names=FALSE)
covs1 <- covs1%>%rbind(data.frame(matrix(c("outcome",rep(1,ncol(covs)-1)),nrow=1)),use.names=FALSE)
colnames(covs)[-1] <- paste(colnames(covs)[-1],"m",sep="_")
colnames(covs1)[-1] <- paste(colnames(covs1)[-1],"p",sep="_")
covs <- covs[covs1,on="var"]
fwrite(covs,"Covariates_all2.txt",sep=" ")
