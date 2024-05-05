library(data.table)
library(tidyverse)
setwd("/lustre/home/acct-clsyzs/clsyzs/SunJianle/singleCellMR/GSE196829/expression_matrix")


idinfo <- map_dfc(unique(cum2$id), function(x) str_split(x,"_")[[1]])
idinfo <- t(idinfo) |> as_tibble()
colnames(idinfo) <- c("FID","IID")
idinfo %>% write_delim("../genotype_chr/imputed/id_slt.txt",delim="\t")

# cum_mat <- map_dfr(c("pace_cum_bm.csv","pace_cum_pl.csv"),fread)

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



id_order <- read_delim("/lustre/home/acct-clsyzs/clsyzs/SunJianle/singleCellMR/GSE196829/genotype_chr/imputed/slt_id_order.txt")
#id_order <- id_order %>% rowwise()%>% mutate(id=paste(FID,IID,sep="_"))
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
