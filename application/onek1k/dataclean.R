library(Seurat)
library(tidyverse)
library(data.table)

seurat_obj1 <- readRDS("local.rds")
df1 <- seurat_obj@meta.data |> as_tibble()
seurat_obj <- readRDS("onek1k.rds")
df <- seurat_obj@meta.data |> rownames_to_column(var = "barcode") |> as_tibble()
df %>% fwrite("cell_info.csv")

eqtl <- fread("esnp_table.tsv")
eqtl$FDR |> max()
snplist <- tibble(SNP=eqtl$SNPID |> unique()) %>% write_delim("genotype_chr/imputed/SNP.txt")


dfgene <- fread("genotype_modified.raw")

info <- fread("GSE196829_series_matrix.txt",header = F,skip = 30, fill=T,sep="\t")
info <- info[-c(37,39),]
ns <- info$V1
info <- t(info[,-1])
info <- as.data.table(info)
colnames(info) <- ns

indi_id <- sub(".*: ", "", info$`!Sample_characteristics_ch1`)
illu_id <- sub(".*_(\\d+_R\\d+C\\d+)_.*", "\\1",info$`!Sample_supplementary_file`)
info_match <- data.table(individual=indi_id,illumina=illu_id)

dfinfo <- merge(dfgene[,1:6],info_match,all.x = T,by.x = "IID",by.y = "illumina")[order(FID)]
dfgene$individual <- dfinfo$individual
dfgene %>% fwrite("genotype_modified.txt")

info <- fread("/lustre/home/acct-clsyzs/clsyzs/SunJianle/singleCellMR/GSE196829/GSE196829_series_matrix.txt",skip=30,header=F)
ins <- t(info)
ins <- as.data.table(ins)
colnames(ins) <- info$V1
ins <- ins[-1,]

ids <- sub(".*/.*_(.*_.*)_Grn\\.idat\\.gz", "\\1", ins$`!Sample_supplementary_file`)
id <- sub(".*: (.*)$","\\1",ins$`!Sample_characteristics_ch1`)
idss <- data.table(id=id,illumina=ids)
tibble(id=id,illumina=ids) %>% write_csv("id_info.csv")
df1 <- idss[dfgene,on=c("illumina"="IID")]
df1[,.(id,individual)]
df1 %>% fwrite("genotype_modified.txt")

snps <- fread("../genotype_modified.txt")
colnames(snps) <- sub("^(rs[^_]+).*", "\\1", colnames(snps))
Mode <- function(x, na.rm = FALSE) {
  if (na.rm) {
    x <- x[!is.na(x)]
  }
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

snps_filled <- snps[, lapply(.SD, function(x) ifelse(is.na(x), Mode(x, na.rm = TRUE), x))]
snps_filled %>% fwrite("../genotype_modified_imputated.txt")

eqtl <- fread("../onek1k_eqtl.tsv")
eqtl1 <- eqtl[FDR<0.05,]
eqtl1 %>% fwrite("../sig_eqtl.tsv")

snps <- fread("../genotype_modified_imputated.txt")
ss <- grep("^rs",colnames(snps),value = T,invert = T)
ss <- gsub("_(.*)", "", ss)
eqtl <- eqtl[,SNP:=gsub("_(.*)", "", SNPID)]
hh <- intersect(ss,eqtl$SNP)
df1 <- eqtl[SNP%in%hh,.(SNP,SNPID,RSID)]|>unique()
col_names <- colnames(snps)
col_names_to_replace <- df1$SNPID
col_names_new <- df1$RSID
for (i in seq_along(col_names_to_replace)) {
  col_names <- ifelse(col_names == col_names_to_replace[i], col_names_new[i], col_names)
}
colnames(snps) <- col_names

mm <- setdiff(eqtl$RSID,colnames(snps))
df2 <- eqtl[GENOTYPED=="Genotyped"&RSID%in%mm,.(SNP,SNPID,RSID)]|>unique()
df2
snps %>% select(starts_with("1:229073815"))
ss <- grep("^rs",colnames(snps),value = T,invert = T)
ss <- grep("^[1-9]",ss,value=T)
sort(ss)

library(SNPlocs.Hsapiens.dbSNP144.GRCh37)
all_snps <- SNPlocs.Hsapiens.dbSNP144.GRCh37
chr_snps <- snpsBySeqname(all_snps, c(as.character(seq(1,22)),"X"))
chr_snps <- as.data.table(chr_snps)

################3
setwd("/lustre/home/acct-clsyzs/clsyzs/SunJianle/singleCellMR/GSE196829/expression_matrix")
library(Seurat)
library(data.table)
library(tidyverse)
library(biomaRt)
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
seurat_obj <- readRDS("../onek1k.rds")
gene_ids <- rownames(seurat_obj)
gene_info <- getBM(attributes = c("ensembl_gene_id","external_gene_name",
                                  "chromosome_name","start_position", 
                                  "end_position","strand"), 
                   filters = "ensembl_gene_id", values = gene_ids, mart = ensembl)
gene_info |> as.data.table() |> fwrite("gene_list.txt",sep = "\t")
