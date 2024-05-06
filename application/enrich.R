
library(clusterProfiler)
library(org.Hs.eg.db)

setwd("/lustre/home/acct-clsyzs/clsyzs/SunJianle/singleCellMR/GSE196829/expression_matrix/results_b/pace_all_new3")

df <- read_csv("re_pace_lasso_outcome.csv")%>%drop_na()
library(qvalue)
df <- df %>%drop_na()%>%  mutate(padj=p.adjust(p,method="BH"),q=qvalue(p,lambda=0)$qvalues) 
df1 <- df %>% filter(padj<0.01)
df1
gene1 <- df1$gene
bitr(gene1, fromType="ENSEMBL", toType=c("SYMBOL", "ENTREZID"), 
     OrgDb="org.Hs.eg.db")

info <- fread("../../gene_list.txt")
info1 <- info[ensembl_gene_id%in%gene1,]

markers <- fread("../../B_markers.csv")
# ms <- bitr(markers$gene, fromType="ENSEMBL", toType=c("SYMBOL", "ENTREZID"), 
#            OrgDb="org.Hs.eg.db")

genes <- info%>%filter(ensembl_gene_id%in%gene1)%>%pull(external_gene_name)

test = bitr(gene1, fromType="ENSEMBL", toType=c("SYMBOL", "ENTREZID"), 
            OrgDb="org.Hs.eg.db")

enrich.go <- enrichGO(gene = gene1, #gene3,#test$SYMBOL,# test$SYMBOL,# genes,  #待富集的基因列表
                      OrgDb = org.Hs.eg.db,  #指定物种的基因数据库，示例物种是绵羊（sheep）
                      keyType = 'ENSEMBL',#'SYMBOL',#'ENSEMBL',  #指定给定的基因名称类型，例如这里以 entrze id 为例
                      ont = 'ALL',  #GO Ontology，可选 BP、MF、CC，也可以指定 ALL 同时计算 3 者
                      pAdjustMethod = 'BH',  #指定 p 值校正方法
                      pvalueCutoff = 0.01,  #指定 p 值阈值（可指定 1 以输出全部）
                      qvalueCutoff = 0.05,  #指定 q 值阈值（可指定 1 以输出全部）
                      readable = FALSE)
enrich.go1 <- enrichGO(gene = test$SYMBOL, #gene3,#test$SYMBOL,# test$SYMBOL,# genes,  #待富集的基因列表
                      OrgDb = org.Hs.eg.db,  #指定物种的基因数据库，示例物种是绵羊（sheep）
                      keyType = 'SYMBOL',#'ENSEMBL',  #指定给定的基因名称类型，例如这里以 entrze id 为例
                      ont = 'ALL',  #GO Ontology，可选 BP、MF、CC，也可以指定 ALL 同时计算 3 者
                      pAdjustMethod = 'BH',  #指定 p 值校正方法
                      pvalueCutoff = 0.01,  #指定 p 值阈值（可指定 1 以输出全部）
                      qvalueCutoff = 0.05,  #指定 q 值阈值（可指定 1 以输出全部）
                      readable = FALSE)

enrich.go
library(enrichplot)
barplot(enrich.go,title = "GO enrichment",font.size=10)

p0 <- barplot(enrich.go,title = "GO enrichment",font.size=10) 
p1 <- dotplot(enrich.go,font.size=10) 
p2 <- cnetplot(enrich.go1,cex_label_category=0.6,cex_label_gene=0.6,cex_category=0.6,cex_gene=0.6) 
emapplot(pairwise_termsim(enrich.go,method="JC"),cex_label_category=0.6,cex_category=0.6,cex_line=0.1) 
heatplot(enrich.go1) 
library(patchwork)
((p0+p1+plot_layout(width=c(1,1)))/(p2))+ plot_layout(heights = c(2, 3))+plot_annotation(tag_levels = "a")
((p0|p1)/p2)+plot_annotation(tag_levels = "a")
((p0/p1)|p2)+ plot_layout(widths = c(1, 2))#+plot_annotation(tag_levels = "a")


enrich.kegg <- enrichKEGG(gene = test$ENTREZID,
                          organism = 'hsa',
                          pvalueCutoff = 0.05,
                          pAdjustMethod = "BH")
dotplot(enrich.kegg,title = "KEGG enrichment")

res <- df1%>%left_join(markers)%>%left_join(test,by=c("gene"="ENSEMBL"))%>%left_join(info1,by=c("gene"="ensembl_gene_id","SYMBOL"="external_gene_name"))
res<-res%>%arrange(q)
res%>%write_csv("results.csv")

###############
setwd("/lustre/home/acct-clsyzs/clsyzs/SunJianle/singleCellMR/Roche/MR/results_oligo/pace_all_new4/")
df <- read_csv("re_pace_lasso_outcome.csv")%>%drop_na()
df <- df %>%drop_na()%>%  mutate(padj=p.adjust(p,method="BH"),q=qvalue(p,lambda=0)$qvalues) 
df1 <- df %>% filter(padj<0.1)
df1
gene1 <- df1$gene

test = bitr(gene1, fromType="ENSEMBL", toType=c("SYMBOL", "ENTREZID"), 
            OrgDb="org.Hs.eg.db")
test

enrich.go <- enrichGO(gene = gene1,#test$SYMBOL,# test$SYMBOL,# genes,  #待富集的基因列表
                      OrgDb = org.Hs.eg.db,  #指定物种的基因数据库，示例物种是绵羊（sheep）
                      keyType = 'ENSEMBL',#'SYMBOL',#'ENSEMBL',  #指定给定的基因名称类型，例如这里以 entrze id 为例
                      ont = 'ALL',  #GO Ontology，可选 BP、MF、CC，也可以指定 ALL 同时计算 3 者
                      pAdjustMethod = 'BH',  #指定 p 值校正方法
                      pvalueCutoff = 0.05,  #指定 p 值阈值（可指定 1 以输出全部）
                      qvalueCutoff = 0.2,  #指定 q 值阈值（可指定 1 以输出全部）
                      readable = FALSE)
enrich.go1 <- enrichGO(gene = test$SYMBOL,# test$SYMBOL,# genes,  #待富集的基因列表
                      OrgDb = org.Hs.eg.db,  #指定物种的基因数据库，示例物种是绵羊（sheep）
                      keyType = 'SYMBOL',#'ENSEMBL',  #指定给定的基因名称类型，例如这里以 entrze id 为例
                      ont = 'ALL',  #GO Ontology，可选 BP、MF、CC，也可以指定 ALL 同时计算 3 者
                      pAdjustMethod = 'BH',  #指定 p 值校正方法
                      pvalueCutoff = 0.05,  #指定 p 值阈值（可指定 1 以输出全部）
                      qvalueCutoff = 0.2,  #指定 q 值阈值（可指定 1 以输出全部）
                      readable = FALSE)

enrich.go
library(enrichplot)
barplot(enrich.go,title = "GO enrichment",font.size=10)

p0 <- barplot(enrich.go,title = "GO enrichment",font.size=12) 
p1 <- dotplot(enrich.go,font.size=12)
p2 <- cnetplot(enrich.go1,cex_label_category=1,cex_label_gene=1,cex_category=1,cex_gene=1)
emapplot(pairwise_termsim(enrich.go,method="JC"),cex_label_category=0.6,cex_category=0.6,cex_line=0.1)
heatplot(enrich.go1)
library(patchwork)
((p0+p1+plot_layout(width=c(1,1)))/(p2))+ plot_layout(heights = c(2, 3))+plot_annotation(tag_levels = "a")
((p0|p1)/p2)+plot_annotation(tag_levels = "a")
# ridgeplot(enrich.go)
((p0/p1)|p2)+ plot_layout(widths = c(1, 2))#+plot_annotation(tag_levels = "a")

markers <- fread("../../../Oligo_Ctr_MS_markers3.csv")
ms <- bitr(markers$gene, fromType="ENSEMBL", toType=c("SYMBOL", "ENTREZID"), 
           OrgDb="org.Hs.eg.db")
info <- fread("../../../eQTL/gene_list.txt")
info1 <- info[ensembl_gene_id%in%gene1,]
res <- df1%>%left_join(markers)%>%left_join(ms,by=c("gene"="ENSEMBL"))%>%left_join(info1,by=c("gene"="ensembl_gene_id","SYMBOL"="external_gene_name"))
res <- res%>%arrange(q)
res
res%>%write_csv("results2.csv")
