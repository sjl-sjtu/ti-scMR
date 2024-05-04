library(MatrixEQTL)

base.dir = getwd()

### avg_eqtl
useModel = modelLINEAR 
SNP_file_name = paste(base.dir, "/SNP.txt", sep="") 
expression_file_name = paste(base.dir, "/GE.txt", sep="") 
covariates_file_name = character() #
# covariates_file_name = paste(base.dir, "/Covariates.txt", sep="") 
output_file_name = tempfile() 
pvOutputThreshold = 1 
errorCovariance = numeric() 

snps = SlicedData$new() 
snps$fileDelimiter = " "       
snps$fileOmitCharacters = "NA" 
snps$fileSkipRows = 1        
snps$fileSkipColumns = 1     
snps$fileSliceSize = 2000      
snps$LoadFile( SNP_file_name ) 

gene = SlicedData$new()
gene$fileDelimiter = " " 
gene$fileOmitCharacters = "NA"
gene$fileSkipRows = 1
gene$fileSkipColumns = 1
gene$fileSliceSize = 2000
gene$LoadFile(expression_file_name)

cvrt = SlicedData$new();
cvrt$fileDelimiter = " ";      # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
if(length(covariates_file_name)>0) {
  cvrt$LoadFile(covariates_file_name);
}


maf.list = vector('list', length(snps))
for(sl in 1:length(snps)) {
  slice = snps[[sl]];
  maf.list[[sl]] = rowMeans(slice,na.rm=TRUE)/2;
  maf.list[[sl]] = pmin(maf.list[[sl]],1-maf.list[[sl]]);
}
maf = unlist(maf.list)

snps$RowReorder(maf>=0.01);

me = Matrix_eQTL_engine(  #进行eQTL分析的主要函数
  snps = snps,  #指定SNP 文件
  gene = gene,  #指定基因表达量文件
  cvrt = cvrt,   #指定协变量文件
  output_file_name = output_file_name,  #指定输出文件
  pvOutputThreshold = pvOutputThreshold,  #指定显著性P值
  useModel = useModel,  #指定使用的计算模型
  errorCovariance = errorCovariance, #指定误差项的协方差矩阵
  verbose = TRUE,
  min.pv.by.genesnp = FALSE,  pvalue.hist = "qqplot",

  noFDRsaveMemory = FALSE)

eqtl <- me$all$eqtls
eqtl$gene |> unique() %>% length()
eqtl %>% write_csv("eqtl_avg_nocov.csv")

dfiv <- read_csv("../dfiv.csv")
ivs <- dfiv %>%
  pivot_longer(-snp, names_to = "gene", values_to = "Value") %>%
  filter(Value == 1) %>%
  select(gene, snp) %>%
  mutate(pair=paste(gene,snp,sep="_"))
eqtl <- as.data.table(eqtl)[FDR<0.05,][,pair:=paste(gene,snps,sep="_")]

pre1 <- (eqtl[pair%in%ivs$pair,]|>nrow())/(eqtl|>nrow()) #TDR
recal1 <- (eqtl[pair%in%ivs$pair,]|>nrow())/(ivs|>nrow()) #Recall

##
### avg_eqtl + 10PCs
useModel = modelLINEAR 
SNP_file_name = paste(base.dir, "/SNP.txt", sep="") 
expression_file_name = paste(base.dir, "/GE.txt", sep="") 
# covariates_file_name = character() #
covariates_file_name = paste(base.dir, "/Covariates.txt", sep="")
output_file_name = tempfile() 
pvOutputThreshold = 1 
errorCovariance = numeric() 

snps = SlicedData$new() 
snps$fileDelimiter = " "       
snps$fileOmitCharacters = "NA" 
snps$fileSkipRows = 1        
snps$fileSkipColumns = 1     
snps$fileSliceSize = 2000      
snps$LoadFile( SNP_file_name ) 

gene = SlicedData$new()
gene$fileDelimiter = " " 
gene$fileOmitCharacters = "NA"
gene$fileSkipRows = 1
gene$fileSkipColumns = 1
gene$fileSliceSize = 2000
gene$LoadFile(expression_file_name)

cvrt = SlicedData$new();
cvrt$fileDelimiter = " ";      # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
if(length(covariates_file_name)>0) {
  cvrt$LoadFile(covariates_file_name);
}


maf.list = vector('list', length(snps))
for(sl in 1:length(snps)) {
  slice = snps[[sl]];
  maf.list[[sl]] = rowMeans(slice,na.rm=TRUE)/2;
  maf.list[[sl]] = pmin(maf.list[[sl]],1-maf.list[[sl]]);
}
maf = unlist(maf.list)

snps$RowReorder(maf>=0.01);

me = Matrix_eQTL_engine(  #进行eQTL分析的主要函数
  snps = snps,  #指定SNP 文件
  gene = gene,  #指定基因表达量文件
  cvrt = cvrt,   #指定协变量文件
  output_file_name = output_file_name,  #指定输出文件
  pvOutputThreshold = pvOutputThreshold,  #指定显著性P值
  useModel = useModel,  #指定使用的计算模型
  errorCovariance = errorCovariance, #指定误差项的协方差矩阵
  verbose = TRUE,
  min.pv.by.genesnp = FALSE,  pvalue.hist = "qqplot",
  
  noFDRsaveMemory = FALSE)

eqtl <- me$all$eqtls
eqtl$gene |> unique() %>% length()
eqtl %>% write_csv("eqtl_avg.csv")

dfiv <- read_csv("../dfiv.csv")
ivs <- dfiv %>%
  pivot_longer(-snp, names_to = "gene", values_to = "Value") %>%
  filter(Value == 1) %>%
  select(gene, snp) %>%
  mutate(pair=paste(gene,snp,sep="_"))
eqtl <- as.data.table(eqtl)[FDR<0.05,][,pair:=paste(gene,snps,sep="_")]

pre2 <- (eqtl[pair%in%ivs$pair,]|>nrow())/(eqtl|>nrow()) #TDR
recal2 <- (eqtl[pair%in%ivs$pair,]|>nrow())/(ivs|>nrow()) #Recall

### pace_eqtl
useModel = modelLINEAR 
SNP_file_name = paste(base.dir, "/SNP.txt", sep="") 
expression_file_name = paste(base.dir, "/GE_pace.txt", sep="") 
covariates_file_name = character() #
# covariates_file_name = paste(base.dir, "/Covariates.txt", sep="") 
output_file_name = tempfile() 
pvOutputThreshold = 1 
errorCovariance = numeric() 

snps = SlicedData$new() 
snps$fileDelimiter = " "       
snps$fileOmitCharacters = "NA" 
snps$fileSkipRows = 1        
snps$fileSkipColumns = 1     
snps$fileSliceSize = 2000      
snps$LoadFile( SNP_file_name ) 

gene = SlicedData$new()
gene$fileDelimiter = " " 
gene$fileOmitCharacters = "NA"
gene$fileSkipRows = 1
gene$fileSkipColumns = 1
gene$fileSliceSize = 2000
gene$LoadFile(expression_file_name)

cvrt = SlicedData$new();
cvrt$fileDelimiter = " ";      # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
if(length(covariates_file_name)>0) {
  cvrt$LoadFile(covariates_file_name);
}

maf.list = vector('list', length(snps))
for(sl in 1:length(snps)) {
  slice = snps[[sl]];
  maf.list[[sl]] = rowMeans(slice,na.rm=TRUE)/2;
  maf.list[[sl]] = pmin(maf.list[[sl]],1-maf.list[[sl]]);
}
maf = unlist(maf.list)

snps$RowReorder(maf>=0.01);

me = Matrix_eQTL_engine(  #进行eQTL分析的主要函数
  snps = snps,  #指定SNP 文件
  gene = gene,  #指定基因表达量文件
  cvrt = cvrt,   #指定协变量文件
  output_file_name = output_file_name,  #指定输出文件
  pvOutputThreshold = pvOutputThreshold,  #指定显著性P值
  useModel = useModel,  #指定使用的计算模型
  errorCovariance = errorCovariance, #指定误差项的协方差矩阵
  verbose = TRUE,
  min.pv.by.genesnp = FALSE,  pvalue.hist = "qqplot",
  
  noFDRsaveMemory = FALSE)

eqtl <- me$all$eqtls
eqtl$gene |> unique() %>% length()
eqtl %>% write_csv("eqtl_pace_nocov.csv")

dfiv <- read_csv("../dfiv.csv")
ivs <- dfiv %>%
  pivot_longer(-snp, names_to = "gene", values_to = "Value") %>%
  filter(Value == 1) %>%
  select(gene, snp) %>%
  mutate(pair=paste(gene,snp,sep="_"))
eqtl <- as.data.table(eqtl)[FDR<0.05,][,pair:=paste(gene,snps,sep="_")]

pre3 <- (eqtl[pair%in%ivs$pair,]|>nrow())/(eqtl|>nrow()) #TDR
recal3 <- (eqtl[pair%in%ivs$pair,]|>nrow())/(ivs|>nrow()) #Recall

### pace_eqtl + 10PCs
useModel = modelLINEAR 
SNP_file_name = paste(base.dir, "/SNP.txt", sep="") 
expression_file_name = paste(base.dir, "/GE_pace.txt", sep="") 
covariates_file_name = character() #
covariates_file_name = paste(base.dir, "/Covariates.txt", sep="")
output_file_name = tempfile() 
pvOutputThreshold = 1 
errorCovariance = numeric() 

snps = SlicedData$new() 
snps$fileDelimiter = " "       
snps$fileOmitCharacters = "NA" 
snps$fileSkipRows = 1        
snps$fileSkipColumns = 1     
snps$fileSliceSize = 2000      
snps$LoadFile( SNP_file_name ) 

gene = SlicedData$new()
gene$fileDelimiter = " " 
gene$fileOmitCharacters = "NA"
gene$fileSkipRows = 1
gene$fileSkipColumns = 1
gene$fileSliceSize = 2000
gene$LoadFile(expression_file_name)

cvrt = SlicedData$new();
cvrt$fileDelimiter = " ";      # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
if(length(covariates_file_name)>0) {
  cvrt$LoadFile(covariates_file_name);
}


maf.list = vector('list', length(snps))
for(sl in 1:length(snps)) {
  slice = snps[[sl]];
  maf.list[[sl]] = rowMeans(slice,na.rm=TRUE)/2;
  maf.list[[sl]] = pmin(maf.list[[sl]],1-maf.list[[sl]]);
}
maf = unlist(maf.list)

snps$RowReorder(maf>=0.01);

me = Matrix_eQTL_engine(  #进行eQTL分析的主要函数
  snps = snps,  #指定SNP 文件
  gene = gene,  #指定基因表达量文件
  cvrt = cvrt,   #指定协变量文件
  output_file_name = output_file_name,  #指定输出文件
  pvOutputThreshold = pvOutputThreshold,  #指定显著性P值
  useModel = useModel,  #指定使用的计算模型
  errorCovariance = errorCovariance, #指定误差项的协方差矩阵
  verbose = TRUE,
  min.pv.by.genesnp = FALSE,  pvalue.hist = "qqplot",
  
  noFDRsaveMemory = FALSE)

eqtl <- me$all$eqtls
eqtl$gene |> unique() %>% length()
eqtl %>% write_csv("eqtl_pace.csv")

dfiv <- read_csv("../dfiv.csv")
ivs <- dfiv %>%
  pivot_longer(-snp, names_to = "gene", values_to = "Value") %>%
  filter(Value == 1) %>%
  select(gene, snp) %>%
  mutate(pair=paste(gene,snp,sep="_"))
eqtl <- as.data.table(eqtl)[FDR<0.05,][,pair:=paste(gene,snps,sep="_")]

pre4 <- (eqtl[pair%in%ivs$pair,]|>nrow())/(eqtl|>nrow()) #TDR
recal4 <- (eqtl[pair%in%ivs$pair,]|>nrow())/(ivs|>nrow()) #Recall

### agg_eqtl
useModel = modelLINEAR 
SNP_file_name = paste(base.dir, "/SNP_sc.txt", sep="") 
expression_file_name = paste(base.dir, "/GE_sc.txt", sep="") 
covariates_file_name = character() #
# covariates_file_name = paste(base.dir, "/Covariates.txt", sep="") 
output_file_name = tempfile() 
pvOutputThreshold = 1 
errorCovariance = numeric() 

snps = SlicedData$new() 
snps$fileDelimiter = " "       
snps$fileOmitCharacters = "NA" 
snps$fileSkipRows = 1        
snps$fileSkipColumns = 1     
snps$fileSliceSize = 2000      
snps$LoadFile( SNP_file_name ) 

gene = SlicedData$new()
gene$fileDelimiter = " " 
gene$fileOmitCharacters = "NA"
gene$fileSkipRows = 1
gene$fileSkipColumns = 1
gene$fileSliceSize = 2000
gene$LoadFile(expression_file_name)

cvrt = SlicedData$new();
cvrt$fileDelimiter = " ";      # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
if(length(covariates_file_name)>0) {
  cvrt$LoadFile(covariates_file_name);
}

maf.list = vector('list', length(snps))
for(sl in 1:length(snps)) {
  slice = snps[[sl]];
  maf.list[[sl]] = rowMeans(slice,na.rm=TRUE)/2;
  maf.list[[sl]] = pmin(maf.list[[sl]],1-maf.list[[sl]]);
}
maf = unlist(maf.list)

snps$RowReorder(maf>=0.01);

me = Matrix_eQTL_engine(  #进行eQTL分析的主要函数
  snps = snps,  #指定SNP 文件
  gene = gene,  #指定基因表达量文件
  cvrt = cvrt,   #指定协变量文件
  output_file_name = output_file_name,  #指定输出文件
  pvOutputThreshold = pvOutputThreshold,  #指定显著性P值
  useModel = useModel,  #指定使用的计算模型
  errorCovariance = errorCovariance, #指定误差项的协方差矩阵
  verbose = TRUE,
  min.pv.by.genesnp = FALSE,  pvalue.hist = "qqplot",
  
  noFDRsaveMemory = FALSE)

eqtl <- me$all$eqtls
eqtl$gene |> unique() %>% length()
eqtl %>% write_csv("eqtl_agg.csv")

dfiv <- read_csv("../dfiv.csv")
ivs <- dfiv %>%
  pivot_longer(-snp, names_to = "gene", values_to = "Value") %>%
  filter(Value == 1) %>%
  select(gene, snp) %>%
  mutate(pair=paste(gene,snp,sep="_"))
eqtl <- as.data.table(eqtl)[FDR<0.05,][,pair:=paste(gene,snps,sep="_")]

pre5 <- (eqtl[pair%in%ivs$pair,]|>nrow())/(eqtl|>nrow()) #TDR
recal5 <- (eqtl[pair%in%ivs$pair,]|>nrow())/(ivs|>nrow()) #Recall

### dyn_eqtl
useModel = modelLINEAR 
SNP_file_name = paste(base.dir, "/SNP_sc.txt", sep="") 
expression_file_name = paste(base.dir, "/GE_sc.txt", sep="") 
covariates_file_name = character() #
covariates_file_name = paste(base.dir, "/Covariates_sc_t.txt", sep="")
output_file_name = tempfile() 
pvOutputThreshold = 1 
errorCovariance = numeric() 

snps = SlicedData$new() 
snps$fileDelimiter = " "       
snps$fileOmitCharacters = "NA" 
snps$fileSkipRows = 1        
snps$fileSkipColumns = 1     
snps$fileSliceSize = 2000      
snps$LoadFile( SNP_file_name ) 

gene = SlicedData$new()
gene$fileDelimiter = " " 
gene$fileOmitCharacters = "NA"
gene$fileSkipRows = 1
gene$fileSkipColumns = 1
gene$fileSliceSize = 2000
gene$LoadFile(expression_file_name)

cvrt = SlicedData$new();
cvrt$fileDelimiter = " ";      # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
if(length(covariates_file_name)>0) {
  cvrt$LoadFile(covariates_file_name);
}


maf.list = vector('list', length(snps))
for(sl in 1:length(snps)) {
  slice = snps[[sl]];
  maf.list[[sl]] = rowMeans(slice,na.rm=TRUE)/2;
  maf.list[[sl]] = pmin(maf.list[[sl]],1-maf.list[[sl]]);
}
maf = unlist(maf.list)

snps$RowReorder(maf>=0.01);

me = Matrix_eQTL_engine(  #进行eQTL分析的主要函数
  snps = snps,  #指定SNP 文件
  gene = gene,  #指定基因表达量文件
  cvrt = cvrt,   #指定协变量文件
  output_file_name = output_file_name,  #指定输出文件
  pvOutputThreshold = pvOutputThreshold,  #指定显著性P值
  useModel = useModel,  #指定使用的计算模型
  errorCovariance = errorCovariance, #指定误差项的协方差矩阵
  verbose = TRUE,
  min.pv.by.genesnp = FALSE,  pvalue.hist = "qqplot",
  
  noFDRsaveMemory = FALSE)

eqtl <- me$all$eqtls
eqtl$gene |> unique() %>% length()
eqtl %>% write_csv("eqtl_dyn.csv")

dfiv <- read_csv("../dfiv.csv")
ivs <- dfiv %>%
  pivot_longer(-snp, names_to = "gene", values_to = "Value") %>%
  filter(Value == 1) %>%
  select(gene, snp) %>%
  mutate(pair=paste(gene,snp,sep="_"))
eqtl <- as.data.table(eqtl)[FDR<0.05,][,pair:=paste(gene,snps,sep="_")]

pre6 <- (eqtl[pair%in%ivs$pair,]|>nrow())/(eqtl|>nrow()) #TDR
recal6 <- (eqtl[pair%in%ivs$pair,]|>nrow())/(ivs|>nrow()) #Recall

### int_eqtl
useModel = modelLINEAR_CROSS 
SNP_file_name = paste(base.dir, "/SNP_sc.txt", sep="") 
expression_file_name = paste(base.dir, "/GE_sc.txt", sep="") 
# covariates_file_name = character() #
covariates_file_name = paste(base.dir, "/Covariates_sc_t.txt", sep="")
output_file_name = tempfile() 
pvOutputThreshold = 1 
errorCovariance = numeric() 

snps = SlicedData$new() 
snps$fileDelimiter = " "       
snps$fileOmitCharacters = "NA" 
snps$fileSkipRows = 1        
snps$fileSkipColumns = 1     
snps$fileSliceSize = 2000      
snps$LoadFile( SNP_file_name ) 

gene = SlicedData$new()
gene$fileDelimiter = " " 
gene$fileOmitCharacters = "NA"
gene$fileSkipRows = 1
gene$fileSkipColumns = 1
gene$fileSliceSize = 2000
gene$LoadFile(expression_file_name)

cvrt = SlicedData$new();
cvrt$fileDelimiter = " ";      # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
if(length(covariates_file_name)>0) {
  cvrt$LoadFile(covariates_file_name);
}


maf.list = vector('list', length(snps))
for(sl in 1:length(snps)) {
  slice = snps[[sl]];
  maf.list[[sl]] = rowMeans(slice,na.rm=TRUE)/2;
  maf.list[[sl]] = pmin(maf.list[[sl]],1-maf.list[[sl]]);
}
maf = unlist(maf.list)

snps$RowReorder(maf>=0.01);

me = Matrix_eQTL_engine(  #进行eQTL分析的主要函数
  snps = snps,  #指定SNP 文件
  gene = gene,  #指定基因表达量文件
  cvrt = cvrt,   #指定协变量文件
  output_file_name = output_file_name,  #指定输出文件
  pvOutputThreshold = pvOutputThreshold,  #指定显著性P值
  useModel = useModel,  #指定使用的计算模型
  errorCovariance = errorCovariance, #指定误差项的协方差矩阵
  verbose = TRUE,
  min.pv.by.genesnp = FALSE,  pvalue.hist = "qqplot",
  noFDRsaveMemory = FALSE)

plot(me)

eqtl <- me$all$eqtls
eqtl$gene |> unique() %>% length()
eqtl %>% write_csv("eqtl_int.csv")

dfiv <- read_csv("../dfiv.csv")
ivs <- dfiv %>%
  pivot_longer(-snp, names_to = "gene", values_to = "Value") %>%
  filter(Value == 1) %>%
  select(gene, snp) %>%
  mutate(pair=paste(gene,snp,sep="_"))
eqtl <- as.data.table(eqtl)[FDR<0.05,][,pair:=paste(gene,snps,sep="_")]

pre7 <- (eqtl[pair%in%ivs$pair,]|>nrow())/(eqtl|>nrow()) #TDR
recal7 <- (eqtl[pair%in%ivs$pair,]|>nrow())/(ivs|>nrow()) #Recall

library(tidyverse)
res <- tibble(Method=c("avg_eqtl","avg_eqtl+10PCs","pace_eqtl",
                       "pace_eqtl+10PCS","agg_eqtl","dyn_eqtl","int_eqtl"),
              Precision=c(pre1,pre2,pre3,pre4,pre5,pre6,pre7),
              Recall=c(recal1,recal2,recal3,recal4,recal5,recal6,recal7))
res %>% write_csv("compare_eqtl.csv")
