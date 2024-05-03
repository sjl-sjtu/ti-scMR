library(MatrixEQTL)

base.dir = getwd()
useModel = modelLINEAR #三种模型可选(modelANOVA or modelLINEAR or modelLINEAR_CROSS)
SNP_file_name = paste(base.dir, "/SNP.txt", sep="") #获取SNP文件位置
expression_file_name = paste(base.dir, "/GE_pace.txt", sep="") #获取基因表达量文件位置
# expression_file = data.table::fread(expression_file_name, header=T) #读取基因表达量文件，可以在R中查看
covariates_file_name = character() #
covariates_file_name = paste(base.dir, "/Covariates.txt", sep="") #读取协变量文件 #
output_file_name = tempfile() # 设置输出文件
pvOutputThreshold = 1 #定义gene-SNP associations的显著性P值
errorCovariance = numeric() #定义误差项的协方差矩阵

#加载基因型文件
snps = SlicedData$new() #创建SNP文件为S4对象（S4对象是该包的默认输入方式）
snps$fileDelimiter = " "       #指定SNP文件的分隔符
snps$fileOmitCharacters = "NA" #定义缺失值
snps$fileSkipRows = 1        #跳过第一行（适用于第一行是列名的情况）
snps$fileSkipColumns = 1     #跳过第一列（适用于第一列是SNP ID的情况）
snps$fileSliceSize = 2000      #每次读取2000条数据
snps$LoadFile( SNP_file_name ) #载入SNP文件

#加载基因表达文件
gene = SlicedData$new()
gene$fileDelimiter = " " 
gene$fileOmitCharacters = "NA"
gene$fileSkipRows = 1
gene$fileSkipColumns = 1
gene$fileSliceSize = 2000
gene$LoadFile(expression_file_name)

#加载协变量
cvrt = SlicedData$new();
cvrt$fileDelimiter = " ";      # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
if(length(covariates_file_name)>0) {
  cvrt$LoadFile(covariates_file_name);
}
#文件的输入部分结束


maf.list = vector('list', length(snps))
for(sl in 1:length(snps)) {
  slice = snps[[sl]];
  maf.list[[sl]] = rowMeans(slice,na.rm=TRUE)/2;
  maf.list[[sl]] = pmin(maf.list[[sl]],1-maf.list[[sl]]);
}
maf = unlist(maf.list)

## Look at the distribution of MAF
# hist(maf[maf<0.1],seq(0,0.1,length.out=100))

cat('SNPs before filtering:',nrow(snps))
# snps$RowReorderSimple(maf>0.05);
snps$RowReorder(maf>=0.01);
cat('SNPs after filtering:',nrow(snps))


#运行分析
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

#结果
cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n'); #查看分析所用时间
cat('Detected eQTLs:', '\n'); 
plot(me)
eqtl <- me$all$eqtls
eqtl$gene |> unique() %>% length()

eqtl %>% write_csv("eqtl_pace.csv")

# 准确率
dfiv <- read_csv("../dfiv.csv")
ivs <- dfiv %>%
  pivot_longer(-snp, names_to = "gene", values_to = "Value") %>%
  filter(Value == 1) %>%
  select(gene, snp) %>%
  mutate(pair=paste(gene,snp,sep="_"))

eqtl <- as.data.table(eqtl)[FDR<0.05,][,pair:=paste(gene,snps,sep="_")]

(eqtl[pair%in%ivs$pair,]|>nrow())/(eqtl|>nrow()) #Precision
(eqtl[pair%in%ivs$pair,]|>nrow())/(ivs|>nrow()) #Recall
