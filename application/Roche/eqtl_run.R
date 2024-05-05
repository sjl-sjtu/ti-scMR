library(MatrixEQTL)
library(data.table)

base.dir = getwd()
useModel = modelLINEAR #三种模型可选(modelANOVA or modelLINEAR or modelLINEAR_CROSS),这里我们选择modelLINEAR
SNP_file_name = paste(base.dir, "/SNP.txt", sep="") #获取SNP文件位置

expression_file_name = paste(base.dir, "/GE.txt", sep="") #获取基因表达量文件位置

# covariates_file_name = character() #
covariates_file_name = paste(base.dir, "/Covariates.txt", sep="") #读取协变量文件 #
# covariates_file = data.table::fread(covariates_file_name, header=T) #读取协变量文件，可在R中查看

snps_location_file_name = paste(base.dir, "/snp_pos.txt", sep="")
gene_location_file_name = paste(base.dir, "/gene_pos.txt", sep="")

# Output file name
output_file_name_cis = tempfile()
output_file_name_tra = tempfile()

# Only associations significant at this level will be saved
pvOutputThreshold_cis = 1e-3;
pvOutputThreshold_tra = 1e-4;
errorCovariance = numeric() #定义误差项的协方差矩阵 (用的很少)
cisDist = 1e6

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

## Run the analysis
snpspos = read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE);
genepos = read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE);

# maf.list = vector('list', length(snps))
# for(sl in 1:length(snps)) {
#   slice = snps[[sl]];
#   maf.list[[sl]] = rowMeans(slice,na.rm=TRUE)/2;
#   maf.list[[sl]] = pmin(maf.list[[sl]],1-maf.list[[sl]]);
# }
# maf = unlist(maf.list)


# cat('SNPs before filtering:',nrow(snps))
# snps$RowReorder(maf>=0.05);
# cat('SNPs after filtering:',nrow(snps))


## Run the analysis
me = Matrix_eQTL_main(
  snps = snps,
  gene = gene,
  cvrt = cvrt,
  output_file_name     = output_file_name_tra,
  pvOutputThreshold     = pvOutputThreshold_tra,
  useModel = useModel,
  errorCovariance = errorCovariance,
  verbose = TRUE,
  output_file_name.cis = output_file_name_cis,
  pvOutputThreshold.cis = pvOutputThreshold_cis,
  snpspos = snpspos,
  genepos = genepos,
  cisDist = cisDist,
  pvalue.hist = "qqplot",
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE)

plot(me)

unlink(output_file_name_tra);
unlink(output_file_name_cis);

## Results:

cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');
cat('Detected local eQTLs:', '\n');
# show(me$cis$eqtls)
cat('Detected distant eQTLs:', '\n');
# show(me$trans$eqtls)
eqtl_cis <- me$cis$eqtls
eqtl_trans <- me$trans$eqtls

eqtl_cis |> fwrite("eqtl_cis.csv")
eqtl_trans |> fwrite("eqtl_trans.csv")

