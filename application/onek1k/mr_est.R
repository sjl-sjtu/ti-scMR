setwd("~/SunJianle/singleCellMR/GSE196829/expression_matrix/")

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

library(doParallel)
library(foreach)
registerDoParallel(10)
getDoParWorkers()

loaded_packages <- search()[grepl("package:", search())]
loaded_packages <- sub("package:", "", loaded_packages)

cum_mat <- map_dfr(c("pace_cum_bm_bc.csv","pace_cum_pl_bc.csv"),fread)
setnames(cum_mat,"sub","id")
# snps <- fread("../genotype_chr/imputed/merged.raw")
# colnames(snps)[-seq(1,6)] <- sapply(colnames(snps)[-seq(1,6)], function(x) strsplit(x, "_")[[1]][1])

eqtl <- fread("../eqtl/eqtl_all_cis2.csv")%>%filter(pvalue<0.001)
# eqtl1 <- fread("../eqtl/eqtl_all_trans.csv")%>%filter(FDR<0.01)
# eqtl <- rbind(eqtl,eqtl1)
# eqtl1 <- fread("../eqtl/eqtl_bm_cis.csv")
# eqtl2 <- fread("../eqtl/eqtl_pl_cis.csv")
# eqtl1 <- eqtl1[pvalue<0.001,]
# eqtl2 <- eqtl2[pvalue<0.001,]
# eqtl <- rbind(eqtl1,eqtl2)

# eqtl3 <- fread("../eqtl/eqtl_bm_trans.csv")
# eqtl4 <- fread("../eqtl/eqtl_pl_trans.csv")
# eqtl3 <- eqtl3[pvalue<0.001,]
# eqtl4 <- eqtl4[pvalue<0.001,]
# eqtl <- rbind(eqtl1,eqtl2,eqtl3,eqtl4)

# eqtl <- eqtl %>% filter(CELL_TYPE%in%c("Plasma Cell","Naïve/Immature B Cell","Memory B Cell" ))
# eqtl[,snp:=sapply(SNPID, function(x) strsplit(x, "_")[[1]][1])]
# rsid <- eqtl$RSID[match(colnames(snps)[-seq(1,6)], eqtl$snp)]
# colnames(snps)[-seq(1,6)] <- rsid
# snps[,id:=paste(FID,IID,sep="_")]
# ids <- intersect(cum_mat$id,snps$id)
# cum_mat <- cum_mat[id %in% ids]
# snps <- snps[id %in% ids]

snps_all <- fread("../eqtl/SNP.txt")
snps_all <- snps_all[IID%in%unique(eqtl$snps)]
# snps <- fread("../eqtl/SNP_all.txt")
ids <- colnames(snps_all)[-1]


outcome <- "outcome"
genes <- colnames(cum_mat)[-c(ncol(cum_mat)-1,ncol(cum_mat))]
# genelist <- intersect(eqtl$GENE_ID|>unique(),genes)
genelist <- genes
outpath <- "results_b/pace_cis_new4"

stepPrune <- function(df,dfgwas,cutoff){
  snplist <- dfgwas %>% pull(snps)
  candidate <- c()
  while(length(snplist)>0){
    dfgwas <- dfgwas%>%filter(snps %in% snplist)
    j <- dfgwas[which.min(dfgwas$FDR),"snps"] %>% as.character()
    candidate <- c(candidate,j)
    snplist <- setdiff(snplist,j)
    calcLD <- function(i,j,df){
      return(cor(df%>%pull(i),df%>%pull(j))^2)  # Rogers and Huff (2008)
    }
    LDs <- sapply(snplist,calcLD,j,df) 
    snplist <- snplist[which(LDs<cutoff)]
  }
  return(candidate)
}

gc()

MR_pace_all <- function(geneName,pca=F,m=5,k=0){
  # IV <- eqtl %>% filter(GENE_ID==geneName) %>% arrange(FDR) %>% 
  #   filter(pvalue<0.001) %>% pull(RSID)
  tryCatch(
     {
        IV <- eqtl %>% filter(gene==geneName) %>% arrange(FDR) %>% 
          filter(pvalue<0.001) %>% pull(snps) %>% unique()
        j <- length(IV)
        if(j<1){
          return(c(j,rep(NA,4)))
        }
        
        snps <- snps_all[IID %in% IV,][order(match(IID,IV))][,-1]
        snps <- as_tibble(t(snps))
        IV <- make.names(IV)
        colnames(snps) <- IV
        snps$id <- ids

        df <- cum_mat %>%
          left_join(snps %>% dplyr::select(id,any_of(IV)),by="id") %>%
          rename(outcome="outcome")
        
        IV <- df%>%dplyr::select(any_of(IV))%>%colnames()
        s <- apply(df%>%dplyr::select(any_of(IV)),2,function(x) unique(x)|>length())
        IV <- IV[which(s!=1)]
        j <- length(IV)
        if(j<1){
          return(c(j,rep(NA,4)))
        }
        if(j>2){
          dfgwas <- eqtl%>%filter(gene==geneName)%>%mutate(snps=make.names(snps))%>%filter(snps%in%IV)%>%dplyr::select(snps,FDR)
          IV <- stepPrune(df,dfgwas,0.1)
        }
        
        if(pca==T&j>1){
          fit_prcomp <- prcomp(df%>%dplyr::select(any_of(IV))%>%as.matrix(), center = T, scale = T)
          IVprom <- as_tibble(fit_prcomp$x)
          colnames(IVprom) <- paste0("IVPC",1:ncol(IVprom))
          df <- df %>% bind_cols(IVprom)
          m1 <- min(m,j-1)
          IVpc <- paste0("IVPC",1:m1)
        }else{
          IVpc <- IV
        }
        
        
        if(k>0){
          fit_prcomp <- prcomp(df%>%dplyr::select(setdiff(genelist,geneName)), center = T, scale = T)
          genepc <- paste0("PC",1:k)
          df <- df %>% bind_cols(fit_prcomp$x%>%as_tibble() %>%dplyr::select(any_of(genepc)))
          formu1 <- as.formula(paste0(geneName,"~",paste(c(IVpc,genepc),collapse="+")))
        }else{
          formu1 <- as.formula(paste0(geneName,"~",paste(IVpc,collapse="+")))
        }
      #   print(formu1)
        
        lm1 <- lm(formu1,data=df)
        df$res <- as.numeric(residuals(lm1))
        
        if(k==0){
          formu2 <- as.formula(paste0("outcome ~",geneName,"+ res"))
        }else{
          formu2 <- as.formula(paste0("outcome ~",geneName,"+ res +", paste(genepc,collapse = "+")))  
        }
      #   print(formu2)
        
        lm2 <- glm(formu2,data=df,family=binomial)
        
        
        Fstat <- unname(summary(lm1)$fstatistic[1])
        
        expWalpha <- fitted.values(lm1)
        expXbeta <- fitted.values(lm2)
        
        W<-cbind(1,df%>%dplyr::select(any_of(IVpc))%>%as.matrix())  #IVpc
        
        if(k==0){
          X<-cbind(1,df%>%pull(geneName),df$res) %>% as.matrix()
        }else{
          X<-cbind(1,df%>%pull(geneName),df$res,df%>%dplyr::select(any_of(genepc))) %>% as.matrix()
        }
        
        beta <- coef(lm2)
        
        if(NA%in%beta){
          return(c(j,Fstat,rep(NA,3)))
        }
        
        Xbeta <- X%*%beta
        bxu<-beta[3]
        
        paJ <- -bxu*expXbeta*expWalpha*W 
        pbJ <- expXbeta*X 
        Bba <- t(pbJ)%*%paJ 
        Bbb <- t(pbJ)%*%pbJ 
        
        
        covalpha <- vcovHC(lm1,type="HC1") 
        covbeta <- vcovHC(lm2,type="HC1") 
        d22<-MASS::ginv(Bbb)%*%Bba%*%covalpha%*%t(Bba)%*%MASS::ginv(Bbb)+covbeta 
        ses<-sqrt(diag(d22))
        b <- beta[2]
        se <- ses[2]
        
        testre <- wald.test(Sigma=d22, b=beta, Terms=2)
        p <- unname(testre$result$chi2["P"])
        return(c(j,Fstat,b,se,p))
     },error=function(e){
      return(rep(NA,5))
     }
  )

  
}


all_re <- foreach(geneName = genelist, .combine = 'rbind',.packages = loaded_packages) %dopar% {
  MR_pace_all(geneName, pca=F, m=0, k=0)
}
colnames(all_re) <- c("IVnum","F","beta","se","p")
all_re <- as_tibble(all_re)
all_re$gene <- genelist

all_re %>% write_csv(paste0(outpath,"/re_pace_all_",outcome,".csv"))



MR_pace_linear <- function(geneName,pca=F,m=5,k=0){
  # IV <- eqtl %>% filter(GENE_ID==geneName) %>% arrange(FDR) %>% 
  #   filter(pvalue<0.001) %>% pull(RSID)
  tryCatch(
    {
      IV <- eqtl %>% filter(gene==geneName) %>% arrange(FDR) %>% 
        filter(pvalue<0.001) %>% pull(snps) %>% unique()
      j <- length(IV)
      
      if(j<1){
        return(c(j,rep(NA,4)))
      }

      snps <- snps_all[IID %in% IV,][order(match(IID,IV))][,-1]
      snps <- as_tibble(t(snps))
      IV <- make.names(IV)
      colnames(snps) <- IV
      snps$id <- ids
      
      df <- cum_mat %>%
        left_join(snps %>% dplyr::select(id,any_of(IV)),by="id") %>%
        rename(outcome="outcome")

      IV <- df%>%dplyr::select(any_of(IV))%>%colnames()
      s <- apply(df%>%dplyr::select(any_of(IV)),2,function(x) unique(x)|>length())
      IV <- IV[which(s!=1)]
      j <- length(IV)
      if(j<1){
        return(c(j,rep(NA,4)))
      }
      if(j>2){
        dfgwas <- eqtl%>%filter(gene==geneName)%>%mutate(snps=make.names(snps))%>%filter(snps%in%IV)%>%dplyr::select(snps,FDR)
        IV <- stepPrune(df,dfgwas,0.1)
      }
      
      if(pca==T&j>1){
        fit_prcomp <- prcomp(df%>%dplyr::select(any_of(IV)), center = T, scale = T)
        IVprom <- as_tibble(fit_prcomp$x)
        colnames(IVprom) <- paste0("IVPC",1:ncol(IVprom))
        df <- df %>% bind_cols(IVprom)
        m1 <- min(m-1,j-1)
        IVpc <- paste0("IVPC",1:m1)
      }else{
        IVpc <- IV
      }
      
      if(k>0){
        fit_prcomp <- prcomp(df%>%dplyr::select(setdiff(genelist,geneName)), center = T, scale = T)
        genepc <- paste0("PC",1:k)
        df <- df %>% bind_cols(fit_prcomp$x[,genepc])
      }
      
      if(k==0){
        formu <- as.formula(paste("outcome ~",geneName,"|",paste(IVpc,collapse = "+")))
      }else{
        formu <- as.formula(paste("outcome ~",geneName,"+", paste(genepc,collapse = "+"), "|",paste(genepc,collapse = "+"),"+",paste(IVpc,collapse = "+")))
      }
      
      ivlm <- AER::ivreg(formu,data = df)
      
      res <- summary(ivlm, vcov = sandwich,diagnostics = TRUE)
      Fstat <- res$diagnostics[1,3]
      b <- res$coefficients[2,1]
      se <- res$coefficients[2,2]
      p <- res$coefficients[2,4]
      
      return(c(j,Fstat,b,se,p))
    },error=function(e){
      return(rep(NA,5))
    }
  )
}


all_linear <- foreach(geneName = genelist,
                      .combine = 'rbind',.packages = loaded_packages) %dopar% {
                        MR_pace_linear(geneName, pca=F, m=0,  k=0)
                      }
colnames(all_linear) <- c("IVnum","F","beta","se","p")
all_linear <- as_tibble(all_linear)
all_linear$gene <- genelist
all_linear %>% write_csv(paste0(outpath,"/re_pace_linear_",outcome,".csv"))


MR_pace_lasso <- function(geneName,m=5,repeats=100){
  # IV <- eqtl %>% filter(GENE_ID==geneName) %>% arrange(FDR) %>% filter(pvalue<0.001) %>% pull(RSID)
  
  tryCatch(
    {

  IV <- eqtl %>% filter(gene==geneName) %>% arrange(FDR) %>% 
    filter(pvalue<0.001) %>% pull(snps) %>% unique()
  j <- length(IV)
  if(j<1){
    return(c(j,rep(NA,6)))
  }

  snps <- snps_all[IID %in% IV,][order(match(IID,IV))][,-1]
  snps <- as_tibble(t(snps))
  IV <- make.names(IV)
  colnames(snps) <- IV
  snps$id <- ids
  
  df <- cum_mat %>%
    left_join(snps %>% dplyr::select(id,any_of(IV)),by="id") %>%
    dplyr::rename(xcum=any_of(geneName),outcome=any_of(outcome))
  
  IV <- df%>%dplyr::select(any_of(IV))%>%colnames()
  s <- apply(df%>%dplyr::select(any_of(IV)),2,function(x) unique(x)|>length())
  IV <- IV[which(s!=1)]
  j <- length(IV)
  if(j<1){
    return(c(j,rep(NA,6)))
  }
  if(j>2){
    dfgwas <- eqtl%>%filter(gene==geneName)%>%mutate(snps=make.names(snps))%>%filter(snps%in%IV)%>%dplyr::select(snps,FDR)
    IV <- stepPrune(df,dfgwas,0.1)
  }
  j <- length(IV)
  
  formu <- as.formula(paste("xcum ~",paste(IV,collapse="+")))
  lm1 <- lm(formu,data=df)
  
  Fstat <- unname(summary(lm1)$fstatistic[1])
  
  df$res <- residuals(lm1)
  y <- df$outcome
  X <- as.matrix(df%>%dplyr::select(xcum,res,any_of(IV)))
  p_fac <- rep(1,ncol(X))
  p_fac[1:2] <- 0
  lm2 <- cv.glmnet(X, y,  
                   family="binomial", 
                   intercept = T, alpha=1,
                   penalty.factor=p_fac)
  lambda <- lm2$lambda.min
  
  getlasso <- function(df,ind){
    dfn <- df[ind,]
    lm1 <- lm(formu,data=dfn)
    dfn$res <- residuals(lm1)
    y <- dfn$outcome
    X <- as.matrix(dfn%>%dplyr::select(xcum,res,any_of(IV)))
    p_fac <- rep(1,ncol(X))
    p_fac[1:2] <- 0
    lm2 <- glmnet(X, y,  
                  family='binomial', lambda=lambda,
                  intercept = T, alpha=1,
                  penalty.factor=p_fac) 
    beta <- coef(lm2)[2]
    return(beta)
  }
  
  samplers <- boot(df,getlasso,R=repeats,parallel = "multicore")
  b <- samplers$t0
  se <- summary(samplers)$bootSE
  cis <- boot.ci(samplers,type = "perc")
  ci_l <- cis$percent[4]
  ci_u <- cis$percent[5]
  pvalues <- boot.pval(samplers,type = "perc")
  
  return(c(j,Fstat,b,ci_l,ci_u,se,pvalues))

    },error=function(e){
      return(rep(NA,7))
    }
  )
  
}  


res_mrl <- foreach(geneName = genelist,.packages = loaded_packages) %dopar% {
  MR_pace_lasso(geneName, m = timeLength, repeats = 1000)
}
res_mrl <- purrr::reduce(res_mrl,rbind)
colnames(res_mrl) <- c("IVnum","F","beta","ci_l","ci_u","se","p")
res_mrl <- as_tibble(res_mrl)
res_mrl$gene <- genelist
res_mrl %>% write_csv(paste0(outpath,"/re_pace_lasso_",outcome,".csv"))



MR_pace_lasso_linear <- function(geneName,repeats=100){
  # IV <- eqtl %>% filter(GENE_ID==geneName) %>% arrange(FDR) %>% filter(pvalue<0.001) %>% pull(RSID)

  tryCatch(
    {
  IV <- eqtl %>% filter(gene==geneName) %>% arrange(FDR) %>% 
    filter(pvalue<0.001) %>% pull(snps) %>% unique()
  j <- length(IV)
  if(j<1){
    return(c(j,rep(NA,6)))
  }

  snps <- snps_all[IID %in% IV,][order(match(IID,IV))][,-1]
  snps <- as_tibble(t(snps))
  IV <- make.names(IV)
  colnames(snps) <- IV
  snps$id <- ids
  
  df <- cum_mat %>%
    left_join(snps %>% dplyr::select(id,any_of(IV)),by="id") %>%
    dplyr::rename(xcum=any_of(geneName),outcome=any_of(outcome))
  
  IV <- df%>%dplyr::select(any_of(IV))%>%colnames()
  s <- apply(df%>%dplyr::select(any_of(IV)),2,function(x) unique(x)|>length())
  IV <- IV[which(s!=1)]
  j <- length(IV)
  if(j<1){
    return(c(j,rep(NA,6)))
  }
  if(j>2){
    dfgwas <- eqtl%>%filter(gene==geneName)%>%mutate(snps=make.names(snps))%>%filter(snps%in%IV)%>%dplyr::select(snps,FDR)
    IV <- stepPrune(df,dfgwas,0.1)
  }
  j <- length(IV)
  
  formu <- as.formula(paste("xcum ~",paste(IV,collapse="+")))
  lm1 <- lm(formu,data=df)
  
  Fstat <- unname(summary(lm1)$fstatistic[1])
  
  df$res <- residuals(lm1)
  y <- df$outcome
  X <- as.matrix(df%>%dplyr::select(xcum,res,any_of(IV)))
  p_fac <- rep(1,ncol(X))
  p_fac[1:2] <- 0
  lm2 <- cv.glmnet(X, y,  
                   family="gaussian", 
                   intercept = T, alpha=1,
                   penalty.factor=p_fac)
  lambda <- lm2$lambda.min
  
  getlasso <- function(df,ind){
    dfn <- df[ind,]
    lm1 <- lm(formu,data=dfn)
    dfn$res <- residuals(lm1)
    y <- dfn$outcome
    X <- as.matrix(dfn%>%dplyr::select(xcum,res,any_of(IV)))
    p_fac <- rep(1,ncol(X))
    p_fac[1:2] <- 0
    lm2 <- glmnet(X, y,  
                  family="gaussian", lambda = lambda,
                  intercept = T, alpha=1,
                  penalty.factor=p_fac) 
    beta <- coef(lm2)[2]
    return(beta)
  }
  
  samplers <- boot(df,getlasso,R=repeats,parallel = "multicore")
  b <- samplers$t0
  se <- summary(samplers)$bootSE
  cis <- boot.ci(samplers,type = "perc")
  ci_l <- cis$percent[4]
  ci_u <- cis$percent[5]
  pvalues <- boot.pval(samplers,type = "perc")
  
  return(c(j,Fstat,b,ci_l,ci_u,se,pvalues))
    },error=function(e){
      return(rep(NA,7))
    }
  )

}  

res_mrl <- foreach(geneName = genelist,.packages = loaded_packages) %dopar% {
  MR_pace_lasso_linear(geneName, repeats = 1000)
}
res_mrl <- purrr::reduce(res_mrl,rbind)
colnames(res_mrl) <- c("IVnum","F","beta","ci_l","ci_u","se","p")
res_mrl <- as_tibble(res_mrl)
res_mrl$gene <- genelist
res_mrl %>% write_csv(paste0(outpath,"/re_pace_lasso_linear_",outcome,".csv"))
