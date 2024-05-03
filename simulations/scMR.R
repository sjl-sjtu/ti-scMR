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

MR_pace_logit <- function(geneName,pca=F,m=5,k=0){
  tryCatch(
    {
      IV <- eqtl %>% filter(gene==geneName) %>% arrange(pvalue) %>%
        filter(FDR<0.05) %>% pull(snps)
      j <- length(IV)
      
      df <- cum_mat %>%
        left_join(snps %>% dplyr::select(id,any_of(IV)),by="id") %>%
        rename(outcome=any_of(outcome)) # %>%
        # mutate_at(vars(any_of(genelist)),function(x) scale(x)|>as.numeric())
      
      if(j<1){
        return(c(j,rep(NA,4)))
      }
      s <- apply(df%>%dplyr::select(any_of(IV)),2,function(x) unique(x)|>length())
      IV <- IV[which(s!=1)]
      j <- length(IV)
      if(j<1){
        return(c(j,rep(NA,4)))
      }
      if(j>2){
        dfgwas <- eqtl%>%filter(gene==geneName)%>%mutate(snps=make.names(snps))%>%filter(snps%in%IV)%>%select(snps,FDR)
        IV <- stepPrune(df,dfgwas,0.1)
      }
      j <- length(IV)
      
      if(pca==T&j>1){
        fit_prcomp <- prcomp(df%>%dplyr::select(any_of(IV))%>%as.matrix(), center = T, scale = T)
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
        df <- df %>% bind_cols(fit_prcomp$x%>%as_tibble() %>%select(any_of(genepc)))
      }
      
      formu1 <- as.formula(paste(geneName,"~",paste(IVpc,collapse="+")))
      lm1 <- lm(formu1,data=df)
      df$res <- as.numeric(residuals(lm1))
      
      if(k==0){
        formu2 <- as.formula(paste("outcome ~",geneName,"+ res"))
      }else{
        formu2 <- as.formula(paste("outcome ~",geneName,"+ res +", paste(genepc,collapse = "+")))
      }
      
      lm2 <- glm(formu2,data=df,family=binomial)
      
      
      Fstat <- unname(summary(lm1)$fstatistic[1])
      
      expWalpha <- fitted.values(lm1)
      expXbeta <- fitted.values(lm2)
      
      W<-cbind(1,df%>%dplyr::select(any_of(IVpc))%>%as.matrix())  #IVpc
      
      if(k==0){
        X<-cbind(1,df%>%pull(geneName),df$res) %>% as.matrix()
      }else{
        X<-cbind(1,df%>%pull(geneName),df$res,df%>%select(any_of(genepc))) %>% as.matrix()
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

MR_pace_linear <- function(geneName,pca=F,m=5,k=0){
  tryCatch(
    {
      IV <- eqtl %>% filter(gene==geneName) %>% arrange(pvalue) %>%
        filter(FDR<0.05) %>% pull(snps)
      j <- length(IV)
      
      df <- cum_mat %>%
        left_join(snps %>% dplyr::select(id,any_of(IV)),by="id") %>%
        rename(outcome=any_of(outcome)) #%>% 
        #mutate_at(vars(any_of(genelist)),function(x) scale(x)|>as.numeric())
      
      if(j<1){
        return(c(j,rep(NA,4)))
      }
      s <- apply(df%>%dplyr::select(any_of(IV)),2,function(x) unique(x)|>length())
      IV <- IV[which(s!=1)]
      j <- length(IV)
      if(j<1){
        return(c(j,rep(NA,4)))
      }
      if(j>2){
        dfgwas <- eqtl%>%filter(gene==geneName)%>%mutate(snps=make.names(snps))%>%filter(snps%in%IV)%>%select(snps,FDR)
        IV <- stepPrune(df,dfgwas,0.1)
      }
      j <- length(IV)
      
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
      
      ivlm <- ivreg(formu,data = df)
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

MR_pace_lasso_logit <- function(geneName,m=5,repeats=5000){
  tryCatch(
    {
      IV <- eqtl %>% filter(gene==geneName) %>% arrange(pvalue) %>% filter(FDR<0.05) %>% pull(snps)
      j <- length(IV)
      
      df <- cum_mat %>%
        left_join(snps %>% dplyr::select(id,any_of(IV)),by="id") %>%
        dplyr::rename(xcum=any_of(geneName),outcome=any_of(outcome)) # %>%
        # mutate_at(vars(any_of(genelist)),function(x) scale(x)|>as.numeric())
      
      if(j<1){
        return(c(j,rep(NA,6)))
      }
      s <- apply(df%>%dplyr::select(any_of(IV)),2,function(x) unique(x)|>length())
      IV <- IV[which(s!=1)]
      j <- length(IV)
      if(j<1){
        return(c(j,rep(NA,6)))
      }
      if(j>2){
        dfgwas <- eqtl%>%filter(gene==geneName)%>%mutate(snps=make.names(snps))%>%filter(snps%in%IV)%>%select(snps,FDR)
        IV <- stepPrune(df,dfgwas,0.1)
      }
      j <- length(IV)
      
      
      formu <- as.formula(paste("xcum ~",paste(IV,collapse="+")))
      lm1 <- lm(formu,data=df)
      
      Fstat <- unname(summary(lm1)$fstatistic[1])
      
      df$res <- residuals(lm1)
      y <- df$outcome
      X <- as.matrix(df%>%select(xcum,res,any_of(IV)))
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
        X <- as.matrix(dfn%>%select(xcum,res,any_of(IV)))
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

MR_pace_lasso_linear <- function(geneName,repeats=5000){
  tryCatch(
    {
      IV <- eqtl %>% filter(gene==geneName) %>% arrange(pvalue) %>% filter(FDR<0.05) %>% pull(snps)
      j <- length(IV)
      
      df <- cum_mat %>%
        left_join(snps %>% dplyr::select(id,any_of(IV)),by="id") %>%
        dplyr::rename(xcum=any_of(geneName),outcome=any_of(outcome)) # %>%
        # mutate_at(vars(any_of(genelist)),function(x) scale(x)|>as.numeric())
      
      if(j<1){
        return(c(j,rep(NA,6)))
      }
      s <- apply(df%>%dplyr::select(any_of(IV)),2,function(x) unique(x)|>length())
      IV <- IV[which(s!=1)]
      j <- length(IV)
      if(j<1){
        return(c(j,rep(NA,6)))
      }
      if(j>2){
        dfgwas <- eqtl%>%filter(gene==geneName)%>%mutate(snps=make.names(snps))%>%filter(snps%in%IV)%>%select(snps,FDR)
        IV <- stepPrune(df,dfgwas,0.1)
      }
      j <- length(IV)
      
      formu <- as.formula(paste("xcum ~",paste(IV,collapse="+")))
      lm1 <- lm(formu,data=df)
      
      Fstat <- unname(summary(lm1)$fstatistic[1])
      
      df$res <- residuals(lm1)
      y <- df$outcome
      X <- as.matrix(df%>%select(xcum,res,any_of(IV)))
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
        X <- as.matrix(dfn%>%select(xcum,res,any_of(IV)))
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

mr_logit <- function(geneName,pca=F,m=5,k=0){
  tryCatch(
    {
      IV <- eqtl %>% filter(gene==geneName) %>% arrange(pvalue) %>%
        filter(FDR<0.05) %>% pull(snps)
      j <- length(IV)
      
      if(j<1){
        return(c(j,rep(NA,4)))
      }
      
      df <- ge %>% select(id,any_of(genelist)) %>%
        group_by(id) %>% summarise_at(vars(any_of(genelist)),mean) %>% ungroup() %>%
        left_join(snps%>%select(any_of(IV),id),by="id") %>%
        left_join(out%>%select(any_of(outcome),id),by="id") %>%
        select(-id) %>%
        rename(outcome=any_of(outcome)) # %>%
        # mutate_at(vars(any_of(genelist)),function(x) scale(x)|>as.numeric())
      
      s <- apply(df%>%dplyr::select(any_of(IV)),2,function(x) unique(x)|>length())
      IV <- IV[which(s!=1)]
      j <- length(IV)
      if(j==0){
        return(c(j,rep(NA,4)))
      }
      if(j>2){
        dfgwas <- eqtl%>%filter(gene==geneName)%>%mutate(snps=make.names(snps))%>%filter(snps%in%IV)%>%select(snps,FDR)
        IV <- stepPrune(df,dfgwas,0.1)
      }
      j <- length(IV)
      
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
      
      
      fulldf <- df %>% select(any_of(geneName),outcome,any_of(IVpc)) %>% rename(expression = any_of(geneName))
      
      if(k>0){
        fit_prcomp <- prcomp(df%>%dplyr::select(setdiff(genelist,geneName)), center = T, scale = T)
        genepc <- paste0("PC",1:k)
        fulldf <- fulldf %>% bind_cols(fit_prcomp$x[,genepc])
      }
      
      
      lm1 <- lm(as.formula(paste("expression~",paste(IVpc,collapse = "+"))),data=fulldf)
      fulldf$res <- residuals(lm1)
      if(k>0){
        lm2 <- glm(as.formula(paste("outcome~expression+res+",paste(genepc,collapse = "+"))),data=fulldf,family=binomial)
      }else{
        lm2 <- glm(outcome~expression+res,data=fulldf,family=binomial)
      }
      
      Fstat <- unname(summary(lm1)$fstatistic[1])
      
      expWalpha <- fitted.values(lm1)
      expXbeta <- fitted.values(lm2)
      
      W<-cbind(1,fulldf%>%dplyr::select(any_of(IVpc))%>%as.matrix())
      if(k>0){
        X<-cbind(1,fulldf%>%pull(expression),fulldf$res,fit_prcomp$x[,genepc])
      }else{
        X<-cbind(1,fulldf%>%pull(expression),fulldf$res)
      }
      
      beta <- coef(lm2)
      if(NA %in% beta){
        return(rep(NA,3))
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

mr_linear <- function(geneName,pca=F,m=5,k=0){
  tryCatch(
    {
      IV <- eqtl %>% filter(gene==geneName) %>% arrange(pvalue) %>% 
        filter(FDR<0.05) %>% pull(snps)
      j <- length(IV)
      
      if(j<1){
        return(c(j,rep(NA,4)))
      }
      
      df <- ge %>% select(id,any_of(genelist)) %>% group_by(id) %>% 
        summarise_at(vars(any_of(genelist)),mean) %>% ungroup() %>%
        left_join(snps%>%select(any_of(IV),id),by="id") %>%
        left_join(out%>%select(any_of(outcome),id),by="id") %>%
        select(-id) %>%
        rename(outcome=any_of(outcome)) # %>%
        # mutate_at(vars(any_of(genelist)),function(x) scale(x)|>as.numeric())
      
      s <- apply(df%>%dplyr::select(any_of(IV)),2,function(x) unique(x)|>length())
      IV <- IV[which(s!=1)]
      j <- length(IV)
      if(j==0){
        return(c(j,rep(NA,4)))
      }
      if(j>2){
        dfgwas <- eqtl%>%filter(gene==geneName)%>%mutate(snps=make.names(snps))%>%filter(snps%in%IV)%>%select(snps,FDR)
        IV <- stepPrune(df,dfgwas,0.1)
      }
      j <- length(IV)
      
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
        df <- df %>% bind_cols(fit_prcomp$x)
        genepc <- paste0("PC",1:k)
        formu <- as.formula(paste("outcome ~",geneName,"+", paste(genepc,collapse = "+"), "|",paste(genepc,collapse = "+"),"+",paste(IVpc,collapse = "+")))
      }else{
        formu <- as.formula(paste("outcome ~ ", geneName, "|",paste(IVpc,collapse = "+")))
      }
      
      ivlm <- ivreg(formu,data = df)
      res <- summary(ivlm, vcov = sandwich,diagnostics=T)
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

MR_lasso_logit <- function(geneName,repeats=5000){
  tryCatch(
    {
      IV <- eqtl %>% filter(gene==geneName) %>% arrange(pvalue) %>% filter(FDR<0.05) %>% pull(snps)
      j <- length(IV)
      
      df <- ge %>% select(id,any_of(genelist)) %>% group_by(id) %>%
        summarise_at(vars(any_of(genelist)),mean) %>% ungroup() %>%
        left_join(snps%>%select(any_of(IV),id),by="id") %>%
        left_join(out%>%select(any_of(outcome),id),by="id") %>%
        select(-id) %>%
        rename(xcum=any_of(geneName),outcome=any_of(outcome)) # %>%
        # mutate_at(vars(any_of(genelist)),function(x) scale(x)|>as.numeric())
      
      if(j<1){
        return(c(j,rep(NA,6)))
      }
      s <- apply(df%>%dplyr::select(any_of(IV)),2,function(x) unique(x)|>length())
      IV <- IV[which(s!=1)]
      j <- length(IV)
      if(j<1){
        return(c(j,rep(NA,6)))
      }
      if(j>2){
        dfgwas <- eqtl%>%filter(gene==geneName)%>%mutate(snps=make.names(snps))%>%filter(snps%in%IV)%>%select(snps,FDR)
        IV <- stepPrune(df,dfgwas,0.1)
      }
      j <- length(IV)
      
      
      formu <- as.formula(paste("xcum ~",paste(IV,collapse="+")))
      lm1 <- lm(formu,data=df)
      
      Fstat <- unname(summary(lm1)$fstatistic[1])
      
      df$res <- residuals(lm1)
      y <- df$outcome
      X <- as.matrix(df%>%select(xcum,res,any_of(IV)))
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
        X <- as.matrix(dfn%>%select(xcum,res,any_of(IV)))
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

MR_lasso_linear <- function(geneName,repeats=5000){
  tryCatch(
    {
      IV <- eqtl %>% filter(gene==geneName) %>% arrange(pvalue) %>% filter(FDR<0.05) %>% pull(snps)
      j <- length(IV)
      
      df <- ge %>% select(id,any_of(genelist)) %>% group_by(id) %>%
        summarise_at(vars(any_of(genelist)),mean) %>% ungroup() %>%
        left_join(snps%>%select(any_of(IV),id),by="id") %>%
        left_join(out%>%select(any_of(outcome),id),by="id") %>%
        select(-id) %>%
        rename(xcum=any_of(geneName),outcome=any_of(outcome)) # %>%
        # mutate_at(vars(any_of(genelist)),function(x) scale(x)|>as.numeric())
      
      if(j<1){
        return(c(j,rep(NA,6)))
      }
      s <- apply(df%>%dplyr::select(any_of(IV)),2,function(x) unique(x)|>length())
      IV <- IV[which(s!=1)]
      j <- length(IV)
      if(j<1){
        return(c(j,rep(NA,6)))
      }
      if(j>2){
        dfgwas <- eqtl%>%filter(gene==geneName)%>%mutate(snps=make.names(snps))%>%filter(snps%in%IV)%>%select(snps,FDR)
        IV <- stepPrune(df,dfgwas,0.1)
      }
      j <- length(IV)
      
      
      formu <- as.formula(paste("xcum ~",paste(IV,collapse="+")))
      lm1 <- lm(formu,data=df)
      
      Fstat <- unname(summary(lm1)$fstatistic[1])
      
      df$res <- residuals(lm1)
      y <- df$outcome
      X <- as.matrix(df%>%select(xcum,res,any_of(IV)))
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
        X <- as.matrix(dfn%>%select(xcum,res,any_of(IV)))
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
