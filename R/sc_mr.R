sc_mr <- function(df,geneName,outcome,IV,method="linear",id_var="id",covs=NULL){
  library(tidyverse)
  library(data.table)
  library(sandwich)
  library(aod)
  library(ivreg)
  library(glmnet)
  library(boot)
  library(boot.pval)
  
  s <- apply(df%>%dplyr::select(any_of(IV)),2,function(x) unique(x)|>length())
  IV <- IV[which(s!=1)]
  j <- length(IV)
  if(j<1){
    return(c(j,rep(NA,4)))
  }
  df <- df %>% dplyr::rename(outcome=get(outcome),id=get(id_var)) %>%
             dplyr::select(id,outcome,any_of(IV),any_of(geneName),any_of(covs))
             
  MR_pace_logit <- function(geneName){
    tryCatch(
      {
        formu1 <- as.formula(paste(geneName,"~",paste(IV,collapse="+")))
        lm1 <- lm(formu1,data=df)
        df$res <- as.numeric(residuals(lm1))
        if(is.null(covs)){
          formu2 <- as.formula(paste("outcome ~",geneName,"+ res"))
        }else{
          formu2 <- as.formula(paste("outcome ~",geneName,"+ res +", paste(covs,collapse = "+")))
        }
        lm2 <- glm(formu2,data=df,family=binomial)
        Fstat <- unname(summary(lm1)$fstatistic[1])
        expWalpha <- fitted.values(lm1)
        expXbeta <- fitted.values(lm2)
        W<-cbind(1,df%>%dplyr::select(any_of(IV))%>%as.matrix())  #IV
        if(is.null(covs)){
          X<-cbind(1,df%>%pull(geneName),df$res) %>% as.matrix()
        }else{
          X<-cbind(1,df%>%pull(geneName),df$res,df%>%select(any_of(covs))) %>% as.matrix()
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
  
  MR_pace_linear <- function(geneName){
    tryCatch(
      {
        if(is.null(covs)){
          formu <- as.formula(paste("outcome ~",geneName,"|",paste(IV,collapse = "+")))
        }else{
          formu <- as.formula(paste("outcome ~",geneName,"+", paste(covs,collapse = "+"), "|",paste(covs,collapse = "+"),"+",paste(IV,collapse = "+")))
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

  MR_pace_lasso_logit <- function(geneName,repeats=5000){
    tryCatch(
      {   
        formu <- as.formula(paste(geneName,"~",paste(IV,collapse="+")))
        lm1 <- lm(formu,data=df)
        Fstat <- unname(summary(lm1)$fstatistic[1])
        df$res <- residuals(lm1)
        y <- df$outcome
        X <- as.matrix(df%>%select(any_of(geneName),res,any_of(IV)))
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
          X <- as.matrix(dfn%>%select(any_of(geneName),res,any_of(IV)))
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
        formu <- as.formula(paste(geneName,"~",paste(IV,collapse="+")))
        lm1 <- lm(formu,data=df)
        Fstat <- unname(summary(lm1)$fstatistic[1])
        df$res <- residuals(lm1)
        y <- df$outcome
        X <- as.matrix(df%>%select(any_of(geneName),res,any_of(IV)))
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
          X <- as.matrix(dfn%>%select(any_of(geneName),res,any_of(IV)))
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

  if(method=="linear"){
    res <- MR_pace_linear(geneName)
    names(res) <- c("IVnum","F","beta","se","p")
  }else if(method="logit"){
    res <- MR_pace_logit(geneName)
    names(res) <- c("IVnum","F","beta","se","p")
  }else if(method=="linear_lasso"){
    res <- MR_pace_lasso_linear(geneName)
    names(res) <- c("IVnum","F","beta","ci_l","ci_u","se","p")
  }else if(method=="logit_lasso"){
    res <- MR_pace_lasso_logit(geneName)
    names(res) <- c("IVnum","F","beta","ci_l","ci_u","se","p")
  }else{
    warning("No method!")
    res <- NULL
  }
  return(res)
}
