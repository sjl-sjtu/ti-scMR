cum_expression <- function(exposureDat,genelist,id_var="id",pseudotime_var="pseudotime"){
  library(tidyverse)
  library(doParallel)
  library(foreach)
  library(fdapace)
  library(pracma)
  registerDoParallel()
  getDoParWorkers()
  loaded_packages <- search()[grepl("package:", search())]
  loaded_packages <- sub("package:", "", loaded_packages)

  exposureDat <- exposureDat %>%
    group_by(id,pseudotime) %>%
    summarise(across(genelist),mean)) %>%
    ungroup() %>%
    drop_na() %>% 
    dplyr::rename(id=get(id_var),pseudotime=get(pseudotime_var))
  
  fullPACE <- function(geneName){
    s <- exposureDat %>% 
      dplyr::select(id,pseudotime,any_of(geneName)) %>%
      dplyr::rename(expression=any_of(geneName)) %>%
      group_by(id)  %>% 
      arrange(pseudotime) %>%
      group_map(~.x)
    N <- length(s)
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
      cums <- trapz(x=Ti_est,y=x_it[i,])
      return(cums)
    }
    xcum <- map_dbl(1:length(s),getCurrEffe)
    return(xcum)
  }

  cum_mat <- foreach(geneName = genelist, .combine = 'cbind',.packages = loaded_packages) %dopar% {
    fullPACE(geneName)
  }

  colnames(cum_mat) <- genelist
  cum_mat <- as_tibble(cum_mat)
  cum_mat$id <- exposureDat %>% group_by(id) %>% group_keys() %>% pull(id)
  cum_mat <- cum_mat %>% left_join(exposureDat %>% dplyr::select(id) %>% 
                                     dplyr::distinct(id,.keep_all = T),by="id")
  return(cum_mat)
}

