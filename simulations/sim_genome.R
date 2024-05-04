#
library(tidyverse)
library(data.table)
library(Seurat)
setwd("/lustre/home/acct-clsyzs/clsyzs/SunJianle/singleCellMR/simulation/n18")
p <- 10000 #genome size
n <- 500 #sample size
g <- 100 #gene num
t <- 20 #total time


###########################
# simulated genomes
eaf <- runif(p,0.05,0.95)
dat <- sapply(eaf,
              function(x) sample(0:2,n,replace=T,prob = c((1-x)^2,2*x*(1-x),x^2)))
colnames(dat) <- paste0("snp",1:p)
dat <- as_tibble(dat)
dat$id <- 1:n #paste0("id",1:n)
dat <- dat %>% arrange(id)
dat %>% fwrite("sim_genome/geno.csv",row.names = FALSE)

# real genomes
dat <- fread("data_10000.csv")
dat <- dat[1:n,7:10006]
dat <- dat[, names(dat) := lapply(.SD, as.integer), .SDcols = names(dat)]
colnames(dat) <- paste0("snp",1:p)
dat$id <- 1:n #paste0("id",1:n)
dat %>% fwrite("true_genome/geno.csv",row.names = FALSE)

