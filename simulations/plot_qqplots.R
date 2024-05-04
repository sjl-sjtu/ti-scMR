library(ggplot2)
library(tidyverse)
library(data.table)
library(ggsci)
library(ggthemes)
library(patchwork)

vec1 <- c("sim_genome","true_genome")
vec2 <- c("stable_eqtl","dyn_eqtl")
vec3 <- c("complete_sample","mixed_sample")
combinations <- expand.grid(vec1, vec2, vec3)
path_list <- apply(combinations, 1, paste, collapse = "/")

names1 <- c("Constant effect","Time-varying effect")
names2 <- c("complete sampling","uneven sampling")
names <- apply(expand.grid(names1,names2),1, paste, collapse = "+")

sim_qq <- function(path,name){
  setwd(paste0("/lustre/home/acct-clsyzs/clsyzs/SunJianle/singleCellMR/simulation/n18/",
               path))
  df1 <- fread("eqtl_avg_nocov.csv")[order(pvalue)]
  df2 <- fread("eqtl_avg.csv")[order(pvalue)]
  df3 <- fread("eqtl_pace_nocov.csv")[order(pvalue)]
  df4 <- fread("eqtl_pace.csv")[order(pvalue)]
  df5 <- fread("eqtl_agg.csv")[order(pvalue)]
  df6 <- fread("eqtl_dyn.csv")[order(pvalue)]
  df7 <- fread("eqtl_int.csv")[order(pvalue)]
  q1 <- ppoints(df1$pvalue)
  q2 <- ppoints(df2$pvalue)
  q3 <- ppoints(df3$pvalue)
  q4 <- ppoints(df4$pvalue)
  q5 <- ppoints(df5$pvalue)
  q6 <- ppoints(df6$pvalue)
  q7 <- ppoints(df7$pvalue)
  
  df <- tibble(
    x = -log10(c(q1, q3, q5, q6, q7)),
    y = -log10(c(df1$pvalue, df3$pvalue, df5$pvalue, df6$pvalue, df7$pvalue)),
    group = c(rep("avg_eqtl", length(q1)), rep("cum_eqtl", length(q3)), 
              rep("agg_eqtl", length(q5)),
              rep("dyn_eqtl", length(q6)), rep("int_eqtl", length(q7)))
  ) %>%
    mutate_at(vars(group),factor,levels=c("avg_eqtl","cum_eqtl","agg_eqtl",
                                          "dyn_eqtl","int_eqtl"))
  p <- ggplot(df, aes(x = x, y = y, color = group)) +
    geom_line() +
    # geom_point(shape = 19, size = 0.3) +
    labs(
      x = "-log10(P-value), expected",
      y = "-log10(P-value), observed",
      title = name
    ) +
    ylim(0, 80) +
    scale_x_continuous(breaks = seq(0,6)) +
    geom_abline(intercept = 0, slope = 1, color = "black") +
    scale_color_lancet() +
    theme_clean() +
    theme(legend.position = c(0.05, 0.95),
          legend.justification = c(0, 1),
          legend.title = element_blank(),
          legend.text = element_text(size = 11),
          axis.text = element_text(size = 11),
          axis.title = element_text(size = 12, face = "bold"),
          plot.title = element_text(size = 13, face = "bold.italic"),
          plot.background = element_blank())
}
ql <- map2(path_list[c(1,3,5,7)],names,sim_qq)

(ql[[1]]|ql[[2]])/(ql[[3]]|ql[[4]])



true_qq <- function(path,name){
  setwd(paste0("/lustre/home/acct-clsyzs/clsyzs/SunJianle/singleCellMR/simulation/n18/",
               path))
  df1 <- fread("eqtl_avg_nocov.csv")[order(pvalue)]
  df2 <- fread("eqtl_avg.csv")[order(pvalue)]
  df3 <- fread("eqtl_pace_nocov.csv")[order(pvalue)]
  df4 <- fread("eqtl_pace.csv")[order(pvalue)]
  df5 <- fread("eqtl_agg.csv")[order(pvalue)]
  df6 <- fread("eqtl_dyn.csv")[order(pvalue)]
  df7 <- fread("eqtl_int.csv")[order(pvalue)]
  q1 <- ppoints(df1$pvalue)
  q2 <- ppoints(df2$pvalue)
  q3 <- ppoints(df3$pvalue)
  q4 <- ppoints(df4$pvalue)
  q5 <- ppoints(df5$pvalue)
  q6 <- ppoints(df6$pvalue)
  q7 <- ppoints(df7$pvalue)
  
  df <- tibble(
    x = -log10(c(q2, q4, q5, q6, q7)),
    y = -log10(c(df2$pvalue, df4$pvalue, df5$pvalue, df6$pvalue, df7$pvalue)),
    group = c(rep("avg_eqtl (+PCs)", length(q2)), rep("cum_eqtl (+PCs)", length(q4)), 
              rep("agg_eqtl", length(q5)),
              rep("dyn_eqtl", length(q6)), rep("int_eqtl", length(q7)))
  )%>%
    mutate_at(vars(group),factor,levels=c("avg_eqtl (+PCs)","cum_eqtl (+PCs)","agg_eqtl",
                                          "dyn_eqtl","int_eqtl"))
  p <- ggplot(df, aes(x = x, y = y, color = group)) +
    geom_line() +
    # geom_point(shape = 19, size = 0.3) +
    labs(
      x = "-log10(P-value), expected",
      y = "-log10(P-value), observed",
      title = name
    ) +
    ylim(0, 80) +
    scale_x_continuous(breaks = seq(0,6)) +
    geom_abline(intercept = 0, slope = 1, color = "black") +
    scale_color_lancet() +
    theme_clean() +
    theme(legend.position = c(0.05, 0.95),
          legend.justification = c(0, 1),
          legend.title = element_blank(),
          legend.text = element_text(size = 11),
          axis.text = element_text(size = 11),
          axis.title = element_text(size = 12, face = "bold"),
          plot.title = element_text(size = 13, face = "bold.italic"),
          plot.background = element_blank())
}
ql <- map2(path_list[c(2,4,6,8)],names,true_qq)

(ql[[1]]|ql[[2]])/(ql[[3]]|ql[[4]])
