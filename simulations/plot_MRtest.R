library(tidyverse)
library(ggsci)

setwd("C:\\Users\\lenovo\\Desktop")

df <- read_csv("quant_rep_mixed_power.csv")
df1 <- df %>% select(-repeats) %>% group_by(outcome) %>% 
  summarise_at(vars(everything()),mean)
df1 <- df1 %>% pivot_longer(cols = -outcome, names_to = "method", values_to = "power")
df2 <- df %>% select(-repeats) %>% group_by(outcome) %>% 
  summarise_at(vars(everything()),function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))))
df2 <- df2 %>% pivot_longer(cols = -outcome, names_to = "method", values_to = "se")
df1 <- df1%>%left_join(df2)
df1$method <- factor(df1$method,levels = c("avg_linear","pace_linear","avg_linear_lasso","pace_linear_lasso"))
# ggplot(df1)+geom_point(aes(outcome,power,color=method))
# ggplot(df1)+geom_point(aes(method,power,color=outcome))
p1 <- ggplot(df1,aes(outcome,power,fill=method))+
  geom_errorbar(aes(ymin = power - se, ymax = power + se), 
                width = 0.2, position = position_dodge(0.9))+
  geom_col(position = "dodge")+
  scale_fill_aaas()+
  ylab("mean power")+
  xlab("outcome")
dfp <- df1

df <- read_csv("quant_rep_mixed_FDR.csv")
df1 <- df %>% select(-repeats) %>% group_by(outcome) %>% summarise_at(vars(everything()),mean,na.rm=T)
df1 <- df1 %>% pivot_longer(cols = -outcome, names_to = "method", values_to = "FDR")
df2 <- df %>% select(-repeats) %>% group_by(outcome) %>% summarise_at(vars(everything()),function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))))
df2 <- df2 %>% pivot_longer(cols = -outcome, names_to = "method", values_to = "se")
df1 <- df1%>%left_join(df2)
# df1$FDR[df1$FDR==0] <- 0.001
df1$method <- factor(df1$method,levels = c("avg_linear","pace_linear","avg_linear_lasso","pace_linear_lasso"))
p2 <- ggplot(df1,aes(outcome,FDR,fill=method))+
  geom_errorbar(aes(ymin = FDR - se, ymax = FDR + se), 
                width = 0.2, position = position_dodge(0.9))+
  geom_col(position = "dodge")+
  scale_fill_aaas()+
  ylab("mean FDR")+
  xlab("outcome")
dff <- df1

df <- left_join(dff,dfp,by=c("outcome","method"))
df %>% write_csv("quant.csv")

df <- read_csv("binary_rep_mixed_power.csv")
df1 <- df %>% select(-repeats) %>% group_by(outcome) %>% summarise_at(vars(everything()),mean)
df1 <- df1 %>% pivot_longer(cols = -outcome, names_to = "method", values_to = "power")
df2 <- df %>% select(-repeats) %>% group_by(outcome) %>% summarise_at(vars(everything()),function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))))
df2 <- df2 %>% pivot_longer(cols = -outcome, names_to = "method", values_to = "se")
df1 <- df1%>%left_join(df2)
df1$method <- factor(df1$method,levels = c("avg_linear","pace_linear",
                                           "pace_linear_lasso","avg_logit",
                                           "pace_logit","pace_logit_lasso"))
p3 <- ggplot(df1,aes(outcome,power,fill=method))+
  geom_errorbar(aes(ymin = power - se, ymax = power + se), 
                width = 0.2, position = position_dodge(0.9))+
  geom_col(position = "dodge")+
  scale_fill_aaas()+
  ylab("mean power")+
  xlab("outcome")
dfp <- df1

df <- read_csv("binary_rep_mixed_FDR.csv")
df1 <- df %>% select(-repeats) %>% group_by(outcome) %>% summarise_at(vars(everything()),mean,na.rm=T)
df1 <- df1 %>% pivot_longer(cols = -outcome, names_to = "method", values_to = "FDR")
df2 <- df %>% select(-repeats) %>% group_by(outcome) %>% summarise_at(vars(everything()),function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))))
df2 <- df2 %>% pivot_longer(cols = -outcome, names_to = "method", values_to = "se")
df1 <- df1%>%left_join(df2)
# df1$FDR[df1$FDR==0] <- 0.001
df1$method <- factor(df1$method,levels = c("avg_linear","pace_linear",
                                           "pace_linear_lasso","avg_logit",
                                           "pace_logit","pace_logit_lasso"))
p4 <- ggplot(df1,aes(outcome,FDR,fill=method))+
  geom_errorbar(aes(ymin = FDR - se, ymax = FDR + se), 
                width = 0.2, position = position_dodge(0.9))+
  geom_col(position = "dodge")+
  scale_fill_aaas()+
  ylab("mean FDR")+
  xlab("outcome")
dff <- df1

df <- left_join(dff,dfp,by=c("outcome","method"))
df %>% write_csv("binary.csv")

library(patchwork)
((p2+p1)/(p4+p3))+plot_annotation(tag_levels = "a")
