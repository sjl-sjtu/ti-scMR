library(tidyverse)
library(ggthemes)
library(ggsci)

vec1 <- c("sim_genome","true_genome")
vec2 <- c("stable_eqtl","dyn_eqtl")
vec3 <- c("complete_sample","mixed_sample")
combinations <- expand.grid(vec1, vec2, vec3)
path_list <- apply(combinations, 1, paste, collapse = "/")

names1 <- c("Constant effect","Time-varying effect")
names2 <- c("complete sampling","uneven sampling")
names <- apply(expand.grid(names1,names2),1, paste, collapse = "+")

getdf <- function(path,name){
  setwd(paste0("/lustre/home/acct-clsyzs/clsyzs/SunJianle/singleCellMR/simulation/n18/",
               path))
  df <- read_csv("compare_eqtl.csv")
  df <- df %>% mutate(Method = case_when(
    Method == "pace_eqtl" ~ "cum_eqtl",
    Method == "pace_eqtl+10PCS" ~ "cum_eqtl+10PCs",
    TRUE ~ Method
  ))
  df$Method <- factor(df$Method,
                      levels = c("avg_eqtl","avg_eqtl+10PCs","cum_eqtl","cum_eqtl+10PCs",
                                 "agg_eqtl","dyn_eqtl","int_eqtl"))
  df$Scenario <- name
  return(df)
}
df_sim <- map2_dfr(path_list[c(1,3,5,7)],names,getdf) %>% 
  filter(Method %in% c("avg_eqtl","cum_eqtl",
         "agg_eqtl","dyn_eqtl","int_eqtl")) %>%
  mutate_at(vars(Scenario),factor,levels=names)
df_true <- map2_dfr(path_list[c(2,4,6,8)],names,getdf) %>%
  mutate_at(vars(Scenario),factor,levels=names)

df1 <- df_sim %>% pivot_longer(-c(Scenario,Method),
                           names_to = "Index",values_to = "Value") %>%
  mutate(Value_label = sprintf("%.2f", Value))
df1 <- df1 %>% mutate(Value = replace_na(Value, 0))
ggplot(df1, aes(x = Method,y = Value,fill = Index)) +
  geom_bar(stat = "identity", position = "dodge") + 
  geom_text(aes(label = Value_label), vjust = -0.5, size = 3, 
            position = position_dodge(width = 1)) +
  facet_wrap(~ Scenario,nrow=2, ncol = 2) +
  labs(x = "Methods", y = "Values") +
  scale_fill_nejm(alpha=0.7)+
  theme_bw()+
  theme(legend.position = "right", 
        axis.text.x = element_text(size=11,angle = 45, hjust = 1),
        axis.title = element_text(size = 12,face = "bold"),
        strip.text = element_text(size = 11,face = "bold.italic"),
        legend.text = element_text(size = 11))+
  ylim(c(0,1))


df1 <- df_true %>% pivot_longer(-c(Scenario,Method),
                               names_to = "Index",values_to = "Value") %>%
  mutate(Value_label = sprintf("%.2f", Value))
df1 <- df1 %>% mutate(Value = replace_na(Value, 0))
ggplot(df1, aes(x = Method,y = Value,fill = Index)) +
  geom_bar(stat = "identity", position = "dodge") + 
  geom_text(aes(label = Value_label), vjust = -0.5, size = 3, 
            position = position_dodge(width = 1)) +
  facet_wrap(~ Scenario,nrow=2, ncol = 2) +
  labs(x = "Methods", y = "Values") +
  scale_fill_nejm(alpha=0.7)+
  theme_bw()+
  theme(legend.position = "right", 
        axis.text.x = element_text(size=11,angle = 45, hjust = 1),
        axis.title = element_text(size = 12,face = "bold"),
        strip.text = element_text(size = 11,face = "bold.italic"),
        legend.text = element_text(size = 11))+
  ylim(c(0,1))
