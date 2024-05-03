library(tidyverse)
library(ggthemes)
library(ggsci)
df <- read_csv("C:\\Users\\lenovo\\Desktop\\sim_eqtl.csv")
df$Method <- factor(df$Method,
                    levels = c("avg_eqtl","cum_eqtl","agg_eqtl","dyn_eqtl","int_eqtl"))
ggplot(df, aes(x = Method)) +
  geom_line(aes(y = Precision, color = "Precision (TDR)"),group=1,size=1) +
  geom_point(aes(y = Precision, color = "Precision (TDR)"), size = 1.5) +
  geom_line(aes(y = Recall, color = "Recall (Power)"),group=1,size=1) +
  geom_point(aes(y = Recall, color = "Recall (Power)"), size = 1.5) +
  facet_wrap(~ Scenario, nrow = 1) +
  labs(x = "Methods", y = "Values",color="") +
  scale_color_jama(alpha=0.7)+
  theme_bw()+
  theme(legend.position = "top", 
        axis.text.x = element_text(size=11,angle = 45, hjust = 1),
        axis.title = element_text(size = 12,face = "bold"),
        strip.text = element_text(size = 11,face = "bold.italic"),
        legend.text = element_text(size = 11))+
  ylim(c(0,1))

df1 <- df %>% pivot_longer(-c(Scenario,Method),
                           names_to = "Index",values_to = "Value") %>%
  mutate(Value_label = sprintf("%.2f", Value))
df1 <- df1 %>% mutate(Value = replace_na(Value, 0))
ggplot(df1, aes(x = Method,y = Value,fill = Index)) +
  geom_bar(stat = "identity", position = "dodge") + 
  geom_text(aes(label = Value_label), vjust = -0.5, size = 3, 
            position = position_dodge(width = 1)) +
  facet_wrap(~ Scenario, ncol = 1) +
  labs(x = "Methods", y = "Values") +
  scale_fill_nejm(alpha=0.7)+
  theme_bw()+
  theme(legend.position = "right", 
        axis.text.x = element_text(size=11,angle = 45, hjust = 1),
        axis.title = element_text(size = 12,face = "bold"),
        strip.text = element_text(size = 11,face = "bold.italic"),
        legend.text = element_text(size = 11))+
  ylim(c(0,1))


df <- read_csv("C:\\Users\\lenovo\\Desktop\\true_eqtl.csv")
df$Method <- factor(df$Method,
                    levels = c("avg_eqtl","avg_eqtl+10PCs","cum_eqtl",
                               "cum_eqtl+10PCs","agg_eqtl","dyn_eqtl","int_eqtl"))
ggplot(df, aes(x = Method)) +
  geom_line(aes(y = Precision, color = "Precision (TDR)"),group=1,size=1) +
  geom_point(aes(y = Precision, color = "Precision (TDR)"), size = 1.5) +
  geom_line(aes(y = Recall, color = "Recall (Power)"),group=1,size=1) +
  geom_point(aes(y = Recall, color = "Recall (Power)"), size = 1.5) +
  facet_wrap(~ Scenario, nrow = 1) +
  labs(x = "Methods", y = "Values",color="") +
  scale_color_jama(alpha=0.7)+
  theme_bw()+
  theme(legend.position = "top", 
        axis.text.x = element_text(size=11,angle = 45, hjust = 1),
        axis.title = element_text(size = 12,face = "bold"),
        strip.text = element_text(size = 11,face = "bold.italic"),
        legend.text = element_text(size = 11))+
  ylim(c(0,1))

df1 <- df %>% pivot_longer(-c(Scenario,Method),
                           names_to = "Index",values_to = "Value") %>%
  mutate(Value_label = sprintf("%.2f", Value))
df1 <- df1 %>% mutate(Value = replace_na(Value, 0))
ggplot(df1, aes(x = Method,y = Value,fill = Index)) +
  geom_bar(stat = "identity", position = "dodge") + 
  geom_text(aes(label = Value_label), vjust = -0.5, size = 3, 
            position = position_dodge(width = 1)) +
  facet_wrap(~ Scenario, ncol = 1) +
  labs(x = "Methods", y = "Values") +
  scale_fill_nejm(alpha=0.7)+
  theme_bw()+
  theme(legend.position = "right", 
        axis.text.x = element_text(size=11,angle = 45, hjust = 1),
        axis.title = element_text(size = 12,face = "bold"),
        strip.text = element_text(size = 11,face = "bold.italic"),
        legend.text = element_text(size = 11))+
  ylim(c(0,1))
