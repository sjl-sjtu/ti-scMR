library(tidyverse)
library(data.table)
library(ggsci)
pal <- "lancet"
alpha=1
n=5
pal <- getFromNamespace(paste0("pal_", pal), "ggsci")
colors <- pal(alpha = alpha)(n)



# df1 <- fread("eqtl_avg_nocov.csv")[order(pvalue)]#[-nrow(df1),]
df2 <- fread("eqtl_avg.csv")[order(pvalue)]#[-nrow(df2),]
# df3 <- fread("eqtl_pace_nocov.csv")[order(pvalue)]#[-nrow(df3),]
df4 <- fread("eqtl_pace.csv")[order(pvalue)]#[-nrow(df4),]
df5 <- fread("eqtl_sc.csv")[order(pvalue)]#[-nrow(df5),]
df6 <- fread("eqtl_sc_t_nointer.csv")[order(pvalue)]#[-nrow(df6),]
df7 <- fread("eqtl_sc_t.csv")[order(pvalue)]#[-nrow(df7),]
q1 <- ppoints(df1$pvalue)
q2 <- ppoints(df2$pvalue)
q3 <- ppoints(df3$pvalue)
q4 <- ppoints(df4$pvalue)
q5 <- ppoints(df5$pvalue)
q6 <- ppoints(df6$pvalue)
q7 <- ppoints(df7$pvalue)

plot.new()
# plot(-log10(q1),-log10(df1$pvalue), type="s",
#      pch=19, col=colors[1],cex=0.3,xlab="-log10(P-value), expected",
#      ylab="-log10(P-value), observed",ylim=c(0,70))
plot(-log10(q2),-log10(df2$pvalue), type="s",
     pch=19, col=colors[1],cex=0.3,xlab="-log10(P-value), expected",
     ylab="-log10(P-value), observed",ylim=c(0,12),xlim=c(0,3.2))
# points(-log10(q3),-log10(df3$pvalue), type="s",
#      pch=19, col=colors[3],cex=0.3,xlab="-log10(P-value), expected",
#      ylab="-log10(P-value), observed",ylim=c(0,70))
points(-log10(q4),-log10(df4$pvalue), type="s",
     pch=19, col=colors[2],cex=0.3,xlab="-log10(P-value), expected",
     ylab="-log10(P-value), observed")#,ylim=c(0,50),xlim=c(0,6))
points(-log10(q5),-log10(df5$pvalue), type="s",
     pch=19, col=colors[3],cex=0.3,xlab="-log10(P-value), expected",
     ylab="-log10(P-value), observed")#,ylim=c(0,50),xlim=c(0,6))
points(-log10(q6),-log10(df6$pvalue), type="s",
     pch=19, col=colors[4],cex=0.3,xlab="-log10(P-value), expected",
     ylab="-log10(P-value), observed")#,ylim=c(0,50),xlim=c(0,6))
abline(a=0, b=1,col="black")

# par(xpd = TRUE)  # 允许图形绘制超出边界
legend("topleft", legend = c(#"avg_eqtl", 
                             "avg_eqtl (+PCs)", 
                             # "cum_eqtl",
                              "cum_eqtl (+PCs)", "dyn_eqtl", "dyn_int_eqtl"), 
       col = colors, pch = 16, inset = c(0, 0), bty = "n", ncol=1,xpd=T,
       xjust = 1, yjust = 1)
# par(xpd = FALSE)  # 恢复图形绘制边界限制


plot.new()
# plot(-log10(q1),-log10(df1$pvalue), type="s",
#      pch=19, col=colors[1],cex=0.3,xlab="-log10(P-value), expected",
#      ylab="-log10(P-value), observed",ylim=c(0,70))
plot(-log10(q2),-log10(df2$pvalue), type="l",
     pch=19, col=colors[1],cex=0.3,xlab="-log10(P-value), expected",
     ylab="-log10(P-value), observed",ylim=c(0,80))#,xlim=c(0,6))
# points(-log10(q3),-log10(df3$pvalue), type="s",
#      pch=19, col=colors[3],cex=0.3,xlab="-log10(P-value), expected",
#      ylab="-log10(P-value), observed",ylim=c(0,70))
points(-log10(q4),-log10(df4$pvalue), type="l",
       pch=19, col=colors[2],cex=0.3,xlab="-log10(P-value), expected",
       ylab="-log10(P-value), observed")#,ylim=c(0,50),xlim=c(0,6))
points(-log10(q5),-log10(df5$pvalue), type="l",
       pch=19, col=colors[3],cex=0.3,xlab="-log10(P-value), expected",
       ylab="-log10(P-value), observed")#,ylim=c(0,50),xlim=c(0,6))
points(-log10(q6),-log10(df6$pvalue), type="l",
       pch=19, col=colors[4],cex=0.3,xlab="-log10(P-value), expected",
       ylab="-log10(P-value), observed")#,ylim=c(0,50),xlim=c(0,6))
points(-log10(q7),-log10(df7$pvalue), type="l",
       pch=19, col=colors[5],cex=0.3,xlab="-log10(P-value), expected",
       ylab="-log10(P-value), observed")#,ylim=c(0,50),xlim=c(0,6))
abline(a=0, b=1,col="black")

# par(xpd = TRUE)  # 允许图形绘制超出边界
legend("topleft", legend = c(
  "avg_eqtl",
  # "avg_eqtl (+PCs)",
  "cum_eqtl",
  # "cum_eqtl (+PCs)",
  "agg_eqtl",
  "dyn_eqtl", yuan
  "int_eqtl"), 
  col = colors, pch = 16, inset = c(0, 0), bty = "n", ncol=1,xpd=T,
  xjust = 1, yjust = 1)
# par(xpd = FALSE)  # 恢复图形绘制边界限制

