#1a
pdf(paste0(cancer,"-volcano_top30.pdf"), width = 10)
p <- ggscatter(deg.data, x = "logFC", y = "log10P",
               color = "Group",
               palette = c("#218359", "#BBBBBB", "#BF3330"),
               alpha = 0.8, size = 2.5,
               label = deg.data$Label,
               
               font.label = 12, repel = TRUE,
               xlab = "logFC", ylab = "-log10(padj)") +
  theme_base() +
  geom_hline(yintercept = 1.30, linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  xlim(c(-10, 10)) +
  labs(x = "logFC", y = "-log10(padj)") +
  scale_x_continuous(breaks = c(-10,-8,-6,-4,-2,-1,-0.5,0,0.5,1,2,4,6,8,10), expand = c(0.01, 0))
p
dev.off()
##############################################################################################
# 绘制小提琴图
outTab <- data.frame()
pdf("20PCDscore_TN_Diff_vioplot.pdf", width = 10, height = 7)
par(las = 1.5, mar = c(15, 6, 3, 3))  # 下边距调整到15以容纳标签
x <- c(1:ncol(rt_log2))
y <- c(1:ncol(rt_log2))
plot(x, y,
     xlim = c(1, ncol(rt_log2) * 3), ylim = c(min(rt_log2), max(rt_log2) + 0.02),  #x,y轴标签
     main = "", xlab = "", ylab = "Fraction",
     pch = 21,
     col = "white",
     xaxt = "n")

# 对每个免疫检查点循环，绘制vioplot，低表达用蓝色表示，高表达用红色表示
bioCol <- c("#7CC767", "#FF0000")
for (i in 1:ncol(rt_log2)) {
  if (sd(rt_log2[1:lowNum, i]) == 0) {
    rt_log2[1, i] <- 0.00001
  }
  if (sd(rt_log2[(lowNum + 1):(lowNum + highNum), i]) == 0) {
    rt_log2[(lowNum + 1), i] <- 0.00001
  }
  lowData <- rt_log2[1:lowNum, i]
  highData <- rt_log2[(lowNum + 1):(lowNum + highNum), i]
  vioplot(lowData, at = 3* (i - 1), lty = 1, add = TRUE, col = bioCol[1])  #间距？
  vioplot(highData, at = 3* (i - 1) + 1, lty = 1, add = TRUE, col = bioCol[2])
  wilcoxTest <- wilcox.test(lowData, highData)
  p <- wilcoxTest$p.value
  cellPvalue <- cbind(Cell = colnames(rt_log2)[i], pvalue = p)
  outTab <- rbind(outTab, cellPvalue)
  mx <- max(c(lowData, highData))
  lines(c(3 * (i - 1) + 0.2, 3 * (i - 1) + 0.8), c(mx, mx))
  text(x = 3 * (i - 1) + 0.5, y = mx + 0.02, labels = ifelse(p < 0.001, "p<0.001", paste0("p=", sprintf("%.03f", p))), cex = 0.8)
}

# 绘制图例
legend("topright",
       legend = c("Normal", "Tumor"),
       lwd = 4.5, bty = "n", cex = 1.5,
       col = bioCol)

# 绘制X轴标签并对齐
axis(1, at = seq(0, ncol(rt_log2) * 3 - 1, 3) + 0.5, labels = colnames(rt_log2), las = 2, cex.axis = 0.9)

dev.off()





##############################################################################################################################

#绘制boxplot
p=ggboxplot(data, x="Gene", y="Expression", color = "Type", 
            ylab="PCD-score",
            xlab="",
            legend.title=x,
            palette = c("#11AA4D","#D20A13"),
            width=0.6, add = "none")
p=p+rotate_x_text(60)
p1=p+stat_compare_means(aes(group=Type),
                        method="wilcox.test",
                        symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                        label = "p.signif")

#输出图片
pdf(file=outFile, width=6, height=7)
print(p1)
dev.off()





#################################################################################################################

library(pheatmap)

pdf(paste("clinical_Diff_all.pdf",sep=''),5,length(rownames(df))/2)

paletteLength = 1000
# #immune
# myColor <- colorRampPalette(c("white", "#FF7C00"))(paletteLength)
# #exp
# myColor <- colorRampPalette(c("white", "red"))(paletteLength)
# #cell
# myColor <- colorRampPalette(c("white","blue"))(paletteLength)
# #drug
# myColor <- colorRampPalette(c("white", "#660BAB"))(paletteLength)
# #yzx_gx
# myColor <- colorRampPalette(c("white", "#C7007D"))(paletteLength)

max(df)

#bk <- c(seq( -max(abs(max(df)),abs(min(df))), -0.1,by=0.01),seq(0,max(abs(max(df)),abs(min(df))),by=0.01))
# bk <- c(seq( -25, -0.1,by=0.01),seq(0,25,by=0.01))
# myColor <- c(colorRampPalette(colors = c("blue","white"))(length(bk)/2),colorRampPalette(colors = c("white","red"))(length(bk)/2))

bk <- seq(0,8, by = 0.01) #,"#ffd5dc"
myColor <- colorRampPalette(colors = c("white", "red"))(length(bk))


#myBreaks <- c(seq( min(df),-min(df), length.out=floor(paletteLength/2)))

#bk <- c(seq(-9,-0.1,by=0.01),seq(0,9,by=0.01))

###
getSig <- function(dc) {
  sc <- ' '
  if (dc < 0.0001) {sc <- '****'}
  else if (dc < 0.001){sc <- '***'}
  else if (dc < 0.01){sc <- '**'}
  else if (dc < 0.05) {sc <- '*'}
  else{sc <- ''}
  return(sc)
}

sig.mat <- matrix(sapply(as.matrix(p_value), getSig), nrow=nrow(as.matrix(p_value)))
str(sig.mat)
###
xx <- pheatmap(df,
               color=myColor,
               breaks=bk,
               clustering_method="average", cluster_rows=F,cluster_cols=F, cellwidth = 20,cellheight = 20,main="-log10(p)",display_numbers=sig.mat)
#print(xx)


###########################################################################################################

my_colors <- c("white", "#00A087B2","#BB362F")

pdf('Entosis_4surv_heatmap_12.6.pdf',width=20,height=20)
High1 <- pheatmap(final_,
                  border_color='black',
                  fontsize_row = 15,
                  fontsize_col = 15,
                  color = my_colors,
                  fontsize_number = 10,
                  cluster_rows=F,
                  cluster_cols=F,
                  cellwidth = 25,
                  cellheight = 25)
dev.off()




dev.off()

