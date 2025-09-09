pdf(paste0(cancer,"-HH-LL-volcano.pdf"), width = 10)
p <- ggscatter(deg.data, x = "logFC", y = "log10P",
               color = "Group",
               palette = c("#2f5688", "#BBBBBB", "#CC0000"),
               alpha = 0.8, size = 2.5,
               # label = deg.data$Label,
               font.label = 12, repel = TRUE,
               xlab = "logFC", ylab = "-log10(padj)") +
  theme_base() +
  geom_hline(yintercept = 1.30, linetype = "dashed") +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed") +
  xlim(c(-10, 10)) +
  labs(x = "logFC", y = "-log10(padj)") +
  scale_x_continuous(breaks = c(-10,-8,-6,-4,-2,-1,-0.5,0,0.5,1,2,4,6,8,10), expand = c(0.01, 0))
p
dev.off()

###################################################################

#绘制boxplot
p=ggboxplot(data, x="PCDscore", y="Expression", color = "group", 
            ylab="PCD-score",
            xlab="",
            legend.title=x,
            palette = c("#A60000","#009999"),
            width=0.6, add = "none")
p=p+rotate_x_text(60)
p1=p+stat_compare_means(aes(group=group),
                        method="wilcox.test",
                        symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                        label = "p.signif")

#输出图片
pdf(file=outFile, width=7, height=3.5)
print(p1)
dev.off()
##########################################################3
pdf(paste("cluster_Tumordryness.pdf",sep=''),width=8,height = 6)
p <- ggboxplot(dat, x = "Gene", y = "value",
               color = "Group", palette = c("#A60000","#009999"), width=0.6,
               add = "jitter") # palette可以按照期刊选择相应的配色，如"npg"等
p <- p+xlab("Type")+ylab("Expression Value")
p <- p + stat_compare_means(aes(group = Group), label = "p.signif",method = "anova")
# p <- p + scale_x_discrete(limits=c("hot","cold")
p

print(p)
dev.off()
#####################################################3
pdf(paste("cor_genedown_Tumordryness.pdf",sep=''),width =3.1,height = 8)
pheatmap(data.r, 
         color=myColor,
         breaks=myBreaks,
         clustering_method="average",cluster_rows = F,cluster_cols=F,display_numbers=sig.mat)

dev.off()
#####################################################33
# 肿瘤微环境差异分析
for(i in colnames(data)[1:(ncol(data)-1)]) {
  # 绘制箱线图
  boxplot <- ggboxplot(data, x = "Cluster", y = i, fill = "Cluster",
                       xlab = "",
                       ylab = i,
                       legend.title = "Cluster",
                       palette = bioCol
  ) + 
    stat_compare_means(comparisons = my_comparisons) +
    geom_jitter(aes(color = Cluster), width = 0.2)
  
  # 输出图形
  pdf(file = paste0(i, ".pdf"), width = 4, height = 4.5)
  print(boxplot)
  dev.off()  
}
#######################################################33
myColor <- colorRampPalette(c("#52ACAC", "white", "#FD5C5C"))(paletteLength)

test <- data.r
# myBreaks <- c(
#   seq(-0.7, 0, length.out = floor(paletteLength/2) + 1)
# )
# myBreaks <- c(seq(min(test), 0, length.out=ceiling(paletteLength/2) + 1),
#               seq(max(test)/paletteLength, max(test), length.out=floor(paletteLength/2)))
myBreaks <- c(seq(-0.3, 0, length.out=ceiling(paletteLength/2) + 1),
              seq(0.8/paletteLength, 0.6, length.out=floor(paletteLength/2)))


#chr [1:6, 1:5] "*" "***" "" "***" "***" "***" "***" "" "***" "**" ...

pdf(paste("cor_genedown_High_cell.pdf",sep=''),width =2.5,height = 7)
pheatmap(data.r, 
         color=myColor,
         breaks=myBreaks,
         clustering_method="average", cluster_rows=F, cluster_cols=F,display_numbers=sig.mat
)

dev.off()
##################################################################3
bioCol=c("#009999","#A60000","#FF9900","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
for(i in 1:ncol(rt)){
  if(sd(rt[1:lowNum,i])==0){
    rt[1,i]=0.00001
  }
  if(sd(rt[(lowNum+1):(lowNum+highNum),i])==0){
    rt[(lowNum+1),i]=0.00001
  }
  lowData=rt[1:lowNum,i]
  highData=rt[(lowNum+1):(lowNum+highNum),i]
  vioplot(lowData,at=3*(i-1),lty=1,add = T,col=bioCol[1])
  vioplot(highData,at=3*(i-1)+1,lty=1,add = T,col=bioCol[2])
  wilcoxTest=wilcox.test(lowData,highData)
  p=wilcoxTest$p.value
  # if(p<pFilter){
  cellPvalue=cbind(Cell=colnames(rt)[i],pvalue=p)
  outTab=rbind(outTab,cellPvalue)
  # }
  mx=max(c(lowData,highData))
  lines(c(x=3*(i-1)+0.2,x=3*(i-1)+0.8),c(mx,mx))
  text(x=3*(i-1)+0.5, y=mx+0.02, labels=ifelse(p<0.001, paste0("p<0.001"), paste0("p=",sprintf("%.03f",p))), cex = 0.8)
}

legend("topright", 
       c("Low", "High"),
       lwd=4.5,bty="n",cex=1.5,
       col=bioCol[1:2])
text(seq(1,64,3),-0.05,xpd = NA,labels=colnames(rt),cex = 0.9,srt = 45,pos=2)
dev.off()
##############################################################333
pdf(paste("cor_genedown_cell.pdf",sep=''),width =7,height = 8)
pheatmap(data.r, 
         color=myColor,
         breaks=myBreaks,
         clustering_method="average", cluster_rows=F, cluster_cols=F,display_numbers=sig.mat
)

dev.off()
#################################################
pdf(paste('linkET_Lowlow.pdf',sep=''),width=length(colnames(CIBER))/2.6,height=length(colnames(CIBER))/2.6)
qcorrplot(correlate(out), type = "upper", diag = FALSE) +#热图绘制
  geom_square() +#热图绘制
  geom_couple(aes(colour = pd, size = rd),data = mantel,curvature = nice_curvature()) +#aes里面是线条格式，data对应的是mantel test 计算结果，curvature控制线条曲率
  scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(11, "RdBu")), limits = c(-0.5, 0.5)) +
  scale_size_manual(values = c(0.5, 1, 2)) +
  scale_colour_manual(values = color_pal(3)) +
  guides(size = guide_legend(title = "r",##guides()函数调整标签样式
                             override.aes = list(colour = "grey35"), 
                             order = 2),
         colour = guide_legend(title = "p", 
                               override.aes = list(size = 3), 
                               order = 1),
         fill = guide_colorbar(title = "R", order = 3))
dev.off()
######################################################3
pdf(paste0("ANN_struct.pdf"),width=40,height=40)
p <- plot(nn,col.text = "black",col.entry='black',col.intercept='#4527A0',col.out.synapse='#DCE775',col.hidden.synapse='#00695C',col.entry.synapse='#D32F2F',col.out='#DCE775',fontsize=30, col.hidden ='#00695C',radius=0.1, show.weights = FALSE,information=F)
print(p)
dev.off()