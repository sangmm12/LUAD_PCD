#绘制boxplot
p=ggboxplot(data, x="Gene", y="Expression", color = "group", 
            ylab="ssGSVA-score",
            xlab="",
            legend.title=x,
            palette = c("#FF4940","#028E9B"),
            width=0.6, add = "none")
p=p+rotate_x_text(60)
p1=p+stat_compare_means(aes(group=group),
                        method="wilcox.test",
                        symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                        label = "p.signif")

#输出图片
pdf(file=outFile, width=6, height=12)
print(p1)
####################################################3
#绘制boxplot
p=ggboxplot(data, x="Gene", y="Expression", color = "group", 
            ylab="ssGSVA-score",
            xlab="",
            legend.title=x,
            palette = c("#FF4940","#028E9B"),
            width=0.6, add = "none")
p=p+rotate_x_text(60)
p1=p+stat_compare_means(aes(group=group),
                        method="wilcox.test",
                        symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                        label = "p.signif")

#输出图片
# pdf(file=outFile, width=6, height=12) #原
pdf(file=outFile, width=8, height=8)
##########################################################
pdf(file=paste(cancer,"_","TIDE.pdf", sep = ""), width=3.5, height=4.5)

ggviolin(data, x="risk", y="TIDE", fill = "risk", 
         xlab="", ylab="TIDE",
         palette=c("#FF4940","#028E9B"),  # 使用红色和绿色"#E64B35B2","#00A087B2"
         legend.title="Risk",
         add = "jitter", add.params = list(fill="white")) + 
  stat_compare_means(comparisons = my_comparisons)

dev.off()
######################################################3
p <- ggplot(var_out, aes(x = reorder(group, -vimp), y = vimp)) + 
  geom_bar(stat = "identity",   
           show.legend = FALSE,   
           width = .7,fill='#D0006E') + aes(fill=vimp)+
  xlab("Gene") + 
  ylab("Vimp")+  theme_classic()+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
print(p)
dev.off()
#################################################3
pdf(paste0("ANN_struct.pdf"),width=40,height=40)
p <- plot(nn,col.text = "black",col.entry='black',col.intercept='#4527A0',col.out.synapse='#DCE775',col.hidden.synapse='#00695C',col.entry.synapse='#D32F2F',col.out='#DCE775',fontsize=30, col.hidden ='#00695C',radius=0.1, show.weights = FALSE,information=F)
print(p)
dev.off()
##################################################3
pdf("ann_rank.pdf",width=15)
p <- ggplot(weight, aes(x = reorder(rownames(weight), -sum), y = sum)) +
  geom_bar(stat = "identity",   
           show.legend = FALSE,   
           width = .8,fill = '#01579B') +
  labs(x = "FeatureName", y = "weight") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
p
#print(p)
dev.off()
################################################################3
pdf("XGboost_rank.pdf",width=13)
p <- ggplot(feature_importance, aes(x = reorder(feature_importance$Feature, -feature_importance$Gain), y = feature_importance$Gain)) +
  geom_bar(stat = "identity",   
           show.legend = FALSE,   
           width = .8,fill = '#F57C00') +
  labs(x = "FeatureName", y = "Gain") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
p
#print(p)
dev.off()
#################################################################3
pdf("svm_rank.pdf",width=13)
p <- ggplot(top.features, aes(x = reorder(FeatureName, AvgRank), y = AvgRank)) +
  geom_bar(stat = "identity",   
           show.legend = FALSE,   
           width = .8,fill = '#D1C4E9') +
  labs(x = "FeatureName", y = "rank") + theme_classic()+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
p
#print(p)
dev.off()
############################################################
pdf('TIDE-ImportanceRank.pdf',width=widelength/3+5,height=10)
ggplot(my_data, aes(x = as.factor(Sample), y = Value, fill = Category)) +
  geom_bar(stat = "identity") +
  facet_grid(Category ~ ., scales = "free_y", space = "free") +
  scale_fill_manual(values = color) +
  labs(title = "", x = "PCD", y = "Important") +
  scale_x_discrete(labels = my_data$name) +  
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45,size=14,color='black', hjust = 1, family = my_font),
        axis.text.y = element_text(size=14,color='black', family = my_font))
dev.off()
############################################################
library(ggpubr)
pdf(paste("cluster_TIP .pdf",sep=''),width=13,height = 6)
p <- ggboxplot(dat, x = "Gene", y = "value",
               color = "Group", palette = c("#A60000","#009999"), 
               add = "none",x.text.angle=60) # palette可以按照期刊选择相应的配色，???"npg"???
p <- p+xlab("Step")+ylab("Expression of TIP")
p <- p + stat_compare_means(aes(group = Group), label = "p.signif",method = "anova")
p
#print(p)
dev.off()
#######################################################
pdf("Stimulaotry_vioplot_adjusted.pdf", width = 13, height = 8)
par(las = 1.5, mar = c(15, 6, 3, 3))  # ?±߾???????15?????ɱ?ǩ
x <- c(1:ncol(rt_log2))
y <- c(1:ncol(rt_log2))
plot(x, y,
     xlim = c(1, ncol(rt_log2) * 3), ylim = c(min(rt_log2), max(rt_log2) + 0.02),  #x,y????ǩ
     main = "", xlab = "", ylab = "Fraction",
     pch = 21,
     col = "white",
     xaxt = "n")