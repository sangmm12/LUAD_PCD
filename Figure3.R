library(survminer)
rt=risk[order(risk$riskscore),]
rt <- rt[order(rt$group), ] #另加
diff=survdiff(Surv(Time, Status) ~ group,data = rt)
pValue=1-pchisq(diff$chisq,df=1)
pValue=signif(pValue,4)
pValue=format(pValue, scientific = TRUE)
fit <- survfit(Surv(Time, Status) ~ group, data = rt)
surPlot=ggsurvplot(fit,
                   data=rt,
                   #font.title = paste(connam[i],sep=''),
                   #ggtitle = paste(connam[i],sep=''),
                   #conf.int=TRUE,
                   legend.labs=c(unique(rt$group)),
                   legend = "top",
                   legend.title="Risk",
                   pval=paste0("p=",pValue),
                   pval.size=5,
                   xlab="Time(years)",
                   # break.time.by = ceiling((max(rt$Time))/4),
                   
                   break.time.by = 1,
                   
                   risk.table.title="",
                   palette=c("red","green"),
                   risk.table=T,
                   risk.table.height=.25,)
pdf(file=paste(name_fold,'/survival.pdf',sep=''),onefile = FALSE,width = 6,height =6)
print(surPlot)
dev.off()
}


#########################################################################


col_color <- list(Class=c(tcga='#ff2e63',GSE30219='#08d9d6',GSE31210='#8785a2',GSE50081='#ffc7c7',GSE68465='#4797b1',GSE72094='#c5ecbe'))

myColor <- colorRampPalette(c("white", "#FF8A65"))(1000)
myBreaks <- seq(0.5,1,length.out=1000)


Low <- pheatmap(final_result, 
                annotation_col=col_groups,
                annotation_colors=col_color,
                color=myColor,
                gaps_row=seq(1,length(rownames(final_result))),
                gaps_col=seq(1,length(colnames(final_result))),
                breaks=myBreaks,
                border_color='black',
                fontsize_row = 15,
                fontsize_col = 15,
                fontsize_number = 10, 
                cluster_rows=F,
                cluster_cols=F,
                cellwidth = 25,
                cellheight = 15,
                display_numbers=sig.mat,
                show_colnames = F,
                number_color='black',
                fontsize=10)
A1 = as.ggplot(Low)

pdf('out_full_减4GSE.pdf',height = 213.1,width=13)









######################################################################
out_put_list <- c()
for(gene_name in names(sorted_roc_curve_list)) {
  roc_curve <- sorted_roc_curve_list[[gene_name]]$roc_curve
  pdf(paste(gene_name, '_TN.pdf', sep = ''), family = "Times", pointsize = 20)
  p <- plot(roc_curve, main = paste("ROC Curve for Gene Expression of ", gene_name, sep=''), col = "#E57373", lwd = 4, print.auc = TRUE, auc.polygon = TRUE, auc.polygon.col = "#FFCDD2")
  print(p)
  dev.off()
  
  plot(roc_curve, main = paste("", gene_name, sep=''), col = "black", lwd = 4, print.auc = TRUE, auc.polygon = TRUE, auc.polygon.col = "#FFCDD2")
  eval(parse(text = paste(gene_name, '<- recordPlot()')))
  out_put_list <- append(out_put_list, c(gene_name, NULL))
  
  auc_normal_cold <- c(auc_normal_cold, as.numeric(p$auc))
}

# 将排序后的ROC曲线拼成网格图
pdf("auc_TLSICD_HHLL_model.pdf", width = 48, height = 810, family = "Times", pointsize = 20)  
eval(parse(text = paste('print(cowplot::plot_grid(', paste(out_put_list, collapse = ","), ', ncol=8))', sep = ''))) #一行8个
dev.off()



