library(timeROC)
predict_full <- risk$riskscore
full_time <- mixed$Time/365
ROC_rt=timeROC(T=full_time,delta=mixed$Status,
               marker=predict_full,cause=1,
               weighting='aalen',
               times=c(1,2,3),ROC=TRUE)
pdf(file=paste(name_fold,'/roc.pdf',sep=''),width = 6,height = 6)
plot(ROC_rt,time=1,col='green3',title=FALSE,lwd=2)
plot(ROC_rt,time=2,col='blue4',add=TRUE,title=FALSE,lwd=2)
plot(ROC_rt,time=3,col='darkred',add=TRUE,title=FALSE,lwd=2)
legend('bottomright',
       c(paste0('AUC at 1 years: ',sprintf("%.03f",ROC_rt$AUC[1])),
         paste0('AUC at 2 years: ',sprintf("%.03f",ROC_rt$AUC[2])),
         paste0('AUC at 3 years: ',sprintf("%.03f",ROC_rt$AUC[3]))),
       col=c("green3",'blue4','darkred'),lwd=2,bty = 'n')
dev.off()
##########################################################33
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
                   break.time.by = ceiling((max(rt$Time))/4),
                   risk.table.title="",
                   palette=c("red","green"),
                   risk.table=T,
                   risk.table.height=.25,)
pdf(file=paste(name_fold,'/survival.pdf',sep=''),onefile = FALSE,width = 6,height =6)
print(surPlot)
dev.off()
}
##################################################################3
pdf("20PCDscore_surv_algorithm_rank.pdf",width=widelength/3+2,height=10)
p <- ggplot(my_data, aes(x = as.factor(Sample), y = Value, fill = Category)) +
  geom_bar(stat = "identity") +
  facet_grid(Category ~ ., scales = "free_y", space = "free") +
  scale_fill_manual(values = color) +
  labs(title = "", x = "Gene", y = "Important") +
  scale_x_discrete(labels = my_data$name) +  
  theme_classic()+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 16),  # 调整x轴标签文本大小
    axis.text.y = element_text(size = 16),  # 调整y轴标签文本大小
    axis.title = element_text(size = 18),  # 调整轴标题文本大小
    plot.title = element_text(size = 20),  # 调整图标题文本大小
    strip.text = element_text(size = 14)  # 调整facet标签文本大小
  )
print(p)
dev.off()
