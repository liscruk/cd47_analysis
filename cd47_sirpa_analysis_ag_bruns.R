library(DESeq2)
library("ggplot2")

# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE178965
counts = read.table("GSE178965_count_matrix.txt.gz",header = T,row.names = 1)
des = as.data.frame(colnames(counts))
des$condition = c("GOOD","RELAPSE","RELAPSE","RELAPSE","RELAPSE","GOOD","GOOD",
                  "GOOD","GOOD","GOOD","GOOD","RELAPSE","RELAPSE","RELAPSE",
                  "RELAPSE","RELAPSE","GOOD","GOOD","GOOD")
colnames(des) = c("sample","condition")
rownames(des) = des$sample

dds = DESeq2::DESeqDataSetFromMatrix(counts,colData = des,design = ~condition)
dds = DESeq2::DESeq(dds)

res = as.data.frame(results(dds))
res = res[complete.cases(res),]

sirpa = plotCounts(dds, gene="SIRPA", intgroup="condition",returnData = T)
cd47 = plotCounts(dds, gene="CD47", intgroup="condition",returnData = T)

Q_sirpa = quantile(sirpa$count, probs=c(.25, .75), na.rm = FALSE)
Q_cd47 = quantile(cd47$count, probs=c(.25, .75), na.rm = FALSE)
iqr_s = IQR(sirpa$count)
iqr_c = IQR(cd47$count)

sirpa_mod = subset(sirpa, sirpa$count > (Q_sirpa[1] - 1.5*iqr_s) & sirpa$count < (Q_sirpa[2]+1.5*iqr_s))

ggplot(sirpa_mod, aes(x=condition, y=count))+
  geom_boxplot()+
  geom_point(position=position_jitter(w=0.05,h=0))+
  ggtitle("SIRPA")+
  theme_classic()

cd47_mod = subset(cd47, cd47$count > (Q_cd47[1] - 1.5*iqr_c) & cd47$count < (Q_cd47[2]+1.5*iqr_c))

rbind(cd47_mod,sirpa_mod)

ggplot(cd47_mod, aes(x=condition, y=count))+
  geom_boxplot()+
  geom_point(position=position_jitter(w=0.05,h=0))+
  ggtitle("CD47")+
  theme_classic()

cd47_mod$gene = "CD47"
sirpa_mod$gene = "SIRPA"

df = rbind(cd47_mod,sirpa_mod)

write.table(df,file = 'DLBCL_normalized_counts.csv',quote = F)

ggplot(df,aes(x=condition,y=count,color=condition))+
  facet_grid(~gene)+
  geom_violin(colour = "black",trim = F)+
  geom_point(colour = "black", size = 2.5)+
  geom_point()+
  scale_color_manual(values=c("red","white"))+
  scale_x_discrete(name ="",labels = c("Good Prognosis","Bad Prognosis"))+
  scale_y_continuous(breaks=c(0,250,500,750,1000,1250,1500),limits = c(0,1500),name="Normalized count")+
  theme_classic()+
  theme(text = element_text(size = 20)) 
