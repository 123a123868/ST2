# 差异分析
rm(list = ls()) 
options(stringsAsFactors = F)
library(limma)
library(edgeR)

# 读入数据
load(file = '00.rawData/expr.GSE60993.Rdata')
expr[1:4,1:4] 
colnames(expr)

# 先把normal的去掉
expr = expr[,grep('Normal',colnames(expr),invert = T)]
expr = trunc(expr)
# write.csv(expr,"00.rawData/expr.csv")
# write.csv(data_gene,"00.rawData/group.csv")
# log化
expr[1:4,1:4] #查看expr这个矩阵的1至4行和1至4列，逗号前为行，逗号后为列
boxplot(expr[,1:8],las=2)  
expr=log10(expr+1)
# boxplot(expr[,1:4],las=2)  
# library(limma)
expr=normalizeBetweenArrays(expr)
boxplot(expr[,1:4],las=2)  


# 再更具高低表达进行分组
expr = as.data.frame(expr)
data_gene = expr["IL1RL1",]
data_gene = as.data.frame(t(data_gene))
data_gene$group = "na"
median(data_gene$IL1RL1)
mean(data_gene$IL1RL1)
data_gene$group = ifelse(data_gene$IL1RL1 > mean(data_gene$IL1RL1), "ST2.H", "ST2.L")   #分组很好用
table(data_gene$group)
group = data_gene$group
save(expr,data_gene,file = '02.DEG//expr.deg.Rdata')

# 设定好分组
identical(rownames(data_gene),colnames(expr))
group_list = data_gene$group
exprSet = expr
# 指定两组进行差异分析
design <- model.matrix(~0+factor(group_list))
colnames(design)=levels(factor(group_list))
rownames(design)=colnames(exprSet)
head(design)
table(group_list)
contrast.matrix<-makeContrasts("ST2.H-ST2.L",levels = design)
contrast.matrix ##这个矩阵声明，，前比后
fit <- lmFit(exprSet,design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)  ## default no trend !!!
deg = topTable(fit2, coef=1, n=Inf,adjust='BH',sort.by="logFC")


# 设置上调下调分组
deg$group=ifelse(deg$P.Value>0.05,'stable', 
                 ifelse( deg$logFC >0.58,'up', 
                         ifelse( deg$logFC < -0.58,'down','stable')))
table(deg$group)
deg["IL1RL1",]



# 检查是否分组分错
library(ggplot2)
library(ggpubr)
deg = deg[order(deg$logFC,decreasing = T,na.last = NA),]
gene.cheak = rownames(deg)[1]
gene.cheak = "IL1RL1"
data.check = as.data.frame(exprSet[gene.cheak,])
data.check = as.data.frame(t(data.check))
identical(row.names(data_gene),rownames(data.check))
data.check$group = data_gene$group
colnames(data.check)[1] = "gene"
# 绘制箱线图
p1 <- ggplot(data.check,aes(x=group,y=gene,fill=group))+
  geom_boxplot(width=0.6,alpha=0.8)+ stat_compare_means()+
  labs(x='',y = "ST2-Expression")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        plot.title = element_text(face = 'bold', hjust = 0.5),
        axis.title = element_text(face = 'bold'),
        legend.position = '',
        text = element_text(size = 20))+
  ggtitle(gene.cheak)
p1

ggsave(p1,filename = '01.boxpolt/boxpolt.ST2HL.pdf',height = 7,width = 5,family = "Times")





ggsave("02.DEG//check.group.pdf",height = 7,width = 6,family="Times")
save(deg,file = '02.DEG//deg.Rdata')
write.csv(deg,file = '02.DEG//deg.csv')
