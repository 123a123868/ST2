# 检查两组之间是否有差异

# 设置
rm(list = ls()) 
options(stringsAsFactors = F)
library(stringr)
library(data.table)
library(data.table)
library(stringr)
library(stringi)
library(data.table)
library(data.table)
library(dplyr)
library(survminer) 
library(survival)
library(loose.rock)
library(futile.logger) 
library(glmSparseNet)
library(ggrisk)  
library(pheatmap) 
library(ggplot2)
library(ggsci)
library(ggstatsplot)

# 读入数据
library("Biobase")
library("affy")
library("limma")

#读入数据
data.raw <- ReadAffy(celfile.path = "00.rawData/GSE66360_RAW/GSM1620819_EA10113_142127_HG-U133_PLUS_2_SC1.CEL.gz/")
data.rma <- rma(data.raw)    # 好像已经经过Mas5或rma方法进行背景校正，标准化（归一化）
expr <- exprs(data.raw)
expr[1:4,1:4]
colnames(expr)
rownames(expr)



# GSE61145
library("GEOquery") #加载GEOquery包


library(lumi)
library(limma)
# 读取 illumina beadchip， 读取校正后的数据
##read.ilmn函数读入原始数据
data <- read.ilmn("00.RawData///GSE60993_non-normalized.txt.gz",probeid="ID_REF",other.columns="Detection Pval")
##2.1 使用neqc函数预处理原始数据
data1 <- neqc(data,detection.p="Detection Pval",offset=100) 
exp1 <- data$E

help(neqc)
# id 转换
# GPL6884 平台 直接点进去平台下载
# 读入数据
GPL6884 =  fread('00.rawData/GPL6884-11607.txt', header = TRUE)
GPL6884[GPL6884$Symbol == "IL1RL1",]


# 先转换
expr = as.data.frame(exp1)
expr = expr[GPL6884$ID,]
kp = intersect(rownames(expr),GPL6884$ID)
expr = expr[kp,]
GPL6884 = as.data.frame(GPL6884)
rownames(GPL6884) = GPL6884$ID
GPL6884 = GPL6884[kp,]
identical(rownames(expr),GPL6884$ID)
expr$Symbol = GPL6884$Symbol

# 多个探针对应一个基因的情况。
# 取平均值，进行去重
colnames()
len = ncol(expr) -1
expr[,2:len] <- apply(expr[,2:len], 2, as.numeric)  # 转换数字格式
expr$mean=apply(expr[,2:len],1,mean) 
expr=expr[order(expr$Symbol,expr$mean,decreasing = T),]
expr=expr[!duplicated(expr$Symbol),]
rownames(expr)=expr$Symbol
expr[1:4,1:4] 
colnames(expr)
expr = as.data.frame(expr[,-c(25:35)])
colnames(expr)

save(expr,file = "00.rawData/expr.GSE60993.Rdata")
# 读入分组数据
# 先看看两组之间的差异。 IL1RL1 就是ST2

# expr[1:4,1:4] #查看expr这个矩阵的1至4行和1至4列，逗号前为行，逗号后为列
# boxplot(expr[,1:4],las=2)  
# # expr=log2(expr+1)
# # boxplot(expr[,1:4],las=2)  
# # library(limma)
# expr=normalizeBetweenArrays(expr)
# boxplot(expr[,1:4],las=2)  
# expr = as.data.frame(expr)

data_gene = expr["IL1RL1",]
data_gene = as.data.frame(t(data_gene))
data_gene$group = "na"
row.names(data_gene)
data_gene$group[1:7] = "Normal"
data_gene$group[8:24] = "AMI"
table(data_gene$group)
data_gene$group = factor(data_gene$group,levels = c("Normal","AMI"))
#### 画图
library(ggplot2)
colour = c("#7CFC00","#DB423E")
colnames(data_gene)[1] = "gene"
data_gene$gene = as.numeric(data_gene$gene )



b2 = ggplot(dat = data_gene, aes(group,gene))+
  geom_boxplot(aes(fill = group))+
  scale_fill_nejm()+
  scale_fill_manual(values= colour)+
  stat_compare_means(method='wilcox.test') + labs(x='', title = 'GSE60993',y = "ST2-Expression")+
  theme_bw()+
  labs()
  theme(panel.grid = element_blank(),
        plot.title = element_text(face = 'bold', hjust = 0.5),
        axis.title = element_text(face = 'bold'),
        legend.position = '',
        text = element_text(size = 20))
b2
ggsave(b2,filename = '01.boxpolt/boxpolt.ST2.pdf',height = 7,width = 5,family = "Times")

colour = c("#7CFC00","#DB423E")
p1 <- ggplot(data_gene,aes(x=group,y=gene,fill=group))+
  geom_boxplot(width=0.6,alpha=0.8)+ stat_compare_means()+
  labs(x='', title = 'GSE60993',y = "ST2-Expression")+
  theme_bw()+
  scale_fill_manual(values= colour)+
  theme(panel.grid = element_blank(),
        plot.title = element_text(face = 'bold', hjust = 0.5),
        axis.title = element_text(face = 'bold'),
        legend.position = '',
        text = element_text(size = 20))+
  ggtitle(data_gene)
p1

ggsave(p1,filename = '01.boxpolt/boxpolt.ST2HL.pdf',height = 7,width = 5,family = "Times")








# ROC 曲线
# 设置基因 和分类 
ROC.data = data.all
colnames(ROC.data)
class.name = "prognosis"
gene = c("ST2")
expr.gene = as.numeric(data_gene$gene)
class = as.character(data_gene$group)

# ROC 分析
# 创建ROC对象
library(pROC)
library(ggsci)
library(ggplot2)
#做分类ROC曲线 可平滑
roc_obj = roc(class,expr.gene,smooth=T)
# 绘制ROC曲线,计算AUC及其置信区间
plot(roc_obj,legacy.axes = TRUE,thresholds="best", # 基于约登指数选择roc曲线最佳阈值点
     print.thres="best")
auc_value <- auc(roc_obj) # 0.8786
auc_value
auc_ci <- ci.auc(roc_obj) # 0.8575-0.8997



#将修改过后的名字替换为roc1对象中基因的名称
ggroc(roc_obj,legacy.axes = T,color="#B2533E",linewidth=0.8)+
  annotate(geom = "segment", x = 0, y = 0, xend =1, yend = 1,
           linetype="dashed",color="#186F65")+
  annotate("text",x=0.8,y=0.3,label="AUC = 0.8117")+
  labs(x="1-Specificity",y="Sensitivity")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        text = element_text(size = 20))

ggsave("01.boxpolt//GSE60993.ST2.ROC.pdf",height = 6,width = 6,family="Times")
