# 免疫评分分析
# estimate评分

# 设置
rm(list=ls())
options(stringsAsFactors = F)
library(stringr) 
library(limma)
library(estimate)
library(tidyverse)
library(ggplot2)
library(ggstatsplot)
library(ggsci)
library(ggpubr)
library(patchwork)
library(data.table)


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
# expr=log10(expr+1)
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



# 
# # 标准化
# # 需要log
# boxplot(expr[,1:4],las=2)  
# expr=log2(expr+1)
# boxplot(expr[,1:4],las=2)
# library(limma)
# expr=normalizeBetweenArrays(expr)
# boxplot(expr[,1:4],las=2)
# colnames(expr)
# boxplot(expr[1:10,]) 

# estimate评分
#  do estimate
estimate_RNAseq <- function(RNAseq_logCPM,pro){
  input.f=paste0("05.immune/estimate_out//",pro,'_estimate_input.txt')
  output.f=paste0("05.immune/estimate_out//",pro,'_estimate_gene.gct')
  output.ds=paste0("05.immune/estimate_out//",pro,'_estimate_score.gct')
  write.table(RNAseq_logCPM,file = input.f,sep = '\t',quote = F)
  
  library(estimate)
  filterCommonGenes(input.f=input.f,
                    output.f=output.f ,
                    id="GeneSymbol")
  estimateScore(input.ds = output.f,
                output.ds=output.ds,
                platform="illumina")
  scores=read.table(output.ds,skip = 2,header = T)
  rownames(scores)=scores[,1]
  scores=t(scores[,3:ncol(scores)])
  scores=data.frame(  scores)
  return(scores)
}
pro <- 'ST2'
dat <- expr
scores <- estimate_RNAseq(dat, pro)
head(scores)
save(scores, file = '05.immune/estimate_out//ST2.estimate_results.Rdata')


# 画图
# 免疫评分分析
# estimate评分

# 设置
# rm(list=ls())
options(stringsAsFactors = F)
load('05.immune/estimate_out/ST2.estimate_results.Rdata')
# 肿瘤纯度计算
scores$purity = cos(0.6049872018+0.0001467884 * scores$ESTIMATEScore)
# TIDE计算：
# 就不算了

# # TIS评分
# #  TIS-基因
gene_TIS = c(
  "TIGIT",
  "CD27",
  "CD8A",
  "PDCD1LG2", #(PD-L2)
  "LAG3",
  "CD274", #(PD-L1)
  "CXCR6",
  "CMKLR1",
  "NKG7",
  "CCL5",
  "PSMB10",
  "IDO1",
  "CXCL9",
  "HLA-DQA1",
  "CD276",
  "STAT1",
  "HLA-DRB1",
  "HLA-E")
# 匹配
expr.gbm = expr
table(gene_TIS %in% rownames(expr.gbm))
gene_TIS = gene_TIS[gene_TIS %in% rownames(expr.gbm)]
expr_TIS = expr.gbm[gene_TIS,]
# 计算平均值
data_TIS = data.frame(sample = colnames(expr_TIS),
                      TIS_mean = apply(expr_TIS, 2,mean)
) 
# 匹配一下数据
identical(rownames(scores),rownames(data_TIS))
data_all = cbind(scores,data_TIS)

# 分组，只计算PEI和TN组
rownames(data_all)
identical(rownames(data_all),rownames(data_gene))
data_all$group = data_gene$group
data_all$group=factor(data_all$group, levels=c("ST2.L", "ST2.H"))

# 标准化便于展示
data_all$StromalScore
# 使用 apply() 函数对每一列数据进行标准化
data_all.st = data_all
data_all.st[,c(1:4,6)] <- apply(data_all.st[,c(1:4,6)], 2, scale)

data_all.st = data_all.st[,c(1,2,7)]



# 转换长数据
colnames(data_all)
data.long <- melt(data_all.st, id = c("sample","group"))
colnames(data.long) = c("sample" ,  "group"  ,  "score.class" ,"score.value" )
table(data.long$group)
# 画图
mycol = c("#008ECA","#DB423E")   #分组的，通用这两个配色吧
ggplot(data.long,aes(x=group, y=score.value)) + 
  geom_boxplot(aes(color=group),alpha=0,width=0.6)+
  scale_color_manual(values = mycol)+
  facet_grid(.~score.class)+
  geom_jitter(aes(color=group),width = 0.2)+
  theme_bw()+
  geom_signif(comparisons = list(c('ST2.L','ST2.H')),map_signif_level = T,
              test = 't.test')+ #或者换成t检验，就是t.test  #wilcox.test
  ylab('score(log)')+
  guides(color=FALSE)+
  theme(title = element_text(size = 10,colour = 'black'),##设置主题
        plot.title = element_text(hjust = 0.5),
        axis.text.y = element_text(size = 11,colour = 'black'),
        axis.text.x = element_text(size = 11,colour = 'black'),
        axis.title.x = element_text(size = 10,colour = 'black'),
        axis.title.y = element_text(size = 10,colour = 'black'),
        strip.text.x = element_text(size = 10,colour = 'black'))
ggsave("05.immune/01.estimate.pdf",width = 8,height = 5)

