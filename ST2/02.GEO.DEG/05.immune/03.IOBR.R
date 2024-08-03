library(IOBR) 
library(tidyverse)
library(IOBR)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(ggsci)
library(pheatmap)
library(cowplot)
library(ComplexHeatmap)

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
# expr = trunc(expr)
# write.csv(expr,"00.rawData/expr.csv")
# write.csv(data_gene,"00.rawData/group.csv")
# log化
expr[1:4,1:4] #查看expr这个矩阵的1至4行和1至4列，逗号前为行，逗号后为列
boxplot(expr[,1:8],las=2)  
expr=log2(expr+1)
# boxplot(expr[,1:4],las=2)  
# library(limma)
expr=normalizeBetweenArrays(expr)
boxplot(expr[,1:4],las=2) 
expr = as.data.frame(expr)
data = expr
data[1:3, 1:3]

# 单一免疫浸润算法的实现与绘图：以CIBERSORT为例
cibersort <- deconvo_tme(eset = data, 
                         method = "cibersort", 
                         arrays = FALSE, 
                         perm = 10) # 理论上来说置换次数越大则结果越准确，同样运行速度也越慢，这里以10次为例
head(cibersort)


# 单一免疫浸润算法的实现与绘图：以CIBERSORT为例
cibersort <- deconvo_tme(eset = data, 
                         method = "cibersort", 
                         arrays = FALSE, 
                         perm = 10) # 理论上来说置换次数越大则结果越准确，同样运行速度也越慢，这里以10次为例
head(cibersort)
write_tsv(cibersort, "cibersort.txt")

# # 可视化
# res <- cell_bar_plot(input = cibersort, 
#                      title = "CIBERSORT Cell Fraction",
#                      legend.position = "top",
#                      coord_filp = TRUE)
# ggsave("cibersort.pdf", width = 10, height = 8)
xcell <- deconvo_tme(eset = data, method = "xcell", arrays = FALSE)
cibersort = xcell
# 读入分组信息
load('02.DEG//expr.deg.Rdata')
identical(rownames(data_gene),cibersort$ID)
cibersort$group = data_gene$group

# ggplot2绘图
head(cibersort)
box_data <- cibersort
# 数据转换
box_data1 <- melt(box_data,
                  id.vars = c("ID", "group"))#转换成长数据
head(box_data1)

ggplot(box_data1, 
       aes(variable, value, fill = group)) +
  geom_boxplot(alpha = 0.85) +
  scale_y_continuous(name = "CIBERSORT Cell Fraction") +
  scale_fill_manual(values=c("red", "blue"))+
  scale_x_discrete(name = "") +
  theme_classic() +
  theme(text = element_text(size = 10),
        legend.position = "top",
        axis.text.x  = element_text(angle=60, vjust=1,hjust=1)) + 
  stat_compare_means(aes(label = ..p.signif..))

ggsave("05.immune/xcell3.pdf", width =20, height = 7)






# 所有免疫浸润分析，非肿瘤推荐xcell：
cibersort <- deconvo_tme(eset = data, method = "cibersort", arrays = FALSE, perm = 10)
epic <- deconvo_tme(eset = data, method = "epic", arrays = FALSE)
mcp <- deconvo_tme(eset = data, method = "mcpcounter")
xcell <- deconvo_tme(eset = data, method = "xcell", arrays = FALSE)
estimate <- deconvo_tme(eset = data, method = "estimate")
timer <- deconvo_tme(eset = data, method = "timer", group_list = rep("uvm", dim(data)[2]))
quantiseq <- deconvo_tme(eset = data, method = "quantiseq", tumor = TRUE, arrays = FALSE, scale_mrna = TRUE)

# 合并所有分析结果
tme_combine <- cibersort %>% 
  inner_join(mcp, "ID") %>% 
  inner_join(xcell, "ID") %>%
  inner_join(epic, "ID") %>% 
  inner_join(estimate, "ID") %>% 
  inner_join(timer, "ID") %>% 
  inner_join(quantiseq, "ID")
head(tme_combine)
write_tsv(tme_combine, "tme_combine.txt")






