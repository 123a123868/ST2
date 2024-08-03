# 设置
rm(list=ls())
options(stringsAsFactors = F) 
# 设置
#准备工作
library(Seurat)
library(SingleR)
library(SeuratWrappers)
library(monocle3)
library(reticulate)
library(Matrix)
library(ggplot2)
library(cowplot)
library(magrittr)
library(patchwork)
library(dplyr)
library(scater)

# 读入数据
sce.all = readRDS ('sce.all.delt.rds')
dim(sce.all)
table(sce.all@meta.data$orig.ident)


# 下载信号通路
library(msigdbr)
library(RColorBrewer)
# msigdbr_show_species() 查看物种
# Homo sapiens
# Mus musculus
gene_sets <- msigdbr(species = "Mus musculus", category = "C2") %>% 
  dplyr::filter(gs_name=="BIOCARTA_IL1R_PATHWAY")


Idents(sce.all) <- sce.all$cell_type
TCR_l <-list(TCR_signaling=unique(gene_sets$gene_symbol)) 
signatures <- TCR_l


# 读入数据
sce.all = readRDS ('sce.all.rds')
# 主函数
library(UCell)
# # 安装包
# install.packages('01.polt/UCell/',repos=NULL, type="source")  # 文件夹
# library(UCell)

sce.all@assays$RNA$counts = sce.all@assays$RNA$scale.data
sce.all@assays$RNA$counts = sce.all@assays$RNA$counts.MI_1D.2

rds <- AddModuleScore_UCell(sce.all, features = signatures)
signature.names <- paste0(names(signatures), "_UCell")

# 小提琴图展示每一个细胞
options(repr.plot.width=6, repr.plot.height=4)
VlnPlot(rds, features = signature.names, 
        group.by = "celltype",
        cols=cuscolors )
library(data.table)
library(ggplot2)
library(tidyverse)
library(ggplot2)
library(ggstatsplot)
library(ggsci)
library(ggpubr)
library(patchwork)



# 标注p值
#画箱形图
require(tidyr)
require(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(tibble)
colour = c("#7CFC00","#DB423E") 

# 宽变长，变成tidy data frame，好用于ggplot2画图
#转换长数据
library(stringi)
data_new = data.frame(
  group = rds$sample.id,
  gene = rds$celltype,
  expression = rds$kp_UCell
)
data_new$group2 = stri_sub(data_new$group,1,1)
table(data_new$group2)
table(data_new$gene)
# data_new = data_new[data_new$gene == "Fibroblasts",]
#############################################################
# 用小提琴图画
# TMB
# 用小提琴图画
data = data_new
colnames(data)
# 设置颜色
jco <- c("#DB423E","#008ECA") 
ggplot(data = data,aes(x = group2, y = expression, fill = group2))+
  scale_fill_manual(values = jco[2:1]) + 
  geom_violin(alpha=0.4, position = position_dodge(width = .75),
              size=0.8, color="black") + # 边框线黑色
  geom_boxplot(notch = FALSE, outlier.size = -1, 
               color="black", lwd=0.8, alpha = 0.7)+ # 背景色透明化
  geom_point(shape = 21, size=1, 
             position = position_jitterdodge(), 
             color="black", alpha=1)+ # 边框线黑色
  theme_classic() +
  ylab(expression("FMOne_mutation_burden_per_MB")) +
  xlab("group")  +
  stat_compare_means(aes(label = ..p.signif..),
                     comparisons = list(c('C',"L")
                     ),
                     method="wilcox.test"
  )+#添加检验
  theme(#panel.border = element_rect(colour = "black", fill=NA, size=0.2), # 原图有框
    axis.ticks = element_line(size=0.2,color="black"),
    axis.ticks.length = unit(0.2,"cm"),
    legend.position = "none",
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10))
# 保存图像
ggsave("9-scRNA_score//score-2.pdf", width = 3, height = 4)



# 相关性
data_cnv = read.csv("10-inferCNV//inferCNV-figures/cnv_scores.csv")
rownames(data_cnv) = data_cnv$X
kp = intersect(rownames(data_cnv),rownames(data_new))
data_cnv = data_cnv[kp,]
data_new = data_new[kp,]
identical(rownames(data_cnv),rownames(data_new))

data_new$cnv_score = data_cnv$cnv_score
colnames(data_new)

boxplot( data_new$cnv_score ~ data_new$group2)  

# 先计算相关性
data_new$cnv_score = log10(data_new$cnv_score)
cor.test(data_new$expression, data_new$cnv_score, method="spearman")

# 简单画图
p1 <- ggplot(data_new, aes(x = expression, y = cnv_score))
p2 <- p1 + geom_point() 
p3 <- p2 + geom_smooth(method="lm")
p3
library(ggpubr)

ggscatter(data_new, x = "expression", y = "cnv_score", 
          color = "red3",fill = "lightgray",
          add = "reg.line", conf.int = TRUE, 
          add.params = list( color = "black",fill = "lightgray",fill = "lightgray"),
          cor.coef = T,
          cor.coeff.args = list(),
          cor.method = "pearson",
          cor.coef.coord = c(10, 15),
          cor.coef.size = 1)


ggsave(p1,"6-scRNA_score/cor.pdf")
ggsave(filename = )
cor.test(data_new$expression, data_new$cnv_score, method="spearman")



# umap图展示
options(repr.plot.width=12, repr.plot.height=4)
FeaturePlot(rds, features = c("kp_UCell"), 
            ncol = 3, order = T,
            min.cutoff = "q03", max.cutoff = "q99",
            cols = c("#F5F5F5", "#333399"), pt.size = 0.1)


u.scores <- ScoreSignatures_UCell(sce.all, features = signatures)
u.scores[1:8, 1:2]