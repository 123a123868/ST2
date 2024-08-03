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
sce.all = readRDS ('sce.all.rds')
dim(sce.all)
table(sce.all@meta.data$orig.ident)


# 删除CTR的亚群。
levels(Idents(sce.all))   #查看细胞亚群 
sce.all = sce.all[,!sce.all@meta.data$orig.ident %in% c("MI_CTR")]
table(sce.all@meta.data$orig.ident)
#保存质控好的数据
saveRDS(sce.all, "sce.all.delt.rds")


# 简单展示
colnames(sce.all@meta.data)
VlnPlot(sce.all,group.by ='orig.ident',features = c("Il1rl1"))


# ump图
# 优化umap图
# devtools::install_github("sajuukLyu/ggunchull", type = "source")
library(ggunchull)
library(ggplot2)
library(randomcoloR)
palette <- distinctColorPalette(40)
palette<-c('#4363d8', '#3cb44b', '#ffe119', '#e6194b', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000')
data <- as.data.frame(sce.all[["umap"]]@cell.embeddings)
data$cell_type <- sce.all$cell_type
colnames(data)
table(data$cell_type)
colnames(data)
ggplot(data, aes(x = umap_1, 
                 y = umap_2, 
                 fill = cell_type,
                 color = cell_type)) +
  stat_unchull(data=subset(data, cell_type=="Fibroblasts"),
               alpha = 0.2, 
               size = 0.5,
               show.legend = F,
               nbin = 300, 
               nsm = 20,
               qval = 0.8,
               sfac = 1.5) +
  geom_point(size = 1) +
  theme_classic()+
  theme(axis.text = element_text(colour = 'black', size = 12))+
  scale_fill_manual(values = palette) +
  scale_color_manual(values = palette)
ggsave('01.polt/ump.cluster.pdf',width = 8,height = 8,family="Times")


# 展示基因的表达。
FeaturePlot(sce.all, features = c("Il1rl1"),dims = c(1, 2),blend.threshold = 4)
ggsave('01.polt/ump.Il1rl1.pdf',width = 5,height = 5,family="Times")

VlnPlot(sce.all, features = "Il1rl1")
colnames(sce.all@meta.data)
VlnPlot(sce.all, features = "Il1rl1",pt.size=0,group.by = "cell_type")+
  #width控制箱体宽度，col控制边框颜色，fill控制填充颜色
  geom_boxplot(width=.2,col="black",fill="white")+ 
  NoLegend()


# 自己画umap展示基因表达

# 提取UMAP坐标
umap_df <- as.data.frame(sce.all@reductions$umap@cell.embeddings)
umap_df$cluster <- as.factor(sce.all@meta.data$cell_type)
head(umap_df)
# 提取基因表达数据并与UMAP坐标合并
gene_df <- as.data.frame(GetAssayData(object = sce.all, slot = "data")[c("Il1rl1"), ])
identical(rownames(umap_df),rownames(gene_df))
merged_df = umap_df
merged_df$Il1rl1 = gene_df$`GetAssayData(object = sce.all, slot = "data")[c("Il1rl1"), ]`
head(merged_df)


#绘图 
library(ggnewscale)
colnames(merged_df)
ggplot(merged_df, vars = c("umap_1", "umap_2", "Il1rl1"), aes(x = umap_1, y = umap_2, colour = Il1rl1)) +
  geom_point(size=0.3, alpha=1) +
  scale_colour_gradientn(colours = c("lightgrey", "red"), limits = c(0, 0.3), oob = scales::squish) +
  new_scale_color()+theme_classic()
ggsave('01.polt/ump2.Il1rl1.pdf',width = 5,height = 5,family="Times")



# V5获取细胞矩阵
# 指定要提取的基因名称
gene_name <- "Il1rl1"
# df1 = as.data.frame(sce.all[["RNA"]]$counts.MI_CTR.1)[gene_name,]
df2 = as.data.frame(sce.all[["RNA"]]$counts.MI_1D.2)[gene_name,]
df3 = as.data.frame(sce.all[["RNA"]]$counts.MI_3D.3)[gene_name,]
df4 = as.data.frame(sce.all[["RNA"]]$counts.MI_5D.4)[gene_name,]
df5 = as.data.frame(sce.all[["RNA"]]$counts.MI_7D.5)[gene_name,]

# 合并多个矩阵
library(dplyr)
combined_df <- bind_cols( df2, df3,df4,df5)
combined_df = as.data.frame(t(combined_df))
# 读入分组数据
meta = sce.all@meta.data
identical(row.names(combined_df),row.names(meta))
meta$geneExpr = combined_df$Il1rl1
table(meta$cell_type)

# ggpolt画图
# Multiple groups with error bars and jitter point
library(ggpubr)
library(ggplot2)
colnames(data.long)
unique(meta$cell_type)
data.long = meta[meta$cell_type == "Macrophages",]
table(data.long$cell_type)
data.long$group = data.long$orig.ident
table(data.long$group)
max(data.long$geneExpr)
ylim = 0.2

unique(data.long$orig.ident)
mean(data.long[data.long$orig.ident == "MI_1D" ,]$geneExpr)
mean(data.long[data.long$orig.ident == "MI_3D" ,]$geneExpr)
mean(data.long[data.long$orig.ident == "MI_5D" ,]$geneExpr)
mean(data.long[data.long$orig.ident == "MI_7D" ,]$geneExpr)


ggbarplot(data.long, x = 'group', y = 'geneExpr', add = c('mean_se'),
          color = "group", fill = "group",palette =  mycolour,
          position = position_dodge(1),width = 0.5) +
  labs(y = "geneExpr") + 
  theme( legend.position = 'none' ) +
  rotate_x_text(angle = 0)+
  coord_cartesian(ylim = c(0, ylim))+  # 设置y轴的坐标范围 
  stat_compare_means(method = 'wilcox.test', aes(group = group),label.y = ylim)
ggsave("01.polt//geom_bar.Il1rl1.pdf",height = 5,width = 5,family="Times")



mycolour<-c('#4363d8', '#3cb44b', '#ffe119', '#e6194b', '#f58231',
            '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', 
            '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', 
            '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', 
            '#ffffff', '#000000')