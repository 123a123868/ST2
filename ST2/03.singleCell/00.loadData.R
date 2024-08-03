# 读入数据

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
options(future.globals.maxSize = 300000 * 1024^2)
#devtools::install_github("hhoeflin/hdf5r")

packageVersion("Seurat")


#准备输入数据
Control.data <- Read10X_h5("00.rawData/GSM4972357_Steady-state_filtered_feature_bc_matrix.h5")
MI_1D.data <- Read10X_h5("00.rawData/GSM4972358_Day1_post-MI_filtered_feature_bc_matrix.h5")
MI_3D.data <- Read10X_h5("00.rawData/GSM4972359_Day3_post-MI_filtered_feature_bc_matrix.h5")
MI_5D.data <- Read10X_h5("00.rawData/GSM4972360_Day5_post-MI_filtered_feature_bc_matrix.h5")
MI_7D.data <- Read10X_h5("00.rawData/GSM4972361_Day7_post-MI_filtered_feature_bc_matrix.h5")

#上述代码也可以写成循环的形式
#fs=list.files(path = "0.单细胞数据集处理/4.心肌梗死_空转_单细胞/GSE163129_scRNA/",pattern = '.h5')
#fs
#lapply(fs, function(x){
#  x=fs[1]
#  print(x)
#  a=Read10X_h5(x)
#  a[1:4,1:4] 
#  library(stringr)
#  (p=str_split(x,'_',simplify = T)[,1])
#  sce <- CreateSeuratObject( a ,project = p )
#  sce
#})




#创建Seurat对象
MI.control <- CreateSeuratObject(Control.data, project = "MI_CTR", min.cells = 10)
MI.1D <- CreateSeuratObject(MI_1D.data, project = "MI_1D", min.cells = 10)
MI.3D <- CreateSeuratObject(MI_3D.data, project = "MI_3D", min.cells = 10)
MI.5D <- CreateSeuratObject(MI_5D.data, project = "MI_5D", min.cells = 10)
MI.7D <- CreateSeuratObject(MI_7D.data, project = "MI_7D", min.cells = 10)

#计算QC指标
MI.control[['percent.mito']] <- PercentageFeatureSet(MI.control, pattern = "^mt-")
MI.1D[['percent.mito']] <- PercentageFeatureSet(MI.1D, pattern = "^mt-")
MI.3D[['percent.mito']] <- PercentageFeatureSet(MI.3D, pattern = "^mt-")
MI.5D[['percent.mito']] <- PercentageFeatureSet(MI.5D, pattern = "^mt-")
MI.7D[['percent.mito']] <- PercentageFeatureSet(MI.7D, pattern = "^mt-")

#计算基因的表达区间
quantile(MI.control$nFeature_RNA, probs = c(0.02,0.98))
quantile(MI.1D$nFeature_RNA, probs = c(0.02,0.98))
quantile(MI.3D$nFeature_RNA, probs = c(0.02,0.98))
quantile(MI.5D$nFeature_RNA, probs = c(0.02,0.98))
quantile(MI.7D$nFeature_RNA, probs = c(0.02,0.98))

#QC
MI.control <- subset(x = MI.control, subset = nFeature_RNA > 371 & nFeature_RNA < 5061 & percent.mito < 10)
MI.1D <- subset(x = MI.1D, subset = nFeature_RNA > 350 & nFeature_RNA < 5199 & percent.mito < 10)
MI.3D <- subset(x = MI.3D, subset = nFeature_RNA > 407 & nFeature_RNA < 7825 & percent.mito < 10)
MI.5D <- subset(x = MI.5D, subset = nFeature_RNA > 493 & nFeature_RNA < 6741 & percent.mito < 10)
MI.7D <- subset(x = MI.7D, subset = nFeature_RNA > 410 & nFeature_RNA < 6132 & percent.mito < 10)


#合并数据
experiment.aggregate <- merge(x = MI.control, y = c(MI.1D,MI.3D,MI.5D,MI.7D))
experiment.aggregate@meta.data


#添加样本来源信息在合并后的Seurat对象中
samplename = experiment.aggregate@meta.data$orig.ident
batchorder = rep("MI",length(samplename))
batchorder[samplename %in% c("MI_CTR")] = "1"
batchorder[samplename %in% c("MI_1D")] = "2"
batchorder[samplename %in% c("MI_3D")] = "3"
batchorder[samplename %in% c("MI_5D")] = "4"
batchorder[samplename %in% c("MI_7D")] = "5"
names(batchorder) = colnames(experiment.aggregate)
table(batchorder)

batchid = rep("MI",length(samplename))
batchid[samplename %in% c("MI_CTR")] = "MI_CTR"
batchid[samplename %in% c("MI_1D")] = "MI_1D"
batchid[samplename %in% c("MI_3D")] = "MI_3D"
batchid[samplename %in% c("MI_5D")] = "MI_5D"
batchid[samplename %in% c("MI_7D")] = "MI_7D"
names(batchid) = colnames(experiment.aggregate)

experiment.aggregate <- AddMetaData(object = experiment.aggregate, metadata = batchorder, col.name = "batchorder")
table(experiment.aggregate@meta.data$batchorder)

experiment.aggregate <- AddMetaData(object = experiment.aggregate, metadata = batchid, col.name = "batchid")
table(experiment.aggregate@meta.data$batchid)



#更改数据格式，方便进行CCA
MI.list <- SplitObject(experiment.aggregate, split.by = "batchid")
#标准化每个数据集并进行高可变基因的筛选
for (i in 1:length(MI.list)) {
  MI.list[[i]] <- NormalizeData(MI.list[[i]], verbose = FALSE)
  MI.list[[i]] <- FindVariableFeatures(MI.list[[i]], selection.method = "vst", 
                                       nfeatures = 2000, verbose = FALSE)
}

# 在这里作者首先对每个子数据集进行标准化以及高可变基因的筛选，然后执行CCA算法
MI.anchors <- FindIntegrationAnchors(object.list = MI.list, dims = 1:30)
MI.integrated <- IntegrateData(anchorset = MI.anchors, dims = 1:30)
MI.integrated = 
# 执行完CCA算法后，作者就开始对整合后的数据进行后续的单细胞标准流程操作
all.genes <- rownames(MI.integrated)
MI.integrated <- ScaleData(MI.integrated, features = all.genes, verbose = FALSE)
MI.integrated <- RunPCA(MI.integrated, verbose = FALSE)
ElbowPlot(object = MI.integrated, ndims = 50)
use.pcs = 1:40
MI.integrated <- RunUMAP(MI.integrated, reduction = "pca", dims = use.pcs)
MI.integrated <- FindNeighbors(object = MI.integrated, reduction = "pca", dims = use.pcs)
MI.integrated <- FindClusters(object = MI.integrated, resolution = seq(0.5,2,0.1))



#可视化
DimPlot(object = MI.integrated, reduction = 'umap', group.by = "integrated_snn_res.1.4", pt.size=0.5, label = TRUE, repel = TRUE) + NoLegend()
DimPlot(object = MI.integrated, reduction = 'umap', group.by = "batchid", pt.size=0.5)
DimPlot(object = MI.integrated, reduction = 'umap', group.by = "orig.ident", split.by = "batchid",  pt.size=0.5, ncol = 3) + NoLegend()




#然后我们选择合适的resolution参数
#并设置后续默认的细胞亚群
MI.integrated <- SetIdent(MI.integrated, value = 'integrated_snn_res.1.4')

#后续我们需要进行差异分析
#但是需要注意我们差异分析的时候就不能使用CCA后合并的数据
#而是需要使用原始的counts数据，因为我们使用CCA算法后抹去了差异
#我们可以简单理解CCA后的目的是为了展示结果
DefaultAssay(MI.integrated) <- "RNA"
MI.integrated <- NormalizeData(MI.integrated, verbose = FALSE)
all.genes <- rownames(MI.integrated)
MI.integrated <- ScaleData(MI.integrated, features = all.genes, verbose = FALSE)
markers_all_RNA <- FindAllMarkers(object = MI.integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = "LR")



#SingleR注释
mrsd.se <- ImmGenData()#调用参考数据集
MI.integrated.se <- as.SingleCellExperiment(MI.integrated, assay = "RNA")#转换数据格式
mrsd.common <- intersect(rownames(MI.integrated.se), rownames(mrsd.se))#选择交集基因
mrsd.se <- mrsd.se[mrsd.common,]
MI.integrated.se <- MI.integrated.se[mrsd.common,]
MI.integrated.se <- logNormCounts(MI.integrated.se)#归一化数据
MI.mrsd.pred <- SingleR(test = MI.integrated.se, ref = mrsd.se, method = "single", labels = mrsd.se$label.main)#注释（大类）
MI.integrated[["mrsd.main"]] <- MI.mrsd.pred$labels
MI.mrsd.pred <- SingleR(test = MI.integrated.se, ref = mrsd.se, method = "single", labels = mrsd.se$label.fine)#注释（小类）
MI.integrated[["mrsd.fine"]] <- MI.mrsd.pred$labels

########################################################################################################################
# Filter Endothelial/Fibroblast, double cell type clusteres and other technical noise[PMID:30471926] (33,35) 
MI.integrated.1st.filtered = subset(MI.integrated, idents = c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,34))

DefaultAssay(MI.integrated.1st.filtered) <- "integrated"

MI.integrated.1st.filtered <- RunPCA(MI.integrated.1st.filtered, verbose = FALSE)

ElbowPlot(object = MI.integrated.1st.filtered, ndims = 50)
use.pcs = 1:40

MI.integrated.1st.filtered <- RunUMAP(MI.integrated.1st.filtered, reduction = "pca", dims = use.pcs)
MI.integrated.1st.filtered <- FindNeighbors(object = MI.integrated.1st.filtered, reduction = "pca", dims = use.pcs)
MI.integrated.1st.filtered <- FindClusters(object = MI.integrated.1st.filtered, resolution = seq(0.5,2,0.1))
MI.integrated.1st.filtered <- SetIdent(MI.integrated.1st.filtered, value = 'integrated_snn_res.1.5')
#因为上述提取了我们感兴趣的细胞亚群
#所以我们需要再次对细胞亚群进行降维和聚类，也就是说再次进行标准流程的操作

#细胞注释
cluster2celltype <- c("0"="Macrophages",
                      "1"="Macrophages", 
                      "2"="Macrophages", 
                      "3"="Neutrophils", 
                      "4"="Macrophages", 
                      "5"="Macrophages",
                      "6"="Macrophages", 
                      "7"= "B cells", 
                      "8"= "Neutrophils",
                      "9"="Macrophages",
                      "10"="Macrophages",
                      "11"="Neutrophils",
                      "12"="Macrophages",
                      "13"="Macrophages",
                      "14"="Macrophages",
                      "15"="Macrophages",
                      "16"="Macrophages",
                      "17"="Macrophages",
                      "18"="Monocytes",
                      "19"="Neutrophils",
                      "20"="Monocytes",
                      "21"="Macrophages",
                      "22"="Macrophages",
                      "23"="NK cells",
                      "24"="Macrophages",
                      "25"="T cells",
                      "26"="Macrophages",
                      "27"="Neutrophils",
                      "28"="Macrophages",
                      "29"="Macrophages",
                      "30"="Monocytes",
                      "31"="Macrophages",
                      "32"="Macrophages",
                      "33"="Endothelial cells",
                      "34"="Endothelial cells",
                      "35"="Fibroblasts",
                      "36"="Macrophages",
                      "37"="ILC",
                      "38"="Macrophages",
                      "39"="Macrophages",
                      "40"="Macrophages")
MI.integrated.1st.filtered[['cell_type']] = unname(cluster2celltype[MI.integrated.1st.filtered@meta.data$seurat_clusters])

#保存质控好的数据
saveRDS(MI.integrated.1st.filtered, "sce.all.rds")




colnames(MI.integrated.1st.filtered@meta.data)
VlnPlot(MI.integrated.1st.filtered,group.by ='cell_type',features = c("Il1rl1"))

cd4_sce2 = MI.integrated.1st.filtered[, Idents(MI.integrated.1st.filtered) %in% c( "Naive CD4 T" , "Memory CD4 T" )]
# subset 函数也可以
sel.clust = "cell_type"
MI.integrated.1st.filtered <- SetIdent(MI.integrated.1st.filtered, value = sel.clust)

Macrophages = MI.integrated.1st.filtered[,MI.integrated.1st.filtered@meta.data$cell_type %in% c("Macrophages")]


VlnPlot(Macrophages,group.by ='orig.ident',features = c("Il1rl1"),pt.size = 0)

Macrophages@

expr <- GetAssayData(object = Macrophages, slot = "counts")
pbmc[["RNA"]]$counts
expr = as.matrix(expr)
Macrophages
# 指定要提取的基因名称
gene_name <- "Il1rl1"
# 提取单个基因的表达矩阵


gene_expression = as.matrix( Macrophages[["RNA"]]$counts) 
rownames(gene_expression)

Macrophages@assays$RNA@counts



gene_expression <- Macrophages[["RNA"]]$counts.MI_1D.2
# 将表达矩阵保存为matrix格式
gene_expression <- as.data.frame(gene_expression)
gene_expression = as.data.frame(t(gene_expression))
meta = Macrophages@meta.data
identical(row.names(gene_expression),row.names(meta))
gene_expression$group = meta$orig.ident
table(gene_expression$group)


# ggpolt画图
data.long = gene_expression
library(ggplot2)
colnames(data.long)
# Multiple groups with error bars and jitter point
library(ggpubr)
table(data.long$group)
colnames(data.long)
max(data.long$Kdr)
ylim = 2
ggbarplot(data.long, x = 'group', y = 'Kdr', add = c('mean_se'),
          color = "group", fill = "group",palette =  c("#008ECA","#DB423E"),
          position = position_dodge(1),width = 0.5) +
  labs(y = "Kdr") + 
  theme( legend.position = 'none' ) +
  rotate_x_text(angle = 0)+
  coord_cartesian(ylim = c(0, ylim))+  # 设置y轴的坐标范围 
  stat_compare_means(method = 'wilcox.test', aes(group = group),label.y = ylim)
ggsave("4.diff.gene/01.CAP/geom_bar.Kdr.pdf",height = 5,width = 5,family="Times")




markers_all_RNA <- FindAllMarkers(object = MI.integrated.1st.filtered, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = "LR")
write.csv(markers_all_RNA, file = "MI_1st_filtered_cluster_allmarker_Genes_RNA_RES15.csv")




