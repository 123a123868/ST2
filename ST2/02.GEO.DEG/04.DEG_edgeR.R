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
expr$gene = row.names(expr)
colnames(expr)
expr = expr[,c("gene",colnames(expr)[1:17])]
expr[1:4,1:4] 
write.table(expr,"05.immune/online.imm.expr.txt")

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


# 创建一个数据框用于存储样本的分组信息，行名为样本名，列名为分组信息
identical(colnames(expr),rownames(data_gene))
colData <- data.frame(row.names = colnames(expr),
                      group = group)
colData$group <- factor(colData$group, levels = c("ST2.H", "ST2.L"))  
head(colData)
# expr 取整数
expr = floor(expr)


library(edgeR)
# 创建 DGEList 对象，用于存储基因表达数据和组信息，还是使用原始计数矩阵
d <- DGEList(counts = expr, group = group)
# 根据每个基因在每个样本中的 CPM（Counts Per Million）值去除低表达基因
keep <- rowSums(cpm(d) > 1) >= 2
# 或者自动过滤，去除低表达基因
# keep <- filterByExpr(d)
table(keep)
# 从 DGEList 对象中筛选出符合条件的基因
d <- d[keep, , keep.lib.sizes = FALSE]
# 更新样本的库大小信息
d$samples$lib.size <- colSums(d$counts)
# # 归一化，TMM 方法
# d <- calcNormFactors(d)
# # 查看归一化后的样本信息
# head(d$samples)
# # 将归一化后的数据赋值给 dge 变量
dge = d

# 创建设计矩阵，用于指定差异分析模型
design <- model.matrix(~0 + factor(group))
rownames(design) <- colnames(dge)
colnames(design) <- levels(factor(group))

# 估计数据的离散度 —— common离散度、trended离散度、tagwise离散度
dge <- estimateGLMCommonDisp(dge, design)
dge <- estimateGLMTrendedDisp(dge, design)
dge <- estimateGLMTagwiseDisp(dge, design)

# 在估计的模型基础上进行 广义线性模型 (GLM) 拟合
fit <- glmFit(dge, design)


# edgrR 涉及到差异表达分析的函数有很多： exactTest、glmFit、glmLRT、glmQLFit、glmQLFTest。 
# qCML估计离散度需要搭配 exact test 进行差异表达分析，对应 exactTest 函数。 
# 而其他四个glm都是与GLM模型搭配使用的函数。其中，glmFit 和 glmLRT 函数是配对使用的，
# 用于 likelihood ratio test (似然比检验)，而 glmQLFit和 glmQLFTest则配对用于 quasi-likelihood F test (拟极大似然F检验)。

# 使用 LRT（Likelihood Ratio Test）计算差异表达
# 注意这里的 contrast 和 DESeq2 不一样，这里我们只需要输入 c(-1, 1) 即可
# -1 对应 normal，1 对应 tumor
lrt <- glmLRT(fit, contrast = c(1, -1))
# 从 LRT 计算结果中获取前 nrow(dge) 个顶部差异表达基因
nrDEG_edgeR <- topTags(lrt, n = nrow(dge))
# 将差异表达基因结果转换为数据框形式
DEG_edgeR_edgeR <- as.data.frame(nrDEG_edgeR)
# 输出差异表达基因结果的前几行
head(DEG_edgeR_edgeR)


# 设置上调下调分组
DEG_edgeR_edgeR$group=ifelse(DEG_edgeR_edgeR$PValue>0.05,'stable', 
                 ifelse( DEG_edgeR_edgeR$logFC >0.58,'up', 
                         ifelse( DEG_edgeR_edgeR$logFC < -0.58,'down','stable')))
table(DEG_edgeR_edgeR$group)
DEG_edgeR_edgeR["IL1RL1",]
DEG_edgeR = DEG_edgeR_edgeR
save(DEG_edgeR,data_gene,file = '02.DEG/expr.DEG_edgeR.Rdata')


# 检查是否分组分错
library(ggplot2)
library(ggpubr)
exprSet= expr
DEG_edgeR = DEG_edgeR_edgeR[order(DEG_edgeR_edgeR$logFC,decreasing = T,na.last = NA),]
gene.cheak = rownames(DEG_edgeR)[1]
gene.cheak = "IL1RL1"
data.check = as.data.frame(exprSet[gene.cheak,])
data.check = as.data.frame(t(data.check))
identical(row.names(data_gene),rownames(data.check))
data.check$group = data_gene$group
colnames(data.check)[1] = "gene"
# 绘制箱线图
p1 <- ggplot(data.check,aes(x=group,y=gene,fill=group))+
  geom_boxplot(width=0.6,alpha=0.8)+ stat_compare_means()+
  ggtitle(gene.cheak)
p1
ggsave("02.DEG_edgeR//check.group.pdf",height = 7,width = 6,family="Times")
save(DEG_edgeR,file = '02.DEG_edgeR//DEG_edgeR.Rdata')
write.csv(DEG_edgeR,file = '02.DEG_edgeR//DEG_edgeR.csv')
