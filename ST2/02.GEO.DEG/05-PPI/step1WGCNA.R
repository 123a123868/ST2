# WGCNA + 差异分析的交集
# 筛选hub基因
# 目的：差异分析
rm(list = ls())  
options(stringsAsFactors = F)
library(AnnoProbe)
library(GEOquery) 
library(openxlsx)

# 读入数据
expr = read.csv("Input_data/AllSamples_Genes_ReadsCount.csv")
dim(expr)
colnames(expr)[1] = "gene_id"

# 转换一下id
ids = read.csv("Input_data/AllSamples_Genes_FPKM.csv")
colnames(ids)[1] = "gene_id"
table(expr$gene_id %in% ids$gene_id)
identical(expr$gene_id,ids$gene_id) # 完全一致，可以转换
expr$SYMBOL = ids$SYMBOL


# 去除重复值
expr$median=apply(expr[,2:11],1,median)   #expr新建median这一列，列名为median，同时对dat这个矩阵按行操作，取每一行的中位数，将结果给到median这一列的每一行
expr=expr[order(expr$SYMBOL,expr$median,decreasing = T),]   #对expr$symbol按照expr$median中位数从大到小排列的顺序排序，将对应的行赋值为一个新的expr
expr=expr[!duplicated(expr$SYMBOL),] #将symbol这一列取取出重复项，'!'为否，即取出不重复的项，去除重复的gene ，保留每个基因最大表达量结果s
dim(expr)
expr = dplyr::filter(expr, !is.na(SYMBOL))
rownames(expr) = expr$SYMBOL
expr = expr[,2:11]
expr[1:4,1:4]
# 手动看过的了，这部分没有搞乱，正确

# 标准化
# 需要log
boxplot(expr[,1:4],las=2) 
library(limma)
expr = log(expr+1)
expr=normalizeBetweenArrays(expr)
boxplot(expr[,1:4],las=2) 


colnames(expr)
# 转置，也就是行变成列，列变成行
datExpr0 = data.frame(t(expr))
colnames(datExpr0) <- rownames(expr)
rownames(datExpr0) <- colnames(expr)
dim(datExpr0)
datExpr0[1:4,1:4]

# 筛选方差前25%的基因
# 这部分可以不做。用前25%的基因会比用全部基因找出来的module数量少。
# 这个还是要筛选的。筛选方差前 90%的
datExpr1<-datExpr0

m.vars=apply(datExpr0,2,var)
expro.upper = datExpr0[,which(m.vars > quantile(m.vars, probs = seq(0, 1, 0.1))[4])]
datExpr1<-data.matrix(expro.upper)
dim(datExpr1)
dim(expr)



# 先检查是否有哪个sample或基因表达量缺失，定义成不好的sample或gene
library(WGCNA)
gsg = goodSamplesGenes(datExpr1, verbose = 5)
gsg$allOK

# 如果你用全部基因作为输入，很有可能返回FALSE，说明存在不好的基因或sample。
# 下面的代码就会去除那些不好的基因或sample。
if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(datExpr1)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr1 = datExpr1[gsg$goodSamples, gsg$goodGenes]
}

# 判断是否有离群样本
# 通过聚类分析，能看出是否有个别sample跟多数sample的距离较远，决定是否需要去除离群样本。
sampleTree = hclust(dist(datExpr1), method = "average")
par(cex = 0.7);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
# 例如：以35000作为cutoff，就会去除最左边的4四个sample，只剩下16个sample。
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2) +
  #想用哪里切，就把“h = 35000”和“cutHeight = 35000”中的“500”换成你的cutoff
  abline(h = 100, col = "red") 
clust = cutreeStatic(sampleTree, cutHeight = 1000000, minSize = 10)
keepSamples = (clust==1)
datExpr = datExpr1[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
dim(datExpr)

# 这里不切割
# 运行下面的代码，用全部10个样本进行后续的分析：
datExpr = as.data.frame(datExpr1)
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)


# 找gene module
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
pdf("06-WGCNA/2Threshold.pdf",width = 10, height = 5)
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence")) +
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red")+
  abline(h=0.90,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity")) +
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()
softPower = sft$powerEstimate 
softPower
# 从上面的结果可以看出，从 10 开始进入“平台期”。因此，我们把下面代码里的power设置为power = 10
# 构建网络，找出gene module
net = blockwiseModules(datExpr, power = 10,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "MyTOM",
                       verbose = 3)
table(net$colors)


dynamicColors = labels2colors (dynamicMods) #在树状图下绘制模块颜色分配。注:灰色为unasunsigned基因。
dynamicColors = labels2colors(dynamicMods) 

# gene module的可视化
mergedColors = labels2colors(net$colors)
pdf("06-WGCNA/2module.pdf",width = 10, height = 5)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]], "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()


MEList = moduleEigengenes (datExpr, colors = dynamicColors) 
#按照模块计算每个module的ME（也就是该模块的第一主成分）





# Calculate eigengenes 
MEList = moduleEigengenes(datExpr, colors = dynamicColors) 
MEs = MEList$eigengenes 
# Calculate dissimilarity of module eigengenes 
MEDiss = 1-cor(MEs); 
# Cluster module eigengenes 
METree = hclust(as.dist(MEDiss), method = "average") 
# Plot the result 
#sizeGrWindow(7, 6) 
pdf(file="6_Clustering of module eigengenes.pdf",width=7,height=6) 
plot(METree, main = "Clustering of module eigengenes", 
     xlab = "", sub = "") 
MEDissThres = 0.4######剪切高度可修改 
# Plot the cut line into the dendrogram 
abline(h=MEDissThres, col = "red") 
dev.off()


# Call an automatic merging function 
merge = mergeCloseModules(datExpr0, dynamicColors, cutHeight = MEDissThres, verbose = 3) 
# The merged module colors 
mergedColors = merge$colors 
# Eigengenes of the new merged modules: 
mergedMEs = merge$newMEs 
table(mergedColors) 

#sizeGrWindow(12, 9) 
pdf(file="7_merged dynamic.pdf", width = 9, height = 6) 
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), 
                    c("Dynamic Tree Cut", "Merged dynamic"), 
                    dendroLabels = FALSE, hang = 0.03, 
                    addGuide = TRUE, guideHang = 0.05) 
dev.off()



yellow = read.csv("06-WGCNA/model/.yellow.csv")

deg = read.csv("02-DEG/TN_PEI//DEG.TN_PEI.csv")
rownames(deg) = deg$X
table(deg$group)

gene_down = rownames(deg[deg$group == "down",])
gene_up = rownames(deg[deg$group == "up",])
gene_all = c(gene_down,gene_up)

intersect(gene_all,yellow$X)
















# module的可重复性reproducible
# 1. 建立训练集和测试集
library(caret)
inTraining <- createDataPartition(datExpr$Zswim5, p = 0.5, list = FALSE)  # 随机分隔数据
inTraining
library(dplyr)
train<- datExpr[inTraining,]
test<-datExpr[-inTraining,]
# 删除全是0的一列
train = train[,which(colSums(train) > 0)]
test = test[,which(colSums(test) > 0)]
kp = intersect(colnames(train),colnames(test))
train = train[,kp]
test = test[,kp]
identical(colnames(train),colnames(test))

setLabels = c("Train", "Test")
multiExpr = list(Train = list(data = train), Test = list(data = test));
moduleColors = labels2colors(net$colors)
multiColor = list(Train = moduleColors)

# 2. preservation分析
# nPermutations官网上给了200，此处为节省时间，nPermutations设置为20
mp = modulePreservation(multiExpr, multiColor,
                        referenceNetworks = 1,
                        nPermutations = 20,
                        randomSeed = 1,
                        quickCor = 0,
                        verbose = 5)











