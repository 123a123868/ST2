# WGCNA + 差异分析的交集
# 筛选hub基因
# 目的：差异分析
rm(list = ls())  
options(stringsAsFactors = F)
library(AnnoProbe)
library(GEOquery) 
library(openxlsx)

# 读入数据
expr = read.csv("Input_data/AllSamples_Genes_FPKM.csv")
dim(expr)
colnames(expr)[1] = "gene_id"

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
expro.upper = datExpr0[,which(m.vars > quantile(m.vars, probs = seq(0, 1, 0.25))[4])]
datExpr1<-data.matrix(expro.upper)
dim(datExpr1)
dim(expr)

# 先检查是否有哪个sample或基因表达量缺失，定义成不好的sample或gene
library(WGCNA)
gsg = goodSamplesGenes(datExpr1, verbose = 5)
gsg$allOK

# 如果你用全部基因作为输入，很有可能返回FALSE，说明存在不好的基因或sample。
# 下面的代码就会去除那些不好的基因或sample。
# if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(datExpr1)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr1 = datExpr1[gsg$goodSamples, gsg$goodGenes]
# }

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
#  datExpr = as.data.frame(datExpr1)
#  nGenes = ncol(datExpr)
#  nSamples = nrow(datExpr)


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
# mergeCutHeight 这个是融合的阈值，调大一点吧
cor <- WGCNA::cor  #，是WGCNA包与其他函数冲突导致的
net = blockwiseModules(datExpr, power = 16,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.3,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "MyTOM",
                       verbose = 3)
table(net$colors)
MEs = net$MEs
cor<-stats::cor  #后面再调回来
# gene module的可视化
mergedColors = labels2colors(net$colors)
pdf("06-WGCNA/2module.pdf",width = 10, height = 5)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]], "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

# 输出每个module内的基因，用于后续分析作图，随便你DIY
moduleColors = labels2colors(net$colors)
color<-unique(moduleColors)
for (i  in 1:length(color)) {
  y=t(assign(paste(color[i],"expr",sep = "."),datExpr[moduleColors==color[i]]))
  write.csv(y,paste("06-WGCNA/model/",color[i],"csv",sep = "."),quote = F)
}

# gene module的深入分析
# 通过计算表型与module的相关性，找出相关性高的gene module，推测可能是因为它们造成了表型差异。
colnames(expr)
samples=read.csv('06-WGCNA/Sam_info_mouse.txt',sep = '\t',row.names = 1)
moduleLabelsAutomatic = net$colors
moduleColorsAutomatic = labels2colors(moduleLabelsAutomatic)
moduleColorsWW = moduleColorsAutomatic
MEs0 = moduleEigengenes(datExpr, moduleColorsWW)$eigengenes
MEsWW = orderMEs(MEs0)
modTraitCor = cor(MEsWW, samples, use = "p")
colnames(MEsWW)


modlues=MEsWW
modTraitP = corPvalueStudent(modTraitCor, nSamples)
textMatrix = paste(signif(modTraitCor, 2), "\n(", signif(modTraitP, 1), ")", sep = "")
dim(textMatrix) = dim(modTraitCor)

pdf("06-WGCNA/4Module-trait.pdf",width = 6, height = 6)
labeledHeatmap(Matrix = modTraitCor, xLabels = colnames(samples), yLabels = names(MEsWW), cex.lab = 0.5,  yColorWidth=0.01, 
               xColorWidth = 0.03,
               ySymbols = colnames(modlues), colorLabels = FALSE, colors = blueWhiteRed(50), 
               textMatrix = textMatrix, setStdMargins = FALSE, cex.text = 0.5, zlim = c(-1,1)
               , main = paste("Module-trait relationships"))
dev.off()

# 找差异基因和 hub基因的交集
hub_gene <-read.csv("06-WGCNA/model/.yellow.csv",row.names = 1)
dim(hub_gene)  # 黄色模块一共有400多个基因

#3.1.首先计算模块特征值(module eigengenes)
MEs = moduleEigengenes(datExpr, moduleColors)$eigengenes
#3.2.计算module membership即MM
datKME =signedKME(datExpr, MEs, outputColumnName="MM.")
#3.3.gene significance即GS
GS = as.data.frame(cor(datExpr, samples, use = "p")) #datTraits为目标性状的矩阵文件
#3.4.循环提取每个模块的hubgene
MMname = colnames(datKME)
for (mm in MMname){
  # mm = MMname[1]
  FilterGenes = abs(GS)> 0.5 & abs(datKME$MM.yellow)> 0.9
  hubgenes = dimnames(data.frame(datExpr))[[2]][FilterGenes]  
}
hubgenes  
hubgenes = na.omit(hubgenes)


# 差异基因
deg = read.csv("02-DEG/TN_PEI//DEG.TN_PEI.csv")
row.names(deg) = deg$X
head(deg)
## 设置阈值
logFC_t= 1
deg$g=ifelse(deg$P.Value>0.05,'stable',
             ifelse( deg$logFC > logFC_t,'UP',
                     ifelse( deg$logFC < -logFC_t,'DOWN','stable') )
)
table(deg$g)
colnames(deg)
gene_up= deg[deg$g == 'UP','X'] 
gene_down=deg[deg$g == 'DOWN','X'] 
gene_diff=c(gene_up,gene_down)
kp = intersect(hubgenes,gene_diff)
kp
----------------------------------------------------------------------------------------------------
# 设置数据集的名称和内容
A <- hubgenes
B <- gene_diff

# 画维恩图
venn.plot = venn.diagram(
  x=list(A=A, B=B), # 数据集
  filename= NULL,
  col=c("blue", "red"), # 设置颜色
  fill=c("dodgerblue", "pink"), # 设置填充色
  cat.col=c("blue", "red"), # 设置分类名称颜色
  cat.cex=1.5, # 设置分类名称大小
  cat.fontface="bold", # 设置分类名称字体
  cat.dist=c(.06, -.06), # 设置分类名称间距
  # cat.pos=c(,), # 设置分类名称位置
  lty=2, # 设置维恩图线型
  main="Venn diagram of A and B", # 设置图标题
  lwd=2, # 设置线宽
  fontfamily="serif", # 设置字体
)
pdf(file="06-WGCNA/venn.pdf")
grid.draw(venn.plot)
dev.off()



















  
moduleColors = labels2colors(net_power$colors) #net_power为一步法构建出的对象
# 基于准备工作，进行下列计算
# 计算KME值
#第二种：直接输入表达矩阵计算
power = 6
kIM <- intramodularConnectivity.fromExpr(datExpr, colors, power = power)

power=6
moduleColors = labels2colors(net_power$colors) #net_power为一步法构建出的对象





#第一种： 先计算邻接矩阵，再计算连通度，推荐此种,邻接矩阵(Adjacency matrix)指基因和基因之间的加权相关性值取power次方即无尺度化之后构成的矩阵。
power = 16
adjacency = abs(cor(datExpr,use="p"))^power #datExpr为表达矩阵，power请与一步法中的软阈值保持一致
kIM <- intramodularConnectivity(adjacency, moduleColors)
head(kIM)
table(moduleColors)
hub = chooseTopHubInEachModule(datExpr,moduleColors, omitColors = "yellow", power = power)
hub



  

MEs = net$colors
datKME = signedKME(datExpr, MEs, outputColumnName="kME_MM.")
write.csv(datKME, "06-WGCNA/kME_MM.csv")
# 提取感兴趣的module，这里以lightyellow为例
FilterGenes = abs(datKME$MM.yellow) > 0.8
# 返回结果
table(FilterGenes)
# 共有41个gene，输出结果
hubgene_blue <- dimnames(data.frame(datExpr))[[2]][FilterGenes]
hubgene_blue

  
  
  
  
  
  
  

 