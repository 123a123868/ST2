# 细胞免疫浸润

# 设置
rm(list = ls())
options(stringsAsFactors = F)
library(GSVA)
library(limma)
library(GSEABase)
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
expr=log2(expr+1)
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
# save(expr,data_gene,file = '02.DEG//expr.deg.Rdata')


# 读入免疫细胞浸润数据
gmtFile <- "05.immune/Immune_cell_infiltration/immune.gmt"
dimnames <- list(rownames(expr),colnames(expr))
mat <- matrix(as.numeric(as.matrix(expr)),nrow=nrow(expr),dimnames=dimnames)
mat <- avereps(mat)
mat <- mat[rowMeans(mat)>0,]
geneSet=getGmt(gmtFile, 
               geneIdType=SymbolIdentifier())

# do gsva
ssgseaScore <- gsva(mat, geneSet, 
                    method='ssgsea', 
                    kcdf='Gaussian', abs.ranking=TRUE)
normalize=function(x){
  return((x-min(x))/(max(x)-min(x)))}
ssgseaOut <- normalize(ssgseaScore)
ssgseaOut <- rbind(id=colnames(ssgseaOut),ssgseaOut)
dim(ssgseaOut)
ssgseaOut[1:4,1:4]

# save data----
write.table(ssgseaOut,file="05.immune//ssgseaOut.txt",sep="\t",quote=F,col.names=F)


################################################################################
# 可视化画图
library(pheatmap) 
library(ggrisk)  
library(pheatmap) 
library(ggplot2)
library(ggsci)
library(ggstatsplot)
library(data.table)

# 绘图数据
# rm(list = ls())  
options(stringsAsFactors = F) 
input <- "05.immune//ssgseaOut.txt"  
immune <- read.table( input ,sep="\t",header=T,row.names=1,check.names=F)
immune[1:4, 1:4]

# 制作ST2.L和ST2.H组的
colnames(immune)


Type = data_gene
identical(rownames(Type),colnames(immune))



#  vioplot----
outTab=data.frame()
table(Type$group)
ST2.L = 7
ST2.H = 10
immune <- data.frame(t(immune),check.names = F)

pdf("05.immune/ssgseaScore.ST2.pdf",height=7,width=15)
par(las=1,mar=c(10,6,3,3))
x=c(1:ncol(immune))
y=c(1:ncol(immune))
plot(x,y,
     xlim=c(0,72),ylim=c(min(immune),max(immune)+0.02),
     main="",xlab="", ylab="Fraction",
     pch=24,
     col="white",
     xaxt="n")

library(vioplot)
pFilter=0.05
for(i in 1:ncol(immune)){
  #i = 1
  if(sd(immune[1:ST2.L,i])==0){
    immune[1,i]=0.001
  }
  if(sd(immune[(ST2.L+1):(ST2.L+ST2.H),i])==0){
    immune[(ST2.L+1),i]=0.001
  }
  lowData=immune[1:ST2.L,i]
  highData=immune[(ST2.L+1):(ST2.L+ST2.H),i]
  vioplot(lowData,at=3*(i-1),lty=1,add = T,drawRect = F,col = '#00A087')
  vioplot(highData,at=3*(i-1)+1,lty=1,add = T,drawRect = F,col = '#E64B35')
  wilcoxTest=t.test(lowData,highData)
  p=wilcoxTest$p.value
  if(p<pFilter){
    cellPvalue=cbind(Cell=colnames(immune)[i],pvalue=p)
    outTab=rbind(outTab,cellPvalue)
  }
  mx=max(c(lowData,highData))
  lines(c(x=3*(i-1)+0.2,x=3*(i-1)+0.8),c(mx,mx))
  text(x=3*(i-1)+0.5, y=mx+0.02, labels=ifelse(p<0.001, paste0("p<0.001"), paste0("p=",sprintf("%.03f",p))), cex = 0.8)
}
legend("topright", 
       c("ST2.L", "ST2.H"),
       lwd=5,bty="n",cex=1.5,
       col=c("#00A087","#E64B35"))
text(seq(1,70,3),-0.1,xpd = NA,labels=colnames(immune),cex = 1,srt = 45,pos=2)
dev.off()



