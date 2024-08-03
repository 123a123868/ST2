####
######
# 火山图
# 功能分析

# 设置
rm(list = ls()) 
options(stringsAsFactors = F)
library(limma)
library(edgeR)

# 读入数据
load("02.DEG/expr.DEG_edgeR.Rdata")
deg = DEG_edgeR
title = "Nornal vs AMI"

# 设置上调下调分组
deg$group=ifelse(deg$PValue>0.05,'stable', 
                 ifelse( deg$logFC > 0.58,'up', 
                         ifelse( deg$logFC < -0.58,'down','stable')))
table(deg$group)

# deg$FC = 2^(deg$logFC)


# 画热图
# 画火山图
library(ggplot2)
library(ggrepel)
library(ggthemes)

# 设置颜色
#plot_mode <- "classic" #经典版
plot_mode <- "advanced" #酷炫版

# 设置分界线
logFCcut <- 0.58 
logFCcut2 <- 0.58 #for advanced mode
logFCcut3 <- 0.58 #for advanced mode

pvalCut <- 0.05
pvalCut2 <- 0.05 #for advanced mode
pvalCut3 <- 0.05 #for advanced mode

adjPcut <- 0.05

x = deg
x$label =  rownames(x)
x$P.Value = x$PValue
x$logFC = x$logFC

head(x)
# re.adj <- x[x$adj.P.Val <adjPcut & (x$logFC > logFCcut | x$logFC < -logFCcut),] 
# 选取p.value < 0.05，且|logFC|>1的基因
sum(x$P.Val < pvalCut)
re.p = x[x$P.Val <pvalCut & (x$logFC > logFCcut | x$logFC < -logFCcut),]
if (plot_mode == "classic"){
  # 簡單的setting for color
  x[,6] <- ifelse((x$P.Value < pvalCut & x$logFC > logFCcut), "red", ifelse((x$P.Value < pvalCut & x$logFC < -logFCcut), "blue","grey30"))
  # 簡單的setting for size
  size <- ifelse((x$P.Value < pvalCut & abs(x$logFC) > logFCcut), 4, 2)
  
} else if (plot_mode == "advanced") {
  # 複雜的的setting for color
  n1 <- length(x[,1])
  cols <- rep("grey30",n1)
  names(cols)<- rownames(x)
  cols[x$P.Value < pvalCut & x$logFC >logFCcut]<- "#FB9A99"
  cols[x$P.Value < pvalCut2 & x$logFC > logFCcut2]<- "#ED4F4F"
  cols[x$P.Value < pvalCut & x$logFC < -logFCcut]<- "#B2DF8A"
  cols[x$P.Value < pvalCut2 & x$logFC < -logFCcut2]<- "#329E3F"
  #cols[names(cols)==genelist]<- "red"
  color_transparent <- adjustcolor(cols, alpha.f = 0.5)  # set color transparence
  x[,6] <- color_transparent
  
  # 複雜的的setting for size
  n1 <- length(x[,1])
  size <- rep(1,n1)
  size[x$P.Value < pvalCut & x$logFC > logFCcut]<- 2
  size[x$P.Value < pvalCut2 & x$logFC > logFCcut2]<- 4
  size[x$P.Value < pvalCut3 & x$logFC > logFCcut3]<- 6
  size[x$P.Value < pvalCut & x$logFC < -logFCcut]<- 2
  size[x$P.Value < pvalCut2 & x$logFC < -logFCcut2]<- 4
  size[x$P.Value < pvalCut3 & x$logFC < -logFCcut3]<- 6
  
} else {
  stop("Unsupport mode")
}

# 畫圖位置x，y軸的最大最小位置
xmin <- (range(x$logFC)[1]-(range(x$logFC)[1]+10))
xmax <- (range(x$logFC)[1]+(10-range(x$logFC)[1]))
ymin <- 0
ymax <- 5
xmin <- -3
xmax <- 3

#ymax <- range(-log(x$P.Val))[2]+1
min(size)
max(size)
size = size/2  # 修改一下size的大小
# Construct the plot object
p1 <- ggplot(data=x, aes(x=logFC, y=-log10(P.Value), colour = group)) +
  geom_point(alpha = 0.6, size=size, colour=x[,6]) +
  scale_color_manual(values = c("lightgrey", "navy", "red")) + 
  
  labs(x=bquote(~Log[2]~"(fold change)"), y=bquote(~-Log[10]~italic("P-value")), title= title) + 
  ylim(c(ymin,ymax)) + 
  scale_x_continuous(
    breaks = c(-5, -logFCcut, 0, logFCcut, 5), #刻度线的位置
    labels = c( -5, -logFCcut, 0, logFCcut, 5),
    limits = c(xmin, xmax) #x轴范围，两侧对称才好看
  ) +
  #或用下面这行：
  #xlim(c(xmin, xmax)) + 
  
  #画阈值分界线
  geom_vline(xintercept = logFCcut, color="grey40", linetype="longdash", size=0.5) +
  geom_vline(xintercept = -logFCcut, color="grey40", linetype="longdash", size=0.5) +
  geom_hline(yintercept = -log10(pvalCut), color="grey40", linetype="longdash", size=0.5) +
  
  guides(colour = guide_legend(override.aes = list(shape=16)))+
  
  theme_bw(base_size = 12, base_family = "Times") +
  theme(legend.position="right",
        panel.grid=element_blank(),
        legend.title = element_blank(),
        legend.text= element_text(face="bold", color="black",family = "Times", size=8),
        plot.title = element_text(hjust = 0.5),  # 调标题 的位置
        axis.text.x = element_text(face="bold", color="black", size=16),
        axis.text.y = element_text(face="bold",  color="black", size=16),
        axis.title.x = element_text(face="bold", color="black", size=16),
        axis.title.y = element_text(face="bold",color="black", size=16))
if (plot_mode == "advanced") {
  p1 <- p1 + 
    geom_vline(xintercept = logFCcut2, color="grey40", linetype="longdash", size=0.5) +
    geom_vline(xintercept = -logFCcut2, color="grey40", linetype="longdash", size=0.5) +
    geom_hline(yintercept = -log10(pvalCut2), color="grey40", linetype="longdash", size=0.5)
}
p1
ggsave("02.DEG//volcano2.pdf",height = 8,width = 8,family="Times",dpi = 100)

