
# 设置
rm(list = ls()) 
options(stringsAsFactors = F)
library(clusterProfiler)
library(dplyr)
library(ggplot2)  
library(org.Hs.eg.db)
library(org.Mm.eg.db)

# 读入数据
load("02.DEG/expr.DEG_edgeR.Rdata")
## 设置阈值
deg = DEG_edgeR
deg$P.Value = deg$PValue
logFC_t= 0.58
deg$g=ifelse(deg$P.Value>0.05,'stable',
             ifelse( deg$logFC > logFC_t,'UP',
                     ifelse( deg$logFC < -logFC_t,'DOWN','stable') )
)
table(deg$g)
head(deg)
# 转换基因
deg$symbol = rownames(deg)
df <- bitr(unique(deg$symbol), fromType = "SYMBOL",
           toType = c( "ENTREZID"),
           OrgDb = org.Hs.eg.db)
head(df)
DEG=deg
head(DEG)
DEG=merge(DEG,df,by.y='SYMBOL',by.x='symbol')
head(DEG)
gene_up= DEG[DEG$g == 'UP','ENTREZID'] 
gene_down=DEG[DEG$g == 'DOWN','ENTREZID'] 
gene_diff=c(gene_up,gene_down)

write.csv(DEG,"02.DEG/DEG.online.csv")

# 画图 kegg
# kegg分析
# 所有的差异
# mmu   has
kk.all <- enrichKEGG(gene         = gene_diff,
                    organism     = 'has',
                    #universe     = gene_all,
                    pvalueCutoff = 0.9,
                    qvalueCutoff =0.9)
head(kk.all)[,1:6]
g_kegg = dotplot(kk.all)
g_kegg
ggsave(g_kegg,filename = "03-go_kegg/02.result/kegg.all.FC1.pdf")
write.csv(kk.all,file="03-go_kegg/02.result/kegg_all.FC1.csv",quote=F,row.names = F)  



###########################################################################################
# kegg分析
kk.up <- enrichKEGG(gene         = gene_up,
                    organism     = 'has',
                    #universe     = gene_all,
                    pvalueCutoff = 0.9,
                    qvalueCutoff =0.9)
head(kk.up)[,1:6]
kk=kk.up
g_kegg = dotplot(kk)
g_kegg
ggsave(g_kegg,filename = "03-go_kegg/02.result/kegg_up.FC1.pdf")
write.csv(kk.up,file="03-go_kegg/02.result/kegg_up.FC1.csv",quote=F,row.names = F)  

#########################################################################
# kegg分析
kk.down <- enrichKEGG(gene         = gene_down,
                    organism     = 'has',
                    #universe     = gene_all,
                    pvalueCutoff = 0.9,
                    qvalueCutoff =0.9)
head(kk.down)[,1:6]
kk=kk.down
g_kegg = dotplot(kk)
g_kegg
ggsave(g_kegg,filename = "03-go_kegg/02.result/kegg_down.FC1.pdf")
write.csv(kk,file="03-go_kegg/02.result/kegg_down.FC1.csv",quote=F,row.names = F)



# 上下调都画在一起。
# 棒棒图
kegg.data =read.csv("03-go_kegg/KEGGall.csv")
#多维度信息富集棒棒糖图绘制：
library(ggplot2)
library(cols4all)
colnames(kegg.data)
factor(kegg.data$Description)
level <- kegg.data$Description
kegg.data$Description = factor(kegg.data$Description,
                               levels = level)
factor(kegg.data$Description)


p4 <- ggplot(kegg.data,aes(x = Count,y = Description)) +
  geom_col(aes(fill = group), width = 0.1) +
  geom_point(aes(size = Count,
                 color = group)) +
  scale_color_manual(values = c("#008ECA","#DB423E"))+
  scale_size_continuous(range = c(2, 7))  +
  theme_classic()+
  scale_fill_manual(values = c("#008ECA","#DB423E"))#设置颜色
p4
ggsave("03-go_kegg/KEGG.all.better.pdf",width = 8,height = 9)

p4




