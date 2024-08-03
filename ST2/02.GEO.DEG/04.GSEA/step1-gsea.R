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

#######
# GSEA分析
geneList = DEG$logFC
names(geneList)=DEG$ENTREZID 
geneList=sort(geneList,decreasing = T) # 降序，按照logFC的值来排序
geneList
## GSEA分析
library(enrichplot)


package.version("clusterProfiler")



kk_gse <- gseKEGG(geneList     = geneList,
                  organism     = 'has',
                  nPerm        = 1000,
                  minGSSize    = 10,
                  pvalueCutoff = 0.9,
                  verbose      = FALSE)
kk_gse=DOSE::setReadable(kk_gse, OrgDb='org.Mm.eg.db',keyType='ENTREZID')
sortkk<-kk_gse[order(kk_gse$enrichmentScore, decreasing = T),]
write.csv(sortkk,"04-gsea///gsea.csv",quote=F,row.names = F)  
## 展示排名前四的通路
head.kp = c(1:4)
tail.kp = c((length(row.names(sortkk))-4):length(row.names(sortkk)))
kp = c("mmu04612","mmu04650","mmu04662","mmu04814","mmu00030")


gseaplot2(kk_gse, kp,subplots = c(1,2))
ggsave("04-gsea/gsea.all.better3.pdf",height = 6,width = 8)
gseaplot2(kk_gse, row.names(sortkk)[tail.kp])
ggsave("02.DEG/02.EMT6/sg12//gsea.tail4.pdf")
write.csv(sortkk,"02.DEG/02.EMT6/sg12//gsea.csv",quote=F,row.names = F)  















# 这里找不到显著下调的通路，可以选择调整阈值，或者其它。
down_kegg<-kk_gse[kk_gse$pvalue<0.01 & kk_gse$enrichmentScore < -0.7,];down_kegg$group=-1
up_kegg<-kk_gse[kk_gse$pvalue<0.01 & kk_gse$enrichmentScore > 0.7,];up_kegg$group=1
dat=rbind(up_kegg,down_kegg)
colnames(dat)
dat$pvalue = -log10(dat$pvalue)
dat$pvalue=dat$pvalue*dat$group 
dat=dat[order(dat$pvalue,decreasing = F),]

g_kegg<- ggplot(dat, aes(x=reorder(Description,order(pvalue, decreasing = F)), y=pvalue, fill=group)) + 
  geom_bar(stat="identity") + 
  scale_fill_gradient(low="blue",high="red",guide = FALSE) + 
  scale_x_discrete(name ="Pathway names") +
  scale_y_continuous(name ="log10P-value") +
  coord_flip() + theme_bw(base_size = 15)+
  theme(plot.title = element_text(hjust = 0.5),  axis.text.y = element_text(size = 15))+
  ggtitle("Pathway Enrichment") 
g_kegg
print(g_kegg)
ggsave(g_kegg,filename = paste0(pro,'_kegg_gsea.png'))


library(enrichplot)
gesa_res=kk@result
gesa_res
lapply(1:nrow(down_kegg), function(i){ 
  gseaplot2(kk,down_kegg$ID[i],
            title=down_kegg$Description[i],pvalue_table = T)
  ggsave(paste0(pro,'_down_kegg_',
                gsub('/','-',down_kegg$Description[i])
                ,'.pdf'))
})
lapply(1:nrow(up_kegg), function(i){ 
  gseaplot2(kk,up_kegg$ID[i],
            title=up_kegg$Description[i],pvalue_table = T)
  ggsave(paste0(pro,'_up_kegg_',
                gsub('/','-',up_kegg$Description[i]),
                '.pdf'))
})






