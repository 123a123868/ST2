
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

write.csv(deg,"02.DEG/edgeR.DEG.online.csv")

# 转换基因
# org.Mm.eg.db
# org.Hs.eg.db
deg$symbol = rownames(deg)
df <- bitr(unique(deg$symbol), fromType = "SYMBOL",
           toType = c( "ENTREZID"),
           OrgDb = org.Hs.eg.db)
head(df)
DEG=deg
head(DEG)
DEG=merge(DEG,df,by.y='SYMBOL',by.x='symbol')
head(DEG)
table(DEG$g)
gene_up= DEG[DEG$g == 'UP','ENTREZID'] 
gene_down=DEG[DEG$g == 'DOWN','ENTREZID'] 
gene_diff=c(gene_up,gene_down)

# 画图
# go和kegg
# 全部基因的
# go分析
# go.all <- enrichGO(gene_diff, OrgDb = "org.Hs.eg.db", ont="all")   
# library(ggplot2)
# library(stringr)
# # barplot(go.all, split="ONTOLOGY")+ facet_grid(ONTOLOGY~., scale="free") 
# # barplot(go.all, split="ONTOLOGY",font.size =10)+ 
# #   facet_grid(ONTOLOGY~., scale="free") + 
# #   scale_y_discrete(labels=function(x) str_wrap(x, width=50))
# # ggsave("03-go_kegg/02.result/go_all.normal.pdf",height = 7,width = 10) 
# # # 保存结果
# # write.csv(go,file="03-go_kegg/02.result/go_all.FC1.csv",quote=F,row.names = F)             
# # 
# 
# 上调基因
gene_up=unique(gene_up)
gene_up
# go分析
go <- enrichGO(gene_up, OrgDb = "org.Hs.eg.db", ont="all")
library(ggplot2)
library(stringr)
barplot(go, split="ONTOLOGY")+ facet_grid(ONTOLOGY~., scale="free")
barplot(go, split="ONTOLOGY",font.size =10)+
  facet_grid(ONTOLOGY~., scale="free") +
  scale_y_discrete(labels=function(x) str_wrap(x, width=50))
ggsave("03.go.kegg/go_up.pdf",height = 7,width = 10)
# 保存结果
write.csv(go,file="03.go.kegg//go_up.csv",quote=F,row.names = F)



# 下调基因
gene_down=unique(gene_down)
gene_down
# go分析
go <- enrichGO(gene_down, OrgDb = "org.Hs.eg.db", ont="all")


library(ggplot2)
library(stringr)
barplot(go, split="ONTOLOGY")+ facet_grid(ONTOLOGY~., scale="free")
barplot(go, split="ONTOLOGY",font.size =10)+
  facet_grid(ONTOLOGY~., scale="free") +
  scale_y_discrete(labels=function(x) str_wrap(x, width=50))
ggsave("03.go.kegg//go_down.pdf",height = 10,width = 15)
# 保存结果
write.csv(go,file="03.go.kegg//go_down.FC1.csv",quote=F,row.names = F)



# 上下调画在一张图上
go.data = read.csv("03.go.kegg//GOall.csv")
colnames(go.data)
go.data$Count.FC = as.numeric(go.data$Count.FC )


ggplot(go.data, aes(x = reorder(Description,Count.FC), Count.FC,fill=group)) + 
  geom_bar(stat = 'identity',alpha = 0.7) + 
  coord_flip() + 
  theme_bw() + #去除背景色
  theme(panel.grid =element_blank())+
  theme(panel.border = element_rect(size = 0.6))+
  labs(x = "",
       y="count")+
  scale_fill_manual(values = c("#008ECA","#DB423E"))#设置颜色
ggsave("03.go.kegg//go.all.better.pdf",width = 8,height = 9,family="Times")

























