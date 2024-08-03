# 上下调都画在一起。
# 棒棒图
kegg.data =read.csv("kegg.csv")
#多维度信息富集棒棒糖图绘制：
library(ggplot2)
library(cols4all)
colnames(kegg.data)
colnames(kegg.data)[1] = "Drug"


barplot(kegg.data$Combined.Score)


colnames(kegg.data)
kegg.data = kegg.data[order(kegg.data$Combined.Score,decreasing = T,na.last = NA),]
level <- kegg.data$Drug
kegg.data$Drug = factor(kegg.data$Drug,
                               levels = level)
factor(kegg.data$Drug)

# 第一张图绘图

ggplot(kegg.data,aes(Drug,Combined.Score,
                     fill= Adjusted.P.value))+
geom_col()+
coord_flip()+
labs(x="Drug",y="Combined.Score")+
scale_fill_gradient(low = "red",high="green")+
scale_y_continuous(limits = c(0,1300),breaks = seq(0,1300,500))+
theme(text = element_text("A",size = 15,face = "bold"))+theme_classic()

ggsave('drug.pdf',width=6,height=7,family = "Times")

















factor(kegg.data$Description)
level <- kegg.data$Description
kegg.data$Description = factor(kegg.data$Description,
                               levels = level)
factor(kegg.data$Description)
kegg.data$group = "up"
colnames(kegg.data)
p4 <- ggplot(kegg.data,aes(x = Combined.Score,y = Drug)) +
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

