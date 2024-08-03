# 先读入数据
# 确定号预后分组
# 画一个小提琴图


# 设置
rm(list = ls()) 
options(stringsAsFactors = F)


#读入数据
library(openxlsx)
data.all = read.xlsx('00.RawData/ST2data.xlsx',sheet=1)
table(data.all$prognosis)

# 先确定预后的筛选标准
# 三点标准
#满足一个条件
library(dplyr)
colnames(data.all)
data.bed = filter(data.all,全因死亡 == "yes" 
                  | 非致死性复发性心肌梗死 ==  "yes"
                  | 重复冠状动脉血运重建和卒中事件 == "yes")  

data.all$MACE.3 = ifelse(data.all$Patient_id %in% data.bed$Patient_id,"bed","good")
table(data.all$MACE.3)
# 四点标准
colnames(data.all)
data.bed = filter(data.all,
                    全因死亡 == "yes" 
                  | 非致死性复发性心肌梗死 ==  "yes"
                  | 重复冠状动脉血运重建和卒中事件 == "yes"
                  | 因不稳定型心绞痛或心力衰竭再住院 == "yes")  
data.all$MACE.4 = ifelse(data.all$Patient_id %in% data.bed$Patient_id,"bed","good")
table(data.all$MACE.4)

# 筛选一些掉一部分数据
write.csv(data.all,"01.ROC/01.filterData/filter.csv")

# 再读进来看看什么情况
data.all = read.csv("01.ROC/01.filterData/filter.3.csv")


table(data.all$MACE性别)


# 先画一个箱形图，看看有没有显著性
# 开始画图
# 循环画多个小提琴图
library(ggpubr)
data = data.all
data$group = data$MACE.4
table(data$group)
data$group=factor(data$group, levels=c("bed", "good"))
group=levels(factor(data$group))
comp=combn(group, 2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

min(data$ST2)
# 开始画图
colnames(data)
score_name = colnames(data)

for(i in score_name){
  # i = "ST2"
  rt=data[,c(i, "group")]
  colnames(rt)=c("KP", "group")
  gg1=ggviolin(rt, x="group", y="KP", fill = "group", 
               xlab="group", ylab=i,
               legend.title="group",
               palette=c("#DB423E", "#008ECA"),
               add = "boxplot", add.params = list(fill="white"))+ 
    stat_compare_means(comparisons = my_comparisons,method = "wilcox.test")
  #stat_compare_means(comparisons = my_comparisons,symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif")
  print(gg1)
  pdf(file=paste0("01.ROC/02.group.result/",i, ".pdf"), width=6, height=5)
  print(gg1)
  dev.off()
}
