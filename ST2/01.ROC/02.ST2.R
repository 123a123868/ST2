# 先读入数据
# 确定号预后分组
# 画一个小提琴图



# 设置
rm(list = ls())
library(pROC)
library(ggsci)
library(ggplot2)
#做分类ROC曲线 可平滑

# 读入数据
# 再读进来看看什么情况
library(readxl)
data.all = read_xlsx("01.ROC/01.filterData/filter.final.xlsx")
data.all = as.data.frame(data.all)
class(data.all$prognosis)
class(data.all$Killip_degree)
table(data.all$prognosis)
colnames(data.all)
data.long = data.all[,c("Killip_degree","prognosis")]
colnames(data.long)[2] = "group"


# ggpolt画图
library(ggplot2)
colnames(data.long)
# Multiple groups with error bars and jitter point
library(ggpubr)
table(data.long$group)
colnames(data.long)
y = 1.5
#data.long$Killip_degree = as.numeric(data.long$Killip_degree)
ggbarplot(data.long, x = 'group', y = 'Killip_degree', add = c('mean_se'),
          color = "group", fill = "group",palette =  c('#000000','#808080'),
          position = position_dodge(1)) +
  labs(x = '', y = 'Killip_degree') + 
  theme( legend.position = 'none' ,
         text=element_text(size=16)) +
  rotate_x_text(angle = 50)+
  coord_cartesian(ylim = c(0, y))+  # 设置y轴的坐标范围 
  stat_compare_means(method = 'wilcox.test', aes(group = group),label.y = y)
ggsave("01.ROC/02.group.result.2/geom_bar.Killip_degree.pdf",height = 5,width = 5,family="Times")

