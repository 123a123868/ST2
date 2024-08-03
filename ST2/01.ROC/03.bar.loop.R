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
data.all = read_xlsx("01.ROC/01.filterData/filter.final.loop.score.xlsx")
data.all = as.data.frame(data.all)
class(data.all$prognosis)
class(data.all$ST2)
table(data.all$prognosis)
colnames(data.all)
# 转置为数字格式
data.all[,2:23] = apply(data.all[,2:23], 2, as.numeric)


kp = colnames(data.all)
kp


i = 6
ylim = 13
kp.name = kp[i]
data.long = data.all[,c(kp.name,"prognosis")]
colnames(data.long) = c("kp.name","group")
# ggpolt画图
library(ggplot2)
colnames(data.long)
# Multiple groups with error bars and jitter point
library(ggpubr)
table(data.long$group)
colnames(data.long)
max(data.long$kp.name)

ggbarplot(data.long, x = 'group', y = 'kp.name', add = c('mean_se'),
          color = "group", fill = "group",palette =  c('#808080','#000000'),
          position = position_dodge(1)) +
  labs(y = kp.name) + 
  theme( legend.position = 'none' ,
         text=element_text(size=20)) +
  rotate_x_text(angle = 50)+
  coord_cartesian(ylim = c(0, ylim))+  # 设置y轴的坐标范围 
  stat_compare_means(method = 'wilcox.test', aes(group = group),label.y = ylim)
ggsave(paste("01.ROC/02.group.result.2//geom_bar.", kp.name, ".pdf"),height = 5,width = 5,family="Times")


















for (i in 10:length(kp)) {
  # i = 3
  kp.name = kp[i]
  data.long = data.all[,c(kp.name,"prognosis")]
  colnames(data.long) = c("kp.name","group")
  # ggpolt画图
  library(ggplot2)
  colnames(data.long)
  # Multiple groups with error bars and jitter point
  library(ggpubr)
  table(data.long$group)
  colnames(data.long)
  max(data.long$kp.name)
  ylim = 15
  ggbarplot(data.long, x = 'group', y = 'kp.name', add = c('mean_se'),
            color = "group", fill = "group",palette =  c('#808080','#000000'),
            position = position_dodge(1)) +
    labs(y = kp.name) + 
    theme( legend.position = 'none' ,
           text=element_text(size=20)) +
    rotate_x_text(angle = 50)+
    coord_cartesian(ylim = c(0, ylim))+  # 设置y轴的坐标范围 
    stat_compare_means(method = 'wilcox.test', aes(group = group),label.y = ylim)
  ggsave(paste("01.ROC/02.group.result.2//geom_bar.", kp.name, ".pdf"),height = 5,width = 5,family="Times")

  
  
  }


ggsave(paste("01.ROC/02.group.result/geom_bar.LVEF." ,".pdf"),height = 5,width = 5)











i = 23
kp.name = kp[i]
data.long = data.all[,c(kp.name,"prognosis")]
colnames(data.long) = c("kp.name","group")
# ggpolt画图
library(ggplot2)
colnames(data.long)
# Multiple groups with error bars and jitter point
library(ggpubr)
table(data.long$group)
colnames(data.long)
max(data.long$kp.name)
data.long$kp.name = as.numeric(data.long$kp.name)

ylim = 4
ggbarplot(data.long, x = 'group', y = 'kp.name', add = c('mean_se'),
          color = "group", fill = "group",palette =  c('#808080','#000000'),
          position = position_dodge(1)) +
  labs(y = kp.name) + 
  theme( legend.position = 'none' ) +
  rotate_x_text(angle = 50)+
  coord_cartesian(ylim = c(0, ylim))+  # 设置y轴的坐标范围 
  stat_compare_means(method = 'wilcox.test', aes(group = group),label.y = ylim)
ggsave(paste("01.ROC/02.group.result/geom_bar.", kp.name, ".pdf"),height = 5,width = 5,family="Times")
max(data.long$kp.name)
























