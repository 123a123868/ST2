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
class(data.all$ST2)

# 设置基因 和分类 
ROC.data = data.all
colnames(ROC.data)
class.name = "prognosis"
gene = c("ST2")
gene = c("NT.proBNP")
gene = c("TnT")
library(dplyr)
#ROC.data = dplyr::filter(ROC.data, !is.na(NT.proBNP))
ROC.data = ROC.data[!ROC.data$TnT == "na",]

expr.gene = as.numeric(ROC.data[,gene])


class = as.character(ROC.data[,class.name])

# ROC 分析
# 创建ROC对象
roc_obj = roc(class,expr.gene,smooth=T)
# 绘制ROC曲线,计算AUC及其置信区间
plot(roc_obj,legacy.axes = TRUE,thresholds="best", # 基于约登指数选择roc曲线最佳阈值点
     print.thres="best")
auc_value <- auc(roc_obj) # 0.8786
auc_ci <- ci.auc(roc_obj) # 0.8575-0.8997
auc_value


#将修改过后的名字替换为roc1对象中基因的名称
ggroc(roc_obj,legacy.axes = T,color="#B2533E",linewidth=0.8)+
  annotate(geom = "segment", x = 0, y = 0, xend =1, yend = 1,
           linetype="dashed",color="#186F65")+
  annotate("text",x=0.8,y=0.3,label="AUC = 0.51")+
  labs(x="1-Specificity",y="Sensitivity")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 20))
ggsave("01.ROC/03ROC/TnT.ROC.pdf",height = 6,width = 6,family="Times")
