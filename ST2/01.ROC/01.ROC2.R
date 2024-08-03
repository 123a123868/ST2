# 设置
rm(list = ls())
library(pROC)
library(ggsci)
library(ggplot2)
#做分类ROC曲线 可平滑

# 读入数据
# 再读进来看看什么情况
data.all = read.csv("01.ROC/01.filterData/filter.3.csv")

# 先区分预后好和预后差的两组
data_good = data.all[data.all$MACE.4 == "good",]
data_bed = data.all[data.all$MACE.4 == "bed",]





colnames(data.all)
fea = c("ST2","Gensini_score")


ROC.data = data.all
formula1<-as.formula(paste("prognosis~",paste(fea,collapse = "+")))
#用pROC得到roc曲线指标
roc1<-roc(formula=as.formula(formula1),data = ROC.data,
          smooth=T)#平滑ROC曲线与否


#把auc取出来放入列表，用来在基因名后面加上auc值
auc_list1<-c()
for (i in names(roc1)) {
  auc_list1[i]<-roc1[[i]]$auc
}
auc_list1<-as.data.frame(auc_list1)

#edit_name就是修改过后的名字
auc_list1$edit_name<-paste(rownames(auc_list1),'AUC:',round(auc_list1$auc_list,3))

#将修改过后的名字替换为roc1对象中基因的名称
names(roc1)<-auc_list1$edit_name


pROC::ggroc(roc1,legacy.axes=T,alpha=1,size=0.8)+
  theme_bw(base_size = 12)+
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color="darkgrey", linetype="dashed")+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+#去掉背景网格线
  theme(legend.position = c(0.8,0.2))+#图例的坐标位置
  scale_color_tron()+
  theme(legend.text = element_text(size = 10),#图例文字
        axis.text = element_text(size = 10))+#刻度文本
  theme(legend.title = element_blank())+
  guides(color=guide_legend(override.aes = list(size=2)))#图例中线条粗细

