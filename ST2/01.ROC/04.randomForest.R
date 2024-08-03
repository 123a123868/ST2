# ST2 和其他标准，做随机深林分析
# R语言ggthemes包 theme_few函数使用说明



# 设置
rm(list = ls())
library(randomForest)
library(ggthemes)
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
data.all$prognosis = factor(data.all$prognosis,levels = c("good","bed"))

# 准备数据框
data = data.all

# 不能纳入缺失值
# 去掉含有缺失值的行
data <- na.omit(data)


# 建立随机森林模型
rf_model <- randomForest(prognosis ~ ., data = data, ntree = 500, importance = TRUE)



library(gridExtra)
# 获取特征重要性并排序
importance_scores <- importance(rf_model)
# 基于基尼指数的特征重要性
importance_gini <- importance(rf_model, type = 2)
# 基于准确度的特征重要性
importance_accuracy <- importance(rf_model, type = 1)
# 基于默认的特征重要性
feature_importance <- data.frame(Feature = rownames(importance_scores), Importance = importance_scores[, "MeanDecreaseGini"])
ordered_features <- feature_importance[order(-feature_importance$Importance), ]
# 基于基尼指数的特征重要性
feature_importance_gini <- data.frame(Feature = rownames(importance_gini), Importance = importance_gini[, "MeanDecreaseGini"])
ordered_features_gini <- feature_importance_gini[order(-feature_importance_gini$Importance), ]
# 基于准确度的特征重要性
feature_importance_accuracy <- data.frame(Feature = rownames(importance_accuracy), Importance = importance_accuracy[, "MeanDecreaseAccuracy"])
ordered_features_accuracy <- feature_importance_accuracy[order(-feature_importance_accuracy$Importance), ]

# 设置渐变颜色范围
color_range <- c("#D8BFD8", "#8B008B")

# 绘制基于基尼指数的特征重要性图
plot_gini <- ggplot(ordered_features_gini, aes(x = reorder(Feature, Importance), y = Importance, fill = Importance, label = round(Importance, 2))) +
  geom_bar(stat = "identity") +
  geom_text(size = 3, position = position_stack(vjust = 0.5), color = "white") +
  scale_fill_gradient(low = color_range[1], high = color_range[2]) +
  theme_few() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),
    plot.caption = element_text(size = 12),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 15),
    text=element_text(size=16,  family="Times")
  ) +
  labs(title = "Feature Importance (Gini Index)", x = "Feature", y = "Importance") +
  coord_flip()

# 绘制基于准确度的特征重要性图
plot_accuracy <- ggplot(ordered_features_accuracy, aes(x = reorder(Feature, Importance), y = Importance, fill = Importance, label = round(Importance, 2))) +
  geom_bar(stat = "identity") +
  geom_text(size = 3, position = position_stack(vjust = 0.5), color = "white") +
  scale_fill_gradient(low = color_range[1], high = color_range[2]) +
  theme_few() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),
    plot.caption = element_text(size = 12),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 15),
    text=element_text(size=16,  family="Times")
  ) +
  labs(title = "Feature Importance (Accuracy)", x = "Feature", y = "Importance") +
  coord_flip()

# 绘制模仿的图
plot_immitation <- ggplot(ordered_features, aes(x = reorder(Feature, Importance), y = Importance, fill = Importance, label = round(Importance, 2))) +
  geom_bar(stat = "identity") +
  geom_text(size = 3, position = position_stack(vjust = 0.5), color = "white") +
  scale_fill_gradient(low = color_range[1], high = color_range[2]) +
  theme_few() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),
    plot.caption = element_text(size = 12),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 15),
    text = element_text(size = 20)
  ) +
  labs(title = "Feature Importance", x = "Feature", y = "Importance") +
  coord_flip()

# 在一张页面上并排展示两个图
grid.arrange(plot_gini, plot_accuracy, ncol = 2)
pdf("01.ROC/03ROC/randomForest.pdf",height = 6,width = 14)
grid.arrange(plot_gini, plot_accuracy, ncol = 2)
dev.off()



