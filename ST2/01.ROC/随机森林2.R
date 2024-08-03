# ST2 和其他标准，做随机深林分析
# R语言ggthemes包 theme_few函数使用说明


rm(list=ls())
options(stringsAsFactors = F) 
library(randomForest)
library(ggthemes)
library(ggplot2)

# 构造数据集
set.seed(123) # 为了结果可重现
num_samples <- 1000
num_genes <- 50
gene_data <- matrix(rnorm(num_samples * num_genes), nrow = num_samples, ncol = num_genes)
colnames(gene_data) <- paste0("Gene_", 1:num_genes)
heart_disease <- factor(sample(c("Yes", "No"), num_samples, replace = TRUE))

# 准备数据框
data <- data.frame(gene_data, Disease = heart_disease)

# 建立随机森林模型
rf_model <- randomForest(Disease ~ ., data = data, ntree = 500, importance = TRUE)



library(gridExtra)
# 获取特征重要性并排序
importance_scores <- importance(rf_model)
# 基于基尼指数的特征重要性
importance_gini <- importance(rf_model, type = 2)
# 基于准确度的特征重要性
importance_accuracy <- importance(rf_model, type = 1)
# 基于默认的特征重要性
feature_importance <- data.frame(Gene = rownames(importance_scores), Importance = importance_scores[, "MeanDecreaseGini"])
ordered_features <- feature_importance[order(-feature_importance$Importance), ]
# 基于基尼指数的特征重要性
feature_importance_gini <- data.frame(Gene = rownames(importance_gini), Importance = importance_gini[, "MeanDecreaseGini"])
ordered_features_gini <- feature_importance_gini[order(-feature_importance_gini$Importance), ]
# 基于准确度的特征重要性
feature_importance_accuracy <- data.frame(Gene = rownames(importance_accuracy), Importance = importance_accuracy[, "MeanDecreaseAccuracy"])
ordered_features_accuracy <- feature_importance_accuracy[order(-feature_importance_accuracy$Importance), ]

# 设置渐变颜色范围
color_range <- c("#D8BFD8", "#8B008B")

# 绘制基于基尼指数的特征重要性图
plot_gini <- ggplot(ordered_features_gini, aes(x = reorder(Gene, Importance), y = Importance, fill = Importance, label = round(Importance, 2))) +
  geom_bar(stat = "identity") +
  geom_text(size = 3, position = position_stack(vjust = 0.5), color = "white") +
  scale_fill_gradient(low = color_range[1], high = color_range[2]) +
  theme_few() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),
    plot.caption = element_text(size = 12),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 15)
  ) +
  labs(title = "Feature Importance (Gini Index)", x = "Gene", y = "Importance") +
  coord_flip()

# 绘制基于准确度的特征重要性图
plot_accuracy <- ggplot(ordered_features_accuracy, aes(x = reorder(Gene, Importance), y = Importance, fill = Importance, label = round(Importance, 2))) +
  geom_bar(stat = "identity") +
  geom_text(size = 3, position = position_stack(vjust = 0.5), color = "white") +
  scale_fill_gradient(low = color_range[1], high = color_range[2]) +
  theme_few() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),
    plot.caption = element_text(size = 12),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 15)
  ) +
  labs(title = "Feature Importance (Accuracy)", x = "Gene", y = "Importance") +
  coord_flip()

# 绘制模仿的图
plot_immitation <- ggplot(ordered_features, aes(x = reorder(Gene, Importance), y = Importance, fill = Importance, label = round(Importance, 2))) +
  geom_bar(stat = "identity") +
  geom_text(size = 3, position = position_stack(vjust = 0.5), color = "white") +
  scale_fill_gradient(low = color_range[1], high = color_range[2]) +
  theme_few() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),
    plot.caption = element_text(size = 12),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 15)
  ) +
  labs(title = "Feature Importance", x = "Gene", y = "Importance") +
  coord_flip()

# 在一张页面上并排展示两个图
grid.arrange(plot_gini, plot_accuracy,plot_immitation, ncol = 3)