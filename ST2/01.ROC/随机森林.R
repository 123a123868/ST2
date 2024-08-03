library(randomForest)
library(ggplot2)
library(dplyr)

# 假设您有一个二分类问题的数据集，其中包含特征和目标变量
# 以下是一个示例数据集
data <- data.frame(
  feature1 = c(1, 2, 3, 4, 5),
  feature2 = c(2, 3, 4, 5, 6),
  feature3 = c(3, 4, 5, 6, 7),
  target = factor(c("A", "A", "B", "B", "A"))
)

# 拆分数据集为特征和目标变量
features <- data %>% select(-target)
target <- data$target

# 训练随机森林模型
rf_model <- randomForest(features, target)

# 提取变量重要性
importance <- importance(rf_model)

# 创建数据框
importance_df <- data.frame(
  variable = row.names(importance),
  importance = importance[, "MeanDecreaseGini"]
)

# 按照重要性排序
importance_df <- importance_df %>%
  arrange(desc(importance))

# 绘制深林图
ggplot(importance_df, aes(x = reorder(variable, importance), y = importance, fill = variable)) +
  geom_bar(stat = "identity", alpha = 0.75) +
  xlab("Variable") +
  ylab("Importance") +
  ggtitle("Variable Importance Plot") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        axis.line = element_line(color = "black", size = 0.5),
        axis.ticks = element_line(color = "black", size = 0.5)) +
  scale_fill_brewer(palette = "Set3")