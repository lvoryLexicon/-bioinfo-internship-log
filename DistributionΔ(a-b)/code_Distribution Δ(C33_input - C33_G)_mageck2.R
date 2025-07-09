# 加载库
library(ggplot2)

# 第一步：读取数据
data <- read.table("~/data/7.8/mageck2_test.normalized.txt", header = TRUE, sep = "\t")

# 第二步：检查列名，确保正确
# tail(names(data), 2) 预期输出应是 "input_C33_73633" 和 "G_C33_73633"
target_cols <- tail(names(data), 2)

# 第三步：计算差值
diff_values <- data[[target_cols[1]]] - data[[target_cols[2]]]

# 第四步：处理极值（两种方法）

## 方法一：剔除上下1%分位数之外的数据
lower_bound <- quantile(diff_values, 0.01)
upper_bound <- quantile(diff_values, 0.99)
filtered_diff <- diff_values[diff_values >= lower_bound & diff_values <= upper_bound]

## 方法二：剔除超出3个标准差之外的数据
# mean_diff <- mean(diff_values)
# sd_diff <- sd(diff_values)
# filtered_diff <- diff_values[abs(diff_values - mean_diff) <= 3 * sd_diff]

# 第五步：画密度图，并叠加标准正态分布曲线
fig <- ggplot(data.frame(diff = filtered_diff), aes(x = diff)) +
  geom_density(fill = "skyblue", alpha = 0.6, color = "black")
fig <- fig + stat_function(fun = dnorm,
                args = list(mean = mean(filtered_diff), sd = sd(filtered_diff)),
                color = "red", linetype = "dashed", size = 1)
fig <- fig + labs(title = "差值密度图（极值处理后）",
       x = "input_C33_73633 - G_C33_73633",
       y = "Density") +   theme_minimal()


print(fig)
pdf("fig.pdf")
print(fig)
dev.off()
