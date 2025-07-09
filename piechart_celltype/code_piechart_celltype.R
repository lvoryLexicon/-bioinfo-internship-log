# 加载库
library(dplyr)
library(ggplot2)
library(RColorBrewer)

# Step 1: 提取 metadata
meta_data <- celltype@meta.data

# Step 2: 每类 celltype 的数量和百分比
celltype_stats <- meta_data %>%
  group_by(celltype) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(
    percent = round(100 * count / sum(count), 1),
    new_label = paste0(celltype, " (", count, ", ", percent, "%)")
  )

# Step 3: 把新标签合并到原始 metadata
meta_with_label <- left_join(meta_data, celltype_stats, by = "celltype")

# Step 4: 用 new_label 分组做统计，用于绘图
plot_data <- meta_with_label %>%
  group_by(new_label) %>%
  summarise(count = n(), .groups = "drop")

# Step 5: 设置配色（根据类别数量自动生成）
num_labels <- nrow(plot_data)
color_palette <- colorRampPalette(brewer.pal(8, "Set3"))(num_labels)

# Step 6: 绘图 - 饼图展示每类细胞的比例
fig <- ggplot(plot_data, aes(x = "", y = count, fill = new_label)) +
  geom_bar(stat = "identity", width = 1, color = "white") 
fig <- fig + coord_polar("y")
fig <- fig + scale_fill_manual(values = color_palette)
fig <- fig + theme_void()
fig <- fig +labs(title = "各细胞类型占比饼图", fill = "细胞类型")
fig <- fig + theme(
    legend.text = element_text(size = 9),
    plot.title = element_text(hjust = 0.5, size = 14)
  )
#显示和保存图
print(fig)
pdf("fig.pdf")
print(fig)
dev.off()
