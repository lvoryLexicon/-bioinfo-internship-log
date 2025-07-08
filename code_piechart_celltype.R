library(dplyr)
library(ggplot2)
library(RColorBrewer)
 
# 统计每类 celltype 的数量和比例，并创建新标签
label_df <- celltype@meta.data %>%
    +     group_by(celltype) %>%
    +     summarise(count = n()) %>%
    +     mutate(percent = round(100 * count / sum(count), 1),
                 +            new_label = paste0(celltype, " (", count, ", ", percent, "%)"))
 
# 把新标签合并回原始 metadata
metadata <- celltype@meta.data %>%
      +     left_join(label_df, by = "celltype")
 
# 统计新标签后的数量（用于画图）
plot_df <- metadata %>%
        +     group_by(new_label) %>%
        +     summarise(count = n())

# 配色
n_colors <- nrow(plot_df)
palette <- colorRampPalette(brewer.pal(8, "Set3"))(n_colors)
  
# 画饼图
ggplot(plot_df, aes(x = "", y = count, fill = new_label)) +
          +     geom_bar(stat = "identity", width = 1, color = "white") +
          +     coord_polar("y") +
          +     scale_fill_manual(values = palette) +
          +     theme_void() +
          +     labs(title = "各细胞类型占比饼图", fill = "细胞类型") +
          +     theme(legend.text = element_text(size = 9),
                      +           plot.title = element_text(hjust = 0.5, size = 14))