library(ggplot2)
library(dplyr)

# 载入数据
dt <- read.csv("C:/Users/SYZ/Documents/data/7.11/Pathway_Gene_Display/GO_enrichment_results.csv", 
               header = TRUE, check.names = FALSE)

# 按pvalue升序排序，仅保留前10个最显著的条目
dt_top10 <- dt %>% 
  arrange(pvalue) %>% 
  head(10)

# 创建新的Description因子（按pvalue从大到小排序）
dt_top10$Description <- factor(
  dt_top10$Description,
  levels = rev(dt_top10$Description)
)

# 基础柱状图
p <- ggplot(dt_top10) +
  geom_bar(
    aes(x = -log10(pvalue), y = Description),
    stat = "identity",
    width = 0.6,
    fill = "#4E79A7"
  ) +
  theme_classic() +
  labs(
    x = "-log10(p-value)", 
    y = "GO Terms",
    title = "Top 10 Significant GO Terms"
  ) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

# 添加GO terms标签（使用筛选后的数据）
p2 <- p + 
  geom_text(
    aes(x = 0.1,
        y = Description,
        label = Description),
    size = 4.5,
    hjust = 0,
    color = "black"
  )

# 添加Marker基因信息（关键修正：使用dt_top10而非dt）
p3 <- p2 +
  geom_text(
    aes(x = 0.1, 
        y = Description, 
        label = geneID),  # 确保列名为geneID
    size = 4,
    fontface = 'italic',
    hjust = 0,
    vjust = 2.3,  # 垂直偏移
    color = "grey30"  # 使用深灰色提高可读性
  )

# 显示最终图形
print(p3)