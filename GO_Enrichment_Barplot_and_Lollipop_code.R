## 加载必要包
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(dplyr)
library(patchwork)

## 读取表达数据
gene_expr_df <- read.delim("C:/Users/SYZ/Documents/data/7.9/readcount_differential.xls", header = TRUE)

## 去除 Ensembl ID 版本号
gene_expr_df$geneID_clean <- sub("\\..*", "", gene_expr_df$geneID)

## 筛选差异表达基因
deg_df <- gene_expr_df %>%
  filter(Padj < 0.05, abs(log2FC) > 1)

## Ensembl → Entrez ID 映射
gene_df <- bitr(deg_df$geneID_clean, fromType = "ENSEMBL", 
                toType = "ENTREZID", OrgDb = org.Hs.eg.db)

entrez_ids <- gene_df$ENTREZID

## GO 富集分析（Biological Process）
go_result <- enrichGO(gene         = entrez_ids,
                      OrgDb        = org.Hs.eg.db,
                      ont          = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.2,
                      readable      = TRUE)

## 保存结果
write.csv(as.data.frame(go_result), "GO_enrichment_results.csv", row.names = FALSE)

## 提取前 20 项
go_top20_df <- as.data.frame(go_result) %>%
  arrange(p.adjust) %>%
  slice_head(n = 20)

## 设置因子顺序
go_top20_df$Description <- factor(go_top20_df$Description, levels = rev(go_top20_df$Description))

## 图1：灰色细柱图 + 气泡（Count 显著性）
bar_bubble_plot <- ggplot(go_top20_df, aes(x = -log10(p.adjust), y = Description)) +
  geom_col(width = 0.1, fill = 'grey') +
  geom_point(aes(size = 12, color = Count)) +
  scale_color_gradient(low = "#66c2a5", high = "#d73027") +
  theme_bw() +
  labs(title = "GO Enrichment - Bubble Overlay",
       x = "-log10(p.adjust)", y = NULL) +
  theme(axis.text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "none")

## 图2：条形图 + 文字标签
bar_label_plot <- ggplot(go_top20_df, aes(x = -log10(p.adjust), y = Description, fill = Count)) +
  geom_bar(stat = "identity", width = 0.5) +
  scale_fill_gradient(low = "#66c2a5", high = "#d73027") +
  geom_text(aes(x = 0.5, label = Description), size = 4, hjust = 0) +
  theme_bw() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, max(-log10(go_top20_df$p.adjust)) + 1)) +
  labs(title = "GO Enrichment - Lollipop Style",
       x = "-log10(p.adjust)", y = "")

## 合并展示（横向）
bar_bubble_plot + bar_label_plot
