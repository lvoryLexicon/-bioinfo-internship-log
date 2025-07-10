#相关R包载入：
library(ggplot2)
library(ggrepel)
library(dplyr)

#测试数据载入：
df <- read.delim("C:/Users/SYZ/Documents/data/7.10/Volcano_Enrichment_BarPlots/readcount_differential.xls", header = TRUE)
head(df) #火山图数据

##用来简单看或是可视化一下数据分布
#summary(df$log2FC)
#hist(df$log2FC, breaks = 50, main = "log2FC分布", col = "skyblue")

#设定阈值：
Pvalue = 0.05
log2FC = 1

#根据阈值添加上下调分组标签：
df$group <- case_when(
  df$log2FC > log2FC & df$Pvalue < Pvalue ~ 'up_group',
  df$log2FC < -log2FC & df$Pvalue < Pvalue ~ 'down_group',
  TRUE ~ 'not_significant'
)
head(df)

#转换为因子,指定绘图顺序；
df$'-log10(Pvalue)' <- -log10(df$Pvalue) #新增-log10p列
df$group <- factor(df$group, levels = c('up_group','down_group','not_significant'))

#自定义配色：
mycol <- c("#8400ff","#dab3ff","#b9b9b9")

#自定义主题：
mytheme <- theme_bw() + 
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 13),
        legend.text = element_text(size = 13),
        panel.grid = element_blank()) 

#ggplot2绘制火山图：
p <- ggplot() + 
  geom_point(data = df,
             aes(x = log2FC, y = -log10(Pvalue), color = group),
             size = 2) + #添加散点
  scale_colour_manual(name = '', values = alpha(mycol, 0.7)) + #自定义散点颜色
  scale_x_continuous(limits = c(-5, 5), 
                     breaks = seq(-5, 5, by = 2)) + #x轴限制
  scale_y_continuous(expand = expansion(add = c(2, 0)), 
                     limits = c(0, 40), 
                     breaks = seq(0, 40, by = 10)) + #y轴限制
  geom_hline(yintercept = c(-log10(Pvalue)),
             size = 0.7,color = "black",lty = "dashed") + #水平阈值线
  geom_vline(xintercept = c(-log2FC, log2FC),
             size = 0.7,color = "black",lty = "dashed") + #垂直阈值线
  mytheme
p

#上下调Top10 id标签获取：
top10_up <- filter(df, group == "up_group") %>%
  distinct(gene_name, .keep_all = T) %>% top_n(10, abs(-log10(Pvalue)))

top10_down <- filter(df, group == "down_group") %>%
  distinct(gene_name, .keep_all = T) %>% top_n(10, abs(-log10(Pvalue)))

#表格合并：
sig <- rbind(top10_up, top10_down)
head(sig)

#目标id标签添加：
p1 <- p +
  geom_text_repel(data = sig,
                  aes(x = log2FC, y = -log10(Pvalue),
                      label = gene_name, color = group),
                  size = 3)
p1

#标签优化对齐：
p2 <- p +
  
  ##上调目标id：
  geom_text_repel(
    data = top10_up,
    aes(x = log2FC, y = -log10(Pvalue),
        label = gene_name, color = group),
    seed = 233,
    size = 3,
    min.segment.length = 0, #始终为标签添加指引线段；若不想添加线段，则改为Inf
    force = 2, #重叠标签间的排斥力
    force_pull = 2, #标签和数据点间的吸引力
    box.padding = 0.1, #标签周边填充量，默认单位为行
    max.overlaps = Inf, ##排斥重叠过多标签，设置为Inf则可以保持始终显示所有标签
    segment.linetype = 3, #线段类型,1为实线,2-6为不同类型虚线
    segment.alpha = 0.5, #线段不透明度
    nudge_x = 3-top10_up$log2FC, #标签x轴起始位置
    direction = "y", #按y轴调整标签位置方向，若想水平对齐则为x
    hjust = 0 #0右对齐，1左对齐，0.5居中
  )+
  
  ##下调目标id：
  geom_text_repel(
    data = top10_down,
    aes(x = log2FC, y = -log10(Pvalue),
        label = gene_name, color = group),
    seed = 233,
    size = 3,
    min.segment.length = 0,
    force = 2, 
    force_pull = 2, 
    box.padding = 0.1,
    max.overlaps = Inf, 
    segment.linetype = 3, 
    segment.alpha = 0.5,
    nudge_x = -3-top10_down$log2FC, 
    direction = "y", 
    hjust = 1
  )
p2

#新增-log10Padj列，注意这里根据上下调分组赋予正负：
dt <- df  # 
dt <- dt %>% filter(!is.na(Padj) & Padj > 0)

dt$`-log10Padj` <- case_when(
  dt$group == "up_group" ~ -log10(dt$Padj),
  dt$group == "down_group" ~ log10(dt$Padj),
  TRUE ~ 0  # 或 NA_real_
)

##转换为因子，指定绘图顺序；
dt$description <- factor(dt$description, levels = rev(unique(dt$description)))

#常规版上下调富集条形图绘制：
# 只保留上下调各10个通路（共20个）
top_up <- dt %>% filter(`-log10Padj` < 0) %>% top_n(-4, `-log10Padj`)
top_down <- dt %>% filter(`-log10Padj` > 0) %>% top_n(3, `-log10Padj`)
dt <- bind_rows(top_up, top_down)

# 确保group为factor
dt$group <- factor(dt$group, levels = c("up_group", "down_group"))

# 拆分数据
up <- dt %>% filter(`-log10Padj` > 0)
down <- dt %>% filter(`-log10Padj` < 0)

# 自定义主题
mytheme1 <- theme(
  legend.position = 'none',
  axis.text.y = element_blank(),
  axis.ticks.y = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.border = element_blank(),
  axis.line.x = element_line(color = 'grey60', size = 1.1),
  axis.text = element_text(size = 12)
)

# 一次性绘图
p_final <- ggplot(data = dt, aes(x = `-log10Padj`, y = description, fill = group)) +
  geom_col() +
  geom_text(data = up, aes(x =- 2, y = description, label = description),
            size = 4, hjust = 1) +
  geom_text(data = down, aes(x =+ 2, y = description, label = description),
            size = 4, hjust = 0) +
  coord_cartesian(xlim = c(min(dt$`-log10Padj`) * 1.1, max(dt$`-log10Padj`) * 1.1)) +
  scale_x_continuous(breaks = seq(floor(min(dt$`-log10Padj`)), ceiling(max(dt$`-log10Padj`)), by = 50)) +
  labs(x = '-log10Padj', y = '', title = 'Descriptions dysgroup in NYHA class IV') +
  theme_bw() +
  mytheme1 +
  theme(plot.title = element_text(hjust = 0.5, size = 14)) +
  scale_fill_manual(values = c("up_group" = "#8400ff", "down_group" = "#dab3ff"))

p_final