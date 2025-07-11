
#相关R包下载：
#install.packages("tidyverse")
#install.packages("ggsankey")
#install.packages("ggplot2")
#install.packages("cols4all")
#install.packages("cowplot")

#相关R包载入：
library(tidyverse)
library(ggsankey)
library(ggplot2)
library(cols4all)
library(cowplot)

#本地数据载入：
##富集气泡图数据：
kegg <- read.csv("C:/Users/SYZ/Documents/data/7.11/Sankey_Bubble_Combined/GO_enrichment_results.csv", header = TRUE)
head(kegg)

##桑基图数据：
sankey <- read.csv("C:/Users/SYZ/Documents/data/7.11/Sankey_Bubble_Combined/GO_enrichment_results.csv",header = T)
head(sankey)

#数据处理：
kegg2 <- kegg[22:1,] #先调转数据框方向
kegg2 <- kegg2 %>%
  mutate(ymax = cumsum(Count)) %>% #ymax为Width列的累加和
  mutate(ymin = ymax -Count) %>%
  mutate(label = (ymin + ymax)/2) #取xmin和xmax的中心位置作为标签位置
head(kegg2)

p2 <- ggplot() +
  geom_point(data = kegg2,
             aes(x = -log10(pvalue),
                 y = label, #替换y轴列为数值
                 size = Count,
                 color = RichFactor)) +
  theme_bw() +
  labs(x = "-log10(pvalue)",
       y = "")
p2

#自定义主题与配色修改：
mytheme <- theme(axis.title = element_text(size = 13),
                 axis.text = element_text(size = 11),
                 axis.text.y = element_blank(),
                 axis.ticks.y = element_blank(),
                 legend.title = element_text(size = 13),
                 legend.text = element_text(size = 11))
p3 <- ggplot() +
  geom_point(data = kegg2,
             aes(x = -log10(pvalue),
                 y = label,
                 size = Count,
                 color = RichFactor)) +
  scale_size_continuous(range=c(2,8)) + #调整气泡大小范围(默认尺寸部分过小)

  scale_colour_distiller(palette = "Reds", direction = 1) + #更改配色
  labs(x = "-log10(pvalue)",
       y = "") +
  theme_bw() +
  mytheme
p3

sankey <- sankey[22:1,]
#将数据转换为绘图所需格式：
df <- sankey %>%
  make_long(Description, geneID)
head(df)

#指定绘图顺序（转换为因子）：
df$node <- factor(df$node,levels = c(sankey$Description %>% unique() %>% rev(),
                                     sankey$geneID %>% unique()%>% rev()
                                     ))
#自定义配色：
c4a_gui()
mycol <- c4a('rainbow_wh_rd',41)
#绘图：
p4 <- ggplot(df, aes(
  x = x,
  next_x = next_x,
  node = node,
  next_node = next_node,
  fill = node   # 这里用 node 做 fill
)) +
  geom_sankey(flow.alpha = 0.5,
              flow.fill = 'grey',
              flow.color = 'grey80', #条带描边色
              node.fill = mycol, #节点填充色
              smooth = 8,
              width = 0.08) +
  scale_fill_manual(values = mycol)

p4 <- ggplot(df, aes(
  x = x,
  next_x = next_x,
  node = node,
  next_node = next_node,
  fill = node,
  label = node
)) +
  geom_sankey(
    flow.alpha = 0.5,
    flow.fill = "grey",
    flow.color = "grey80",
    smooth = 8,
    width = 0.08
  ) +
  geom_sankey_text(size = 3.2, color = "black") +
  scale_fill_manual(values = mycol) +
  theme_void() +
  theme(legend.position = "none")

p4

#通过调整页边距在桑基图右侧先留出空白：
p5 <- p4 + theme(plot.margin = unit(c(0,10,0,0),units="cm"))
p5
#拼图(在p5的空白位置中插入p3)：
ggdraw() + draw_plot(p5) + draw_plot(p3, scale = 0.5, x = 0.62, y=-0.21, width=0.48, height=1.3)
