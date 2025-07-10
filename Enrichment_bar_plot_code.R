#相关R包载入：
library(ggplot2)
library(RColorBrewer)

#载入本地差异富集结果表：
dt <- read.csv("C:/Users/SYZ/Documents/data/7.10/Volcano_Enrichment_BarPlots/GO_enrichment_results.csv", header = TRUE)

#查看列名信息：
colnames(dt)

#按照通路指定绘图顺序(转换为因子)：
dt <- dt[order(dt$pvalue), ][1:20, ]
dt$Description <- factor(dt$Description, levels = rev(unique(dt$Description)))

#先自定义主题：
mytheme <- theme(
  axis.title = element_text(size = 13),
  axis.text = element_text(size = 11),
  plot.title = element_text(size = 14,
                            hjust = 0.5,
                            face = "bold"),
  legend.title = element_text(size = 13),
  legend.text = element_text(size = 11),
  plot.margin = margin(t = 5.5,
                       r = 10,
                       l = 5.5,
                       b = 5.5)
)
#常规富集条形图绘图：
p <- ggplot(data = dt, aes(x = Count, y = Description, fill = -log10(pvalue))) +
  scale_fill_distiller(palette = "RdPu",direction = 1) + #更改配色
  geom_bar(stat = "identity", width = 0.8) + #绘制条形图
  labs(x = "Number of Gene", y = "", title = "KEGG enrichment barplot") + #修改/添加各标题
  theme_bw() + mytheme #主题更改
p

#先分别绘制:
#绘制富集条形图展示-log10p：
p2 <- ggplot(data = dt,aes(x = Description, y = -log10(pvalue))) +
  geom_bar(stat = "identity", width = 0.8, fill = 'pink') +
  scale_y_continuous(limits = c(0, 25)) +
  labs(x = "", y = "-log10(pvalue)", title = "KEGG enrichment barplot") +
  theme_classic() + mytheme +
  coord_flip()
p2

#绘制折线散点图展示count数目：
p3 <- ggplot(data = dt, aes(x = Description, y = Count, group = 1)) +
  geom_line(color = "grey", cex = 0.5) + #绘制折线
  geom_point(size = 2) + #绘制Count散点
  scale_y_continuous(limits = c(20, 60)) +
  labs(x = "", y = "Count", title = "KEGG enrichment barplot") +
  theme_classic() + mytheme +
  coord_flip()
p3

#将不同坐标系图表进行组合：
p4 <- ggplot() +
  geom_bar(data = dt, aes(x = Description, y = -log10(pvalue)),
           stat = "identity", width = 0.8, fill = 'pink') +
  scale_y_continuous(limits = c(0,25), #主y轴刻度限制
                     breaks = seq(0,25,5), #主y轴刻度线断点
                     expand = c(0,0),
                     sec.axis = sec_axis(~. *2.4, #创建与主轴相对的辅助轴，通过函数计算实现1对1的转换。如主轴范围为0-25，这里我们通过*2.4实现辅助轴范围为0-60
                                         name = 'Count', #辅助轴标题
                                         breaks = seq(0,60,10))) + #辅助轴刻度线断点
  geom_line(data = dt, aes(x = Description, y = Count/2.4, group = 1),
            color = "grey", cex = 0.5) +
  geom_point(data = dt, aes(x = Description, y = Count/2.4),
             size = 2) +
  labs(x = "", y = "-log10(pvalue)", title = "KEGG enrichment barplot") +
  theme_classic() +
  mytheme +
  coord_flip() #由于geom_line函数绘制折线是按照变量在x轴水平方向顺序进行连接，因此这里把x/y轴对应变量交换了一下，最后再通过坐标翻转回来
p4