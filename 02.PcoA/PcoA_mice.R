library(openxlsx)
library(ade4)
library(ggplot2)
library(RColorBrewer)
library(plyr)
library(vegan)#用于计算距离
library(ggsci)

setwd("C:\\Users\\dengwei\\Google Drive\\HKU\\01.Project\\03.20200307_NAFLD\\scale_reads_to_35M\\microbitas\\PcoA\\mice")
data <- read.xlsx("merged_mice_genus.xlsx", sheet = 1)
group <- read.table('group.txt', sep = '\t', header = T, stringsAsFactors = FALSE)
row.names(data) <- data[,1]
data <- data[,-1]

distance <- vegdist(t(data), method = 'bray')
pcoa <- cmdscale(distance, k = (nrow(t(data)) - 1), eig = TRUE)

# save the distance info
write.csv(as.matrix(distance), 'distance.csv', quote = F)
write.csv(as.matrix(pcoa$points), 'pcoa_value.csv', quote = F)
#-------------------------------------------------------------------
# plot using 'ordiplot', a in-house function in vegan
ordiplot(scores(pcoa)[ ,c(1, 2)], type = 't')
# check info
summary(pcoa)

# check the coordinate value
pcoa$eig
point <- data.frame(pcoa$point)

# save the coordinate info
write.csv(point, 'pcoa.sample.csv')
#-------------------------------------------------------------------

# compute first two coordinate value
pcoa_eig <- (pcoa$eig)[1:2] / sum(pcoa$eig)

# extract first two coordinate value
sample_site <- data.frame({pcoa$point})[1:2]
sample_site$names <- rownames(sample_site)
names(sample_site)[1:2] <- c('PCoA1', 'PCoA2')

#为样本点坐标添加分组信息
sample_site <- merge(sample_site, group, by = 'names', all.x = FALSE)
#可选输出，例如输出为 csv 格式
write.csv(sample_site, 'sample_site.csv', quote = F)

#sample_site$site <- factor(sample_site$species, levels = c('Human', 'Mice'))
sample_site$type <- factor(sample_site$type, levels = c('healthy_FMT', 'lean_FMT', 'obese_FMT'))


pcoa_plot <- ggplot(sample_site, aes(PCoA1, PCoA2,colour=type, shape = type)) +
  theme_classic()+#去掉背景框
  geom_vline(xintercept = 0, color = 'gray', size = 0.4) + 
  geom_hline(yintercept = 0, color = 'gray', size = 0.4) +
  geom_point(size = 3)+  #可在这里修改点的透明度、大小
  theme(panel.grid = element_line(color = 'gray', linetype = 2, size = 0.1), 
        panel.background = element_rect(color = 'black', fill = 'transparent'), 
        legend.title=element_blank())+
  labs(x = paste('PCoA1: ', round(100 * pcoa_eig[1], 2), '%'), y = paste('PCoA2: ', round(100 * pcoa_eig[2], 2), '%'))+
  theme(axis.text = element_text(size=14), axis.title = element_text(size=16)) +
  theme(plot.title = element_text(hjust = 0.5))+
  scale_color_npg()

pcoa_plot

ggsave('PCoA_2.png', pcoa_plot, width = 7, height = 5)
ggsave('PCoA_2.pdf', pcoa_plot, width = 7, height = 5)

