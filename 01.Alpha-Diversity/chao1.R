library(openxlsx)
library(vegan)
library(ggplot2)
library(ggpubr)

rm(list = ls())
setwd("C:\\Users\\dengwei\\Google Drive\\HKU\\01.Project\\03.20200307_NAFLD\\scale_reads_to_35M\\microbitas\\chao_1")
dt <- read.xlsx("merged_reads_genus.xlsx", sheet = 1)
rownames(dt) <- dt$clade_name
dt <- dt[,-1]
dt_human <- dt[,1:24]
dt_mice <- dt[,25:48]

chao1_human <- data.frame(estimateR(t(dt_human))[2,])
chao1_mice <- data.frame(estimateR(t(dt_mice))[2,])

colnames(chao1_human) <- c("chao1")
colnames(chao1_mice) <- c("chao1")
chao1_human$names <- rownames(chao1_human)
chao1_mice$names <- rownames(chao1_mice)

group <- read.table("group.txt", header = T)

box_human <- merge(chao1_human, group, by="names")
box_mice <- merge(chao1_mice, group, by="names")

write.csv(box, file = "alpha_chao1.csv")

p_human<-ggplot(data=box_human, aes(x=type,y=chao1))+
  geom_boxplot(aes(fill=type))+
  labs(x="",y="Chao 1 richness", main = "")+
  theme_classic()+
  guides(fill=FALSE)+
  scale_x_discrete(labels=c("Healthy", "Lean NAFLD", "Obese NAFLD"))+
  theme(axis.text.x = element_text(size = 18, vjust = 0.5, hjust = 0.5, angle = 45))+
  theme(axis.text.y = element_text(size = 20))+
  theme(axis.title.y = element_text(size = 20, vjust = 4))+
  theme(plot.margin=unit(rep(2,4),'lines'))

p_mice<-ggplot(data=box_mice, aes(x=type,y=chao1))+
  geom_boxplot(aes(fill=type))+
  labs(x="",y="Chao 1 richness", main = "")+
  theme_classic()+
  guides(fill=FALSE)+
  scale_x_discrete(labels=c("FMT-Healthy", "FMT-Lean NAFLD", "FMT-Obese NAFLD"))+
  theme(axis.text.x = element_text(size = 18, vjust = 0.5, hjust = 0.5, angle = 45))+
  theme(axis.text.y = element_text(size = 20))+
  theme(axis.title.y = element_text(size = 20, vjust = 4))+
  theme(plot.margin=unit(rep(2,4),'lines'))

#p <- ggarrange(p_human, p_mice, 
#               labels = c("A", "B"),
#               ncol = 1, nrow = 2)
ggsave(p_human,filename = "chao1_human.pdf",width = 6.5,height = 8)
ggsave(p_mice,filename = "chao1_mice.pdf",width = 6.5,height = 8)

# multiple correction test
# 1) t test, adjust by BH
pairwise.t.test(box$chao1, box$type, p.adjust="BH", pool.sd=T)
# 2) anova test, adjust by TukeyHSD
aov_fit <- aov(chao1~type, data=box)
TukeyHSD(aov_fit, conf.level = .95)
# 3) mann-whitney test, adjust by BH
p_H_hl <- wilcox.test(box[box$type=="human_healthy", 1], box[box$type=="human_lean_NAFLD", 1], exact = FALSE)$p.value
p_H_ho <- wilcox.test(box[box$type=="human_healthy", 1], box[box$type=="human_obese_NAFLD", 1], exact = FALSE)$p.value
p_H_lo <- wilcox.test(box[box$type=="human_lean_NAFLD", 1], box[box$type=="human_obese_NAFLD", 1], exact = FALSE)$p.value
p_M_hl <- wilcox.test(box[box$type=="mice_healthy_FMT", 1], box[box$type=="mice_lean_FMT", 1], exact = FALSE)$p.value
p_M_ho <- wilcox.test(box[box$type=="mice_healthy_FMT", 1], box[box$type=="mice_obese_FMT", 1], exact = FALSE)$p.value
p_M_lo <- wilcox.test(box[box$type=="mice_obese_FMT", 1], box[box$type=="mice_lean_FMT", 1], exact = FALSE)$p.value
p_HM <- wilcox.test(box[box$species=="Human", 1], box[box$species=="Mice", 1], exact = FALSE)$p.value


p_value <- c(p_H_hl, p_H_ho, p_H_lo, p_M_hl, p_M_ho, p_M_lo, p_HM)
q_value <- p.adjust(p_value, method = "BH")

dt_value <- data.frame(p_value, q_value)
rownames(dt_value) <- c("H_hl", "H_ho", "H_lo", "M_hl", "M_ho", "M_lo", "HM")
write.table(dt_value, 'pvalue.txt', row.names = TRUE, sep = '\t', quote = FALSE, na = '')


