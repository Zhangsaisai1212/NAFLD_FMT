
library(openxlsx)
library(vegan)
library(ggplot2)

setwd("C:\\Users\\dengwei\\Google Drive\\HKU\\01.Project\\03.20200307_NAFLD\\scale_reads_to_35M\\microbitas\\shannon")

rm(list = ls())
workbook <- "merged_genus.xlsx"
mydayaframe <- read.xlsx(workbook, 1)
rownames(mydayaframe) <- mydayaframe$clade_name
dt <- t(mydayaframe[,-1]) # remove the sample id

dt_human <- dt[1:24,]
dt_mice <- dt[25:48,]

#alpha_shannon <- diversity(dt, index = "shannon")
#write.csv(alpha_shannon, file = "alpha_shannon.csv")

alpha_shannon_human <- diversity(dt_human, index = "shannon")
alpha_shannon_mice <- diversity(dt_mice, index = "shannon")

group <- read.table("group.txt", header = T)

alph_human <- data.frame(names=names(alpha_shannon_human),value=alpha_shannon_human)
alph_mice <- data.frame(names=names(alpha_shannon_mice),value=alpha_shannon_mice)

human <- merge(alph_human, group, by="names")
mice <- merge(alph_mice, group, by="names")


## the following precedures were performed after editing the "alpha_shannon.csv" 
shannon <- read.csv("alpha_shannon.csv", header = T)

p<-ggplot(data=shannon, aes(x=Group,y=shannon))+
  geom_boxplot(aes(fill=Group))+
  labs(x="",y="Shannon index", main = "")+
  guides(fill=FALSE)+
  theme(axis.text.x = element_text(size = 15, vjust = 0.5, hjust = 0.5, angle = 45))+
  theme(axis.text.y = element_text(size = 15))+
  theme(axis.title.y = element_text(size = 20, vjust = 4))+
  theme(plot.margin=unit(rep(2,4),'lines'))
  
p

ggsave(p,filename = "shannon.pdf",width = 10,height = 9)


p_human<-ggplot(data=human, aes(x=type,y=value))+
  geom_boxplot(aes(fill=type))+
  labs(x="",y="Shannon index", main = "")+
  theme_classic()+
  guides(fill=FALSE)+
  scale_x_discrete(labels=c("Healthy", "Lean NAFLD", "Obese NAFLD"))+
  theme(axis.text.x = element_text(size = 18, vjust = 0.5, hjust = 0.5, angle = 45))+
  theme(axis.text.y = element_text(size = 20))+
  theme(axis.title.y = element_text(size = 20, vjust = 4))+
  theme(plot.margin=unit(rep(2,4),'lines'))

p_mice<-ggplot(data=mice, aes(x=type,y=value))+
  geom_boxplot(aes(fill=type))+
  labs(x="",y="Shannon index", main = "")+
  theme_classic()+
  guides(fill=FALSE)+
  scale_x_discrete(labels=c("FMT-Healthy", "FMT-Lean NAFLD", "FMT-Obese NAFLD"))+
  theme(axis.text.x = element_text(size = 18, vjust = 0.5, hjust = 0.5, angle = 45))+
  theme(axis.text.y = element_text(size = 20))+
  theme(axis.title.y = element_text(size = 20, vjust = 4))+
  theme(plot.margin=unit(rep(2,4),'lines'))

ggsave(p_human,filename = "shannon_human.pdf",width = 6.5,height = 8)
ggsave(p_mice,filename = "shannon_mice.pdf",width = 6.5,height = 8)
