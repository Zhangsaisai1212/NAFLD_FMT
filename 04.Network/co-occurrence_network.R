# 设置工作目录：请修改下方目录或在Rstudio的Session菜单中选择下载测试数据所在的目录
# setwd("~/Downloads/chenliang")

# 安装需要的包，默认不安装，没安装过的请取消如下注释
# install.packages("igraph")
# install.packages("psych")

# 加载包
library(igraph)
library(psych)
library(RColorBrewer) #颜色包
library(readr)
library(openxlsx)

setwd("/Volumes/My Passport/zhangdw/Lenovo computer/20190716 backup/Bioinformatics/Project/Hein/20200307/scale_reads_to_35M/network_in_mice/taxa/species")

healthy <- read.xlsx("healthy_FMT.xlsx", sheet = 1)
rownames(healthy) <- healthy$clade_name
healthy <- healthy[,-1]

# 计算OTU间两两相关系数矩阵
# 数据量小时可以用psych包corr.test求相关性矩阵，数据量大时，可应用WGCNA中corAndPvalue, 但p值需要借助其他函数矫正
occor <- corr.test(t(healthy),use="pairwise",method="spearman",adjust="BH",alpha=.05)
occor.r <- occor$r # 取相关性矩阵R值
occor.p <- occor$p # 取相关性矩阵p值
occor.r[is.na(occor.r)] <- 0
occor.p[is.na(occor.p)] <- 0

# 确定物种间存在相互作用关系的阈值，将相关性R矩阵内不符合的数据转换为0
occor.r[occor.p>0.05|abs(occor.r)<0.3] = 0 

# 将occor.r保存为csv文件
write.csv(occor.r,file="occurrence_healthy_0.05_0.3.csv")


######################################################################
lean <- read.xlsx("lean_NAFLD_FMT.xlsx", sheet = 1)
rownames(lean) <- lean$clade_name
lean <- lean[,-1]

# 计算OTU间两两相关系数矩阵
# 数据量小时可以用psych包corr.test求相关性矩阵，数据量大时，可应用WGCNA中corAndPvalue, 但p值需要借助其他函数矫正
occor <- corr.test(t(lean),use="pairwise",method="spearman",adjust="BH",alpha=.05)
occor.r <- occor$r # 取相关性矩阵R值
occor.p <- occor$p # 取相关性矩阵p值
occor.r[is.na(occor.r)] <- 0
occor.p[is.na(occor.p)] <- 0

# 确定物种间存在相互作用关系的阈值，将相关性R矩阵内不符合的数据转换为0
occor.r[occor.p>0.05|abs(occor.r)<0.3] = 0 

# 将occor.r保存为csv文件
write.csv(occor.r,file="occurrence_lean_0.05_0.3.csv")


#################################################################################
obese <- read.xlsx("obese_NAFLD_FMT.xlsx", sheet = 1)
rownames(obese) <- obese$clade_name
obese <- obese[,-1]

# 计算OTU间两两相关系数矩阵
# 数据量小时可以用psych包corr.test求相关性矩阵，数据量大时，可应用WGCNA中corAndPvalue, 但p值需要借助其他函数矫正
occor <- corr.test(t(obese),use="pairwise",method="spearman",adjust="BH",alpha=.05)
occor.r <- occor$r # 取相关性矩阵R值
occor.p <- occor$p # 取相关性矩阵p值
occor.r[is.na(occor.r)] <- 0
occor.p[is.na(occor.p)] <- 0

# 确定物种间存在相互作用关系的阈值，将相关性R矩阵内不符合的数据转换为0
occor.r[occor.p>0.05|abs(occor.r)<0.3] = 0 

# 将occor.r保存为csv文件
write.csv(occor.r,file="occurrence_obese_0.05_0.3.csv")

