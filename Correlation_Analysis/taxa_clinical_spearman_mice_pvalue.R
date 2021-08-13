library(corrplot)
library(openxlsx)
library(ggplot2)
library(psych)

setwd("C:\\Users\\dengwei\\Google Drive\\HKU\\01.Project\\03.20200307_NAFLD\\scale_reads_to_35M\\correlation_taxa_function_clinical\\species_level\\mice")

dt <- read.xlsx("species_clinical_partial.xlsx")
rownames(dt) <- dt$clade_name
dt <- dt[,-1]
dt <- t(dt)
dt_taxt <- dt[,1:30]
dt_clinical <- dt[,31:42]

#cor_s <- cor(dt_taxt, dt_clinical, method = "spearman", use="pairwise.complete.obs")
#cor_s[is.na(cor_s)] <- 0 

# calculate significance
#taxa_clinical_corr <- corr.test(dt_taxt, dt_clinical, method = 'spearman', adjust = "fdr")
taxa_clinical_corr <- corr.test(dt_taxt, dt_clinical, method = 'spearman', adjust = "none")
# correlation
taxa_clinical_corr$r
write.csv(taxa_clinical_corr$r, file="taxa_clinical_correlation.csv")
# p-value
taxa_clinical_corr$p
write.csv(taxa_clinical_corr$p, file="taxa_clinical_pvalue.csv")
# filter: r>=0.3 or r<=-0.3
taxa_clinical_corr$r[abs(taxa_clinical_corr$r) < 0.3] <- 0
taxa_clinical_corr$r[is.na(taxa_clinical_corr$r)] <- 0 
taxa_clinical_corr$p[is.na(taxa_clinical_corr$p)] <- 1

#corrplot(t(taxa_clinical_corr$r), method="circle", type="full",order="original",tl.cex = 0.8,tl.col = "black")
col3 <- colorRampPalette(c("blue", "white", "red")) 
corrplot(taxa_clinical_corr$r, p.mat = taxa_clinical_corr$p, method = 'circle',sig.level = 0.05,  
         order="original", insig = "blank", tl.cex=0.8, tl.col = "black",col = col3(200))

pdf(file = "spearman_taxa_clinical_mice_pvalue.pdf",height = 8.2,width = 8)
#par(mar=c(4,4,4,4)+0.5,xpd=TRUE)
corrplot(taxa_clinical_corr$r, p.mat = taxa_clinical_corr$p, method = 'circle',sig.level = 0.05,  
         order="original", insig = "blank", tl.cex=0.8, tl.col = "black",col = col3(200))
dev.off()

pdf(file = "spearman_taxa_clinical_mice_without_pvalue.pdf",height = 8,width = 8)
par(mar=c(4,4,4,4)+0.5,xpd=TRUE)
#corrplot(t(cor_s), method="circle", type="full",order="original",tl.cex = 0.8,tl.col = "black")
corrplot(taxa_clinical_corr$r, method = 'circle',
         order="original", tl.cex=0.8, tl.col = "black",col = col3(200))
dev.off()
