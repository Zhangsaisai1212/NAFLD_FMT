library(corrplot)
library(openxlsx)
library(ggplot2)
library(psych)

setwd("/Users/cq/Google 云端硬盘/HKU/01.Project/03.20200307_NAFLD/scale_reads_to_35M/correlation_taxa_function_clinical/species_level/mice")


dt <- read.xlsx("species_pathway_pvalue_0.05.xlsx")
rownames(dt) <- dt$clade_name
dt <- dt[,-1]
dt <- t(dt)
dt_pathway <- dt[,1:46]
dt_taxa <- dt[,47:76]

# calculate significance
#taxa_clinical_corr <- corr.test(dt_taxt, dt_clinical, method = 'spearman', adjust = "fdr")
pathway_taxa_corr <- corr.test(dt_pathway, dt_taxa, method = 'spearman', adjust = "none")
# correlation
pathway_taxa_corr$r
write.csv(pathway_taxa_corr$r, file = "pathway_taxa_correlation.csv")
# p-value
pathway_taxa_corr$p
write.csv(pathway_taxa_corr$p, file = "pathway_taxa_pvalue.csv")
# filter: r>=0.3 or r<=-0.3
pathway_taxa_corr$r[abs(pathway_taxa_corr$r) < 0.3] <- 0
pathway_taxa_corr$r[is.na(pathway_taxa_corr$r)] <- 0 
pathway_taxa_corr$p[is.na(pathway_taxa_corr$p)] <- 1


pdf(file = "spearman_pathway_pvalue_0.05_taxa_mice_pvalue.pdf",height = 10,width = 12)
#corrplot(t(taxa_clinical_corr$r), method="circle", type="full",order="original",tl.cex = 0.8,tl.col = "black")
col3 <- colorRampPalette(c("blue", "white", "red")) 
corrplot(pathway_taxa_corr$r, p.mat = pathway_taxa_corr$p, method = 'circle',sig.level = 0.05,  
         order="original", insig = "blank", tl.cex=0.8, tl.col = "black",col = col3(200))
dev.off()


pdf(file = "spearman_pathway_taxa_mice_pvalue.pdf",height = 10,width = 12)
par(mar=c(4,4,4,4)+0.5,xpd=TRUE)
#corrplot(t(cor_s), method="circle", type="full",order="original",tl.cex = 0.8,tl.col = "black")
corrplot(pathway_taxa_corr$r, p.mat = pathway_taxa_corr$p, method = 'circle',sig.level = 0.05,  
         order="original", insig = "blank", tl.cex=0.8, tl.col = "black",col = col3(200))
dev.off()

pdf(file = "spearman_pathway_taxa_mice_without_pvalue.pdf",height = 8,width = 8)
par(mar=c(4,4,4,4)+0.5,xpd=TRUE)
#corrplot(t(cor_s), method="circle", type="full",order="original",tl.cex = 0.8,tl.col = "black")
corrplot(pathway_taxa_corr$r, method = 'circle',
         order="original", tl.cex=0.8, tl.col = "black",col = col3(200))
dev.off()
