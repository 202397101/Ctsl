#################################
# Figure 4L
#################################
# NSCLC
dir = "E:/Dropbox/PNU/시스템생물학연구실/data/cachexia/"

library(biomaRt)
library(rtracklayer)
library(dplyr)
library(stringr)
library(ggfortify)
library(DESeq2)
library(edgeR)
library(openxlsx)
library(ggplot2)
library(ggrepel)
library(gridExtra)
library(ComplexHeatmap)
library(circlize)

load(file = sprintf("%s/Rdata/5-5-3_human-reference.Rdata", dir)) # ginfo
load(file = sprintf("%s/Rdata/5-5-3_RCM_TPM.Rdata", dir)) # rcm, tpm, sinfo
load(file = sprintf("%s/Rdata/5-5-3_DEG.Rdata", dir)) # ddf


cor.tpm = tpm[match(c('CTSL','HLA-B'), rownames(tpm)),]
cor.tpm = data.frame(t(log2(cor.tpm+1)))

rownames_vector = rownames(cor.tpm)
cor.tpm$sample[grepl("Healthy", rownames_vector)] = "CON"
cor.tpm$sample[grepl("NSCLC", rownames_vector)] = "CAC"
cor.tpm$sample = factor(cor.tpm$sample)

s.col = c("CON" = "darkcyan",
          "CAC" = "darkslategray")
cor.tpm$color = s.col[as.character(cor.tpm$sample)]


# PDAC
load(file = sprintf("%s/Rdata/5-8_GSE133523_rcm_tpm_deg.Rdata", dir)) # rcm1, tpm1, gin, ddf

rownames(tpm1) = gin$HGNC[match(rownames(tpm1), gin$Ensembl)]
cor.tpm1 = tpm1[match(c('CTSL','HLA-B'), rownames(tpm1)),]
cor.tpm1 = data.frame(t(log2(cor.tpm1+1)))

rownames_vector = rownames(cor.tpm1)
cor.tpm1$sample[grepl("control", rownames_vector)] = "CON"
cor.tpm1$sample[grepl("cachexia", rownames_vector)] = "CAC"
cor.tpm1$sample = factor(cor.tpm1$sample)

s.col1 = c("CON" = "orchid",
           "CAC" = "maroon")
cor.tpm1$color = s.col1[as.character(cor.tpm1$sample)]

mergedTpm = rbind(cor.tpm, cor.tpm1)
pch_values = ifelse(grepl("CON", mergedTpm$sample), 16, 17)

par(mar=c(5,0,0,0))
plot(1, type="n", axes=FALSE, xlab="", ylab="")
legend("bottom", legend=c("NSCLC_CON", "NSCLC_CAC", "PDAC_CON", "PDAC_CAC"), col= c("darkcyan","darkslategray", "orchid","maroon"), pch=c(16,17,16,17), pt.cex=1.5, horiz=F, x.intersp=1, ncol = 2) 


corrt = cor.test(x = mergedTpm[[2]], y = mergedTpm[[1]], method = "pearson")

tiff(filename = sprintf("%s/figure_revision/Figure4K.tiff", dir), width = 10, height = 10, units = 'cm', res = 300)
par(mar=c(3,3,2,1), mgp = c(2,1,0))
plot(mergedTpm[[2]] ~ mergedTpm[[1]], xlab = "", ylab = "", pch=pch_values, cex = 1.5, col = mergedTpm$color, main = "", cex.main = 1.5,cex.lab=1, xlim=c(4,5.5), ylim = c(5,8), axes = F)
axis(1, lwd=1, cex.axis=1.2)
axis(2, lwd=1, cex.axis=1.2)
box(lwd=1.2)
abline(lm(mergedTpm[[2]] ~ mergedTpm[[1]]), col = "black", lwd=2, lty = 3)
dev.off()


# for source data
df = data.frame(sample = rownames(mergedTpm), CTSL = mergedTpm$CTSL, 'HLA-B' = mergedTpm$HLA.B)
write.xlsx(df, file = sprintf("%s/excel/Fig4K.xlsx", dir), rowNames=T)



#################################
# Figure 4M tumor
#################################
dir = "E:/Dropbox/PNU/시스템생물학연구실/data/cachexia"

load(file = sprintf("%s/Rdata/Figure6_RCM_TPM_DEGL.Rdata", dir)) #ta.rcm, ta.tpm, ta.ginfo, ta.sinfo, ta.degl, ga.lung.rcm, ga.lung.tpm, ga.lung.ginfo, ga.lung.sinfo, ga.lung.degl

library(igraph)
library(RCy3)
library(RandomWalkRestartMH)
library(ggplot2)
library(gridExtra)
library(clusterProfiler)
library(org.Mm.eg.db)
library(gridExtra)
library(stringr)
library(biomaRt)
library(GOfuncR)
library(ggvenn)
library(corrplot)


cor.tpm = ga.lung.tpm[match(c('Ctsl','H2.D1'), rownames(ga.lung.tpm)),]
cor.tpm = data.frame(t(log2(cor.tpm+1)))
cor.tpm = cor.tpm[c(9:11,1:5),]

rownames_vector = rownames(cor.tpm)
cor.tpm$sample[grepl("Ctrl", rownames_vector)] = "CON"
cor.tpm$sample[grepl("CAC", rownames_vector)] = "TB_Tumor"
cor.tpm$sample = factor(cor.tpm$sample)

s.col = c("CON" = "grey",
          "TB_Tumor" = "darkorange")
cor.tpm$color = s.col[as.character(cor.tpm$sample)]
pch_values = ifelse(grepl("CON", cor.tpm$sample), 16, 17)

par(mar=c(5,0,0,0))
plot(1, type="n", axes=FALSE, xlab="", ylab="")
legend("bottom", legend=c("CON", "TB_Tumor"), col= c("grey","darkorange"), pch=c(16,17), pt.cex=1.5, horiz=F, x.intersp=1, ncol = 1) 

corrt = cor.test(x = cor.tpm[[2]], y = cor.tpm[[1]], method = "pearson")

tiff(filename = sprintf("%s/figure_revision/Figure4L-Tumor.tiff", dir), width = 10, height = 10, units = 'cm', res = 300)
par(mar=c(3,3,2,1), mgp = c(2,1,0))
plot(cor.tpm[[2]] ~ cor.tpm[[1]], xlab = "Ctsl", ylab = colnames(cor.tpm)[2], pch=pch_values, cex = 1.5, col = cor.tpm$color, main = "", cex.main = 1.5, cex.axis=1, cex.lab=1, xlim=c(6,10), ylim = c(8,11.5), axes = F)
axis(1, lwd=1, cex.axis=1.2)
axis(2, lwd=1, cex.axis=1.2, at=seq(8,11,1), labels=seq(8,11,1))
box(lwd=1)
abline(lm(cor.tpm[[2]] ~ cor.tpm[[1]]), col = "black", lwd=2, lty = 3)
dev.off()



#################################
# Figure 4M muscle
#################################
dir = "E:/Dropbox/PNU/시스템생물학연구실/data/cachexia"

load(file = sprintf("%s/Rdata/Figure6_RCM_TPM_DEGL.Rdata", dir)) #ta.rcm, ta.tpm, ta.ginfo, ta.sinfo, ta.degl, ga.lung.rcm, ga.lung.tpm, ga.lung.ginfo, ga.lung.sinfo, ga.lung.degl

library(igraph)
library(RCy3)
library(RandomWalkRestartMH)
library(ggplot2)
library(gridExtra)
library(clusterProfiler)
library(org.Mm.eg.db)
library(gridExtra)
library(stringr)
library(biomaRt)
library(GOfuncR)
library(ggvenn)
library(corrplot)


cor.tpm = ta.tpm[match(c('Ctsl','H2-D1'), rownames(ta.tpm)),]
cor.tpm = data.frame(t(log2(cor.tpm+1)))

rownames_vector = rownames(cor.tpm)
cor.tpm$sample[grepl("CON", rownames_vector)] = "CON"
cor.tpm$sample[grepl("IgG", rownames_vector)] = "IgG"
cor.tpm$sample[grepl("PD-L1", rownames_vector)] = "aPD-L1"
cor.tpm$sample[grepl("CD8", rownames_vector)] = "aCD8"
cor.tpm$sample = factor(cor.tpm$sample)

s.col = c("CON" = "grey",
          "IgG" = "red",
          "aPD-L1" = "blue",
          "aCD8" = "purple")
cor.tpm$color = s.col[as.character(cor.tpm$sample)]
pch_values = ifelse(grepl("CON", cor.tpm$sample), 16, 17)

par(mar=c(5,0,0,0))
plot(1, type="n", axes=FALSE, xlab="", ylab="")
legend("bottom", legend=c("CON", "IgG", "aPD-L1", "aCD8"), col= c("grey","red","blue","purple"), pch=c(16,17,17,17), pt.cex=1.5, horiz=F, x.intersp=1, ncol = 1) 


corrt = cor.test(x = cor.tpm[[2]], y = cor.tpm[[1]], method = "pearson")

tiff(filename = sprintf("%s/figure_revision/Figure4L-muscle.tiff", dir), width = 10, height = 10, units = 'cm', res = 300)
par(mar=c(3,3,2,1), mgp = c(2,1,0))
plot(cor.tpm[[2]] ~ cor.tpm[[1]], xlab = "Ctsl", ylab = colnames(cor.tpm)[2], pch=pch_values, cex = 1.5, col = cor.tpm$color, main = "", cex.main = 1.5, cex.axis=1, cex.lab=1, xlim=c(4,9.5), ylim = c(4.7,6.3), axes = F)
axis(1, lwd=1, cex.axis=1.2)
axis(2, lwd=1, cex.axis=1.2)
box(lwd=1)
abline(lm(cor.tpm[[2]] ~ cor.tpm[[1]]), col = "black", lwd=2, lty = 3)
#x.pos = par("usr")[1] + diff(par("usr")[1:2])*0.8
#y.pos = par("usr")[4] - diff(par("usr")[3:4])*0.8
#text(x.pos, y.pos, paste0("p=", ifelse(corrt$p.value<=1e-6, "1e-6", round(corrt$p.value, 5))), cex=1.2, col = "black")
#x.pos1 = par("usr")[1] + diff(par("usr")[1:2])*0.8
#y.pos1 = par("usr")[4] - diff(par("usr")[3:4])*0.9
#text(x.pos1, y.pos1, paste0("r =", round(corrt$estimate, 3)), cex=1.2, col = "black")  
dev.off()


# for source data
df.tumor = data.frame(sample = rownames(cor.tpm), Ctsl = cor.tpm$Ctsl, 'H2-D1' = cor.tpm$H2.D1)
write.xlsx(df.tumor, file = sprintf("%s/excel/Fig4L_tumor.xlsx", dir), rowNames=T)

df.muscle= data.frame(sample = rownames(cor.tpm), Ctsl = cor.tpm$Ctsl, 'H2-D1' = cor.tpm$H2.D1)
write.xlsx(df.muscle, file = sprintf("%s/excel/Fig4L_muscle.xlsx", dir), rowNames=T)


