############################
# Figure 1A
############################

library(readxl)
library(ggplot2)
library(purrr)
library(lubridate)
#BiocManager::install('ggpubr')
library(ggpubr)
library(tidyr)
library(broom)
library(ggrepel)
#install.packages('ggbreak')
library(ggbreak)
library(dplyr)
library(writexl)


ici = as.data.frame(read_excel(path = sprintf("%s/0_reviewer/ICI-weight/PMID39681653.xlsx", dir), col_names = T))
ici = ici[,c(1,2,3,4,6,10,11,12)]
ici = ici[ici$`Grp #` %in% c("1", "2"),]
ici$`Grp #`[ici$`Grp #` == "1"] = "Control"
ici$`Grp #`[ici$`Grp #` == "2"] = "aPD-L1"
colnames(ici)[3] = "Group"
ici = ici[ici$Day %in% c('0','4','8','11','13','15','18'),]
ici$response = NA
ici$response[ici$ID %in% c(408, 433, 462, 500, 412)] = "Responder"
ici$response[ici$ID %in% c(406, 437, 453, 484)] = "Nonresponder"
ici$response[ici$Group == "Control"] = "Control"
ici = ici[!is.na(ici$response),]
ici$TumorWeight = ici$Vol/1000
ici$TumorFreeWeight = ici$Weight - ici$TumorWeight
ici = ici[!ici$response == "Nonresponder",]


# tumor volume
baseline_by_group = ici %>%
  filter(Day == 0) %>%
  group_by(response) %>%
  summarise( baseline = median(Vol, na.rm = TRUE),.groups  = "drop")
  

tm = ici %>%
  filter(Day %in% c(13)) %>%         
  left_join(baseline_by_group, by = "response") %>%
  mutate(delta = (Vol - baseline) / baseline * 100) %>%
  select(-baseline)

#tm = ici %>%
  filter(Day %in% c(13)) %>%         
  left_join(baseline_by_group, by = "response") %>%
  mutate(delta = (Vol - baseline) / baseline * 100) %>%
  dplyr::select(-baseline)

tm_pct = tm %>%
  group_by(response, Day) %>%
  summarise(n = n(),
            median_pct = median(delta, na.rm = TRUE),
            sem = sd(delta, na.rm = TRUE) / sqrt(n), .groups = "drop")

tm = tm %>%
  mutate(Day_f = factor(Day, levels = c(13), labels = c("Day 13")))
tm_pct = tm_pct %>%
  mutate(Day_f = factor(Day, levels = c(13), labels = c("Day 13")))


ggplot() +
  stat_summary(data = tm, aes(x = Day_f, y = delta, fill = response), fun = median, geom = "col",
               position = position_dodge2(width = 0.4, padding = 0.4), width = 0.4, alpha = 0.8) +
  geom_errorbar(data = tm_pct, aes(x = Day_f, ymin = median_pct - sem, ymax = median_pct + sem, group= response), position = position_dodge(width = 0.4), width = 0.1, size = 0.2) +  
  geom_point(data = tm, aes(x = Day_f, y = delta, color = response), position = position_dodge(width = 0.4), size = 0.5, stroke = 0.4, alpha = 1, shape = 1) +
  scale_fill_manual(values = c("Control" = "dimgrey", "Responder" = "darkgreen")) +
  scale_color_manual(values = c("Control" = "black", "Responder" = "green3")) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.3) +
  labs(x = NULL, y = "Tumor Volume (mm³)", fill = "Response") +
  scale_x_discrete(expand = expansion(add = c(0.25,0.25))) +
  coord_cartesian(ylim = c(-500, 1150)) +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 9, color = "black"),
    axis.title = element_text(size = 9, color = "black"),
    legend.position = "none",
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", size = 0.2))


# tumor free body weight
baseline_by_group = ici %>%
  filter(Day == 0) %>%
  group_by(response) %>%
  summarise(baseline = median(TumorFreeWeight, na.rm = TRUE), .groups  = "drop")

bw = ici %>%
  filter(Day %in% c(13)) %>%         
  left_join(baseline_by_group, by = "response") %>%
  mutate( delta = (TumorFreeWeight - baseline) / baseline * 100) %>%
  select(-baseline)

#bw = ici %>%
  filter(Day %in% c(13)) %>%         
  left_join(baseline_by_group, by = "response") %>%
  mutate( delta = (TumorFreeWeight - baseline) / baseline * 100) %>%
  dplyr::select(-baseline)

bw_pct = bw %>%
  group_by(response, Day) %>%
  summarise(n = n(),
            median_pct = median(delta, na.rm = TRUE),
            sem = sd(delta, na.rm = TRUE) / sqrt(n), .groups = "drop")

bw = bw %>%
  mutate(Day_f = factor(Day, levels = c(13), labels = c("Day 13")))
bw_pct = bw_pct %>%
  mutate(Day_f = factor(Day, levels = c(13), labels = c("Day 13")))


ggplot() +
  stat_summary(data = bw, aes(x = Day_f, y = delta, fill = response), fun = median, geom = "col",
               position = position_dodge2(width = 0.4, padding = 0.4), width = 0.4, alpha = 0.8) +
  geom_errorbar(data = bw_pct, aes(x = Day_f, ymin = median_pct - sem, ymax = median_pct + sem, group= response), position = position_dodge(width = 0.4), width = 0.1, size = 0.2) + 
  geom_point(data = bw, aes(x = Day_f, y = delta, color = response), position = position_dodge(width = 0.4), size = 0.5, stroke = 0.4, alpha = 1, shape = 1) +
  scale_fill_manual(values = c("Control" = "dimgrey", "Responder" = "darkgreen")) +
  scale_color_manual(values = c("Control" = "black", "Responder" = "green3")) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.3) +
  labs(x = NULL, y = "% of IBW", fill = "Response") +
  scale_x_discrete(expand = expansion(add = c(0.25, 0.25))) +
  coord_cartesian(ylim = c(-5, 9)) +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 9, color = "black"),
    axis.title = element_text(size = 9, color = "black"),
    legend.position = "none",
    panel.grid = element_blank(),
    axis.line = element_line(color = "black",size = 0.2))



############################
# Figure 1J
############################

library(ggplot2)
library(cluster)
library(ComplexHeatmap)
library(colorRamp2)
library(ggfortify)
library(pheatmap)
library(RColorBrewer)

m = log2(ta.tpm+1)
sds = apply(m, 1, sd)
hist(sds, breaks = 50, col = "skyblue")
uidx = sds>=0.2
m = m[uidx, ]  #13587 genes
xx = t(m)

pca = prcomp(xx)
plot(pca, type = "l")
summary(pca)

xx1 = as.data.frame(xx)
xx1$samples = rep(c("Con", "Anti-CD8", "IgG", "Anti-PD-L1"), c(5,3,4,3))


autoplot(pca, data = xx1, colour = 'samples', frame = T, label = F, label.size = 3) + 
  ggtitle(paste0("PCA ", nrow(m)," genes"))+
  scale_color_manual(values = c("Con" = "black", "Anti-CD8" = "purple", "IgG" = "red2", "Anti-PD-L1" = "blue1"))+
  scale_fill_manual(values = c("Con" = "white", "Anti-CD8" = "white", "IgG" = "white", "Anti-PD-L1" = "white"))+
  theme_bw()+
  theme(legend.position = "none", axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12))



############################
# Figure 1K
############################

library(ggplot2)
library(gridExtra)
library(clusterProfiler)
library(org.Mm.eg.db)
library(gridExtra)
library(stringr)
library(biomaRt)
library(GOfuncR)
library(ComplexHeatmap)
library(colorRamp2)
library(writexl)

ta.deg = ta.degl$`IgG vs WT`
ta.updeg = ta.deg[ta.deg$log2FoldChange > 1 & ta.deg$padj<0.01 & !is.na(ta.deg$padj),] 
ta.downdeg = ta.deg[ta.deg$log2FoldChange < -1 & ta.deg$padj<0.01 & !is.na(ta.deg$padj),] 

#### up pathway
up.ora = enrichGO(gene = ta.updeg$Genes, OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "BP", pvalueCutoff = 1, pAdjustMethod = "fdr", qvalueCutoff = 1, minGSSize = 10, maxGSSize = 500, readable = T)
upgo = as.data.frame(gofilter(up.ora, level = 6))
head(upgo, 10)
upgo$p.adjust = -log10(upgo$p.adjust)
colnames(upgo)[6] = 'logFDR'
upgo = upgo[order(upgo$logFDR),]

par(mai = c(0.5,0.1,0.1,1))
bp = barplot(tail(upgo$logFDR, 5), xlim = c(0, 24), horiz = T, xaxt = 'n', yaxt = 'n', 
             xlab = "", names.arg = NA, 
             width = 0.7, border = NA, col = "lightsalmon", main = " ", cex.main = 1.5)
axis(1, at=seq(0,24,4), labels = seq(0,24,4), cex.axis=1, las=1)
abline(v=seq(4,24,4), lty=3, col= "dimgrey")
bp = barplot(tail(upgo$logFDR, 5), xlim = c(0, 24), horiz = T, xaxt = 'n', yaxt = 'n', 
             xlab = "", names.arg = NA, 
             width = 0.7, border = NA, col = "lightsalmon", main = " ", cex.main = 1.5, add = T)
abline(v=0, lty=1)
text(x=0.2, y=bp ,labels = tail(upgo$Description, 5), col = "black", xpd=T, cex=1.2, adj=0)


#### down pathway
down.ora = enrichGO(gene = ta.downdeg$Genes, OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "BP", pvalueCutoff = 1, pAdjustMethod = "fdr", qvalueCutoff = 1, minGSSize = 10, maxGSSize = 500, readable = T)
downgo = as.data.frame(gofilter(down.ora, level = 6))
head(downgo, 10)
downgo$p.adjust = -log10(downgo$p.adjust)
colnames(downgo)[6] = 'logFDR'
downgo = downgo[order(downgo$logFDR),]

par(mai = c(0.5,0.1,0.1,1))
bp = barplot(tail(downgo$logFDR, 5), xlim = c(0, 24), horiz = T, xaxt = 'n', yaxt = 'n', 
             xlab = "", names.arg = NA, 
             width = 0.7, border = NA, col = "cornflowerblue", main = " ", cex.main = 1.5)
axis(1, at=seq(0,24,4), labels = seq(0,24,4), cex.axis=1, las=1)
abline(v=seq(4,24,4), lty=3, col= "dimgrey")
bp = barplot(tail(downgo$logFDR, 5), xlim = c(0, 24), horiz = T, xaxt = 'n', yaxt = 'n', 
             xlab = "", names.arg = NA, 
             width = 0.7, border = NA, col = "cornflowerblue", main = " ", cex.main = 1.5, add = T)
abline(v=0, lty=1)
text(x=0.2, y=bp ,labels = tail(downgo$Description, 5), col = "black", xpd=T, cex=1.2, adj=0)



############################
# Figure 1L
############################

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
library(ComplexHeatmap)
library(colorRamp2)
library(msigdbr)
library(writexl)


pd.deg = ta.degl$`Anti-PD-L1 vs IgG`
pd.updeg = pd.deg[pd.deg$log2FoldChange > 0 & pd.deg$pvalue<0.005 & !is.na(pd.deg$padj),] 
pd.downdeg = pd.deg[pd.deg$log2FoldChange < 0 & pd.deg$pvalue<0.005 & !is.na(pd.deg$padj),] 

#### up pathway
up.ora = enrichGO(gene = pd.updeg$Genes, OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "BP", pvalueCutoff = 1, pAdjustMethod = "fdr", qvalueCutoff = 1, minGSSize = 10, maxGSSize = 500, readable = T)
up.ora.df = as.data.frame(up.ora)
up.ora.df$p.adjust = -log10(up.ora.df$p.adjust)
colnames(up.ora.df)[6] = 'logFDR'
up.ora.df = up.ora.df[order(up.ora.df$logFDR),]

par(mai = c(0.5,0.1,0.1,1))
bp = barplot(tail(up.ora.df$logFDR, 5), xlim = c(0, 32), horiz = T, xaxt = 'n', yaxt = 'n', 
             xlab = "", names.arg = NA, 
             width = 0.7, border = NA, col = "lightsalmon", main = " ", cex.main = 1.5)
axis(1, at=seq(0,32,4), labels = seq(0,32,4), cex.axis=1, las=1)
abline(v=seq(4,32,4), lty=3, col= "dimgrey")
bp = barplot(tail(up.ora.df$logFDR, 5), xlim = c(0, 32), horiz = T, xaxt = 'n', yaxt = 'n', 
             xlab = "", names.arg = NA, 
             width = 0.7, border = NA, col = "lightsalmon", main = " ", cex.main = 1.5, add = T)
abline(v=0, lty=1)
text(x=0.2, y=bp ,labels = tail(up.ora.df$Description, 5), col = "black", xpd=T, cex=1.2, adj=0)


#### down pathway
down.ora = enrichGO(gene = pd.downdeg$Genes, OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "BP", pvalueCutoff = 1, pAdjustMethod = "fdr", qvalueCutoff = 1, minGSSize = 10, maxGSSize = 500, readable = T)
down.ora.df = as.data.frame(down.ora)
down.ora.df$p.adjust = -log10(down.ora.df$p.adjust)
colnames(down.ora.df)[6] = 'logFDR'
down.ora.df = down.ora.df[order(down.ora.df$logFDR),]

par(mai = c(0.5,0.1,0.1,1))
bp = barplot(tail(down.ora.df$logFDR, 5), xlim = c(0, 5), horiz = T, xaxt = 'n', yaxt = 'n', 
             xlab = "", names.arg = NA, 
             width = 0.7, border = NA, col = "cornflowerblue", main = " ", cex.main = 1.5)
axis(1, at=seq(0,5,1), labels = seq(0,5,1), cex.axis=1, las=1)
abline(v=seq(1,5,1), lty=3, col= "dimgrey")
bp = barplot(tail(down.ora.df$logFDR, 5), xlim = c(0,5), horiz = T, xaxt = 'n', yaxt = 'n', 
             xlab = "", names.arg = NA, 
             width = 0.7, border = NA, col = "cornflowerblue", main = " ", cex.main = 1.5, add = T)
abline(v=0, lty=1)
text(x=0.1, y=bp ,labels = tail(down.ora.df$Description, 5), col = "black", xpd=T, cex=1.2, adj=0)


############################
# Figure 1M
############################

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
library(ComplexHeatmap)
library(colorRamp2)

cd8.deg = ta.degl$`Anti-CD8 vs IgG`
cd8.updeg = cd8.deg[cd8.deg$log2FoldChange > 1 & cd8.deg$padj<0.01 & !is.na(cd8.deg$padj),] 
cd8.downdeg = cd8.deg[cd8.deg$log2FoldChange < -1 & cd8.deg$padj<0.01 & !is.na(cd8.deg$padj),] 

dim(cd8.updeg)
dim(cd8.downdeg)


#### up pathway
up.ora = enrichGO(gene = cd8.updeg$Genes, OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "BP", pvalueCutoff = 1, pAdjustMethod = "fdr", qvalueCutoff = 1, minGSSize = 10, maxGSSize = 500, readable = T)
up.ora.df = as.data.frame(up.ora)

up.top10 = head(up.ora.df, 10)
up.top10$p.adjust = -log10(up.top10$p.adjust)
colnames(up.top10)[6] = "logFDR"
up.top10 = up.top10[order(up.top10$logFDR),]

par(mai = c(0.5,0.1,0.1,1))
bp = barplot(tail(up.top10$logFDR, 5), xlim = c(0, 8), horiz = T, xaxt = 'n', yaxt = 'n', 
             xlab = "", names.arg = NA, 
             width = 0.7, border = NA, col = "lightsalmon", main = " ", cex.main = 1.5)
axis(1, at=seq(0,8,2), labels = seq(0,8,2), cex.axis=1, las=1)
abline(v=seq(2,8,2), lty=3, col= "dimgrey")
bp = barplot(tail(up.top10$logFDR, 5), xlim = c(0, 8), horiz = T, xaxt = 'n', yaxt = 'n', 
             xlab = "", names.arg = NA, 
             width = 0.7, border = NA, col = "lightsalmon", main = " ", cex.main = 1.5, add = T)
abline(v=0, lty=1)
text(x=0.2, y=bp ,labels = tail(up.top10$Description, 5), col = "black", xpd=T, cex=1.2, adj=0)

#### down pathway
down.ora = enrichGO(gene = cd8.downdeg$Genes, OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "BP", pvalueCutoff = 1, pAdjustMethod = "fdr", qvalueCutoff = 1, minGSSize = 10, maxGSSize = 500, readable = T)
down.ora.df = as.data.frame(down.ora)

down.top10 = head(down.ora.df, 10)
down.top10$p.adjust = -log10(down.top10$p.adjust)
colnames(down.top10)[6] = "logFDR"
down.top10$Description[7] = gsub("proteasome-mediated ubiquitin-dependent protein catabolic process", "proteasome-mediated Ub-dependent protein catabolic process", down.top10$Description[7])
down.top10 = down.top10[order(down.top10$logFDR),]

par(mai = c(0.5,0.1,0.1,1))
bp = barplot(tail(down.top10$logFDR, 5), xlim = c(0, 30), horiz = T, xaxt = 'n', yaxt = 'n', 
             xlab = "", names.arg = NA, 
             width = 0.7, border = NA, col = "cornflowerblue", main = " ", cex.main = 1.5)
axis(1, at=seq(0,30,5), labels = seq(0,30,5), cex.axis=1, las=1)
abline(v=seq(5,30,5), lty=3, col= "dimgrey")
bp = barplot(tail(down.top10$logFDR, 5), xlim = c(0,30), horiz = T, xaxt = 'n', yaxt = 'n', 
             xlab = "", names.arg = NA, 
             width = 0.7, border = NA, col = "cornflowerblue", main = " ", cex.main = 1.5, add = T)
abline(v=0, lty=1)
text(x=0.1, y=bp ,labels = tail(down.top10$Description, 5), col = "black", xpd=T, cex=1.2, adj=0)





