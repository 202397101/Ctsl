############################
# Figure2A
############################

library('Seurat')
library('tidyverse')
library('ggplot2')
library('patchwork')
library(stringr)
library(dplyr)
library(ggplot2)
#remotes::install_version("Matrix", "1.6.1")
library(Matrix)
library(celldex)
library(dbplyr)
library(dittoSeq)
library(colorRamp2)
library(viridis)
library(randomcoloR)
library(GPTCelltype)
library(openai)
library(cowplot)
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
library(pheatmap)
library(RColorBrewer)
library(GEOquery)
library(stringr)
library(tidyverse)
library(clusterProfiler)
library(org.Mm.eg.db)
library(gridExtra)
library(stringr)
library(ggvenn)
library(msigdbr)
library(fgsea)
library(rWikiPathways)
library(GOfuncR)
library(biomaRt)
#library(devEMF)

mycol = c("mediumpurple","orange","hotpink","cadetblue","gold","maroon","#A9E469","dodgerblue")
grcol = c('lightgoldenrod1', 'mediumturquoise')

DimPlot(s.integrated, reduction='umap', group.by='cachexia', cols = grcol) + ggtitle('') + 
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5), 
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none") 

DimPlot(s.integrated, group.by = 'cellType', label=F, cols = mycol) + ggtitle('')+ 
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5), 
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none") 


############################
# Figure2B, C
############################

library('Seurat')
library('tidyverse')
library('ggplot2')
library('patchwork')
library(stringr)
library(dplyr)
library(ggplot2)
#remotes::install_version("Matrix", "1.6.1")
library(Matrix)
library(celldex)
library(dbplyr)
library(dittoSeq)
library(colorRamp2)
library(viridis)
library(randomcoloR)
library(GPTCelltype)
library(openai)
library(cowplot)
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
library(pheatmap)
library(RColorBrewer)
library(GEOquery)
library(stringr)
library(tidyverse)
library(clusterProfiler)
library(org.Mm.eg.db)
library(gridExtra)
library(stringr)
library(ggvenn)
library(msigdbr)
library(fgsea)
library(rWikiPathways)
library(GOfuncR)
library(biomaRt)
#library(devEMF)

mycol = c("mediumpurple","orange","hotpink","cadetblue","gold","maroon","#A9E469","dodgerblue")
grcol = c('lightgoldenrod1', 'mediumturquoise')

lung.pop = s.integrated@meta.data %>% 
  group_by(group, cellType) %>% 
  summarise(counts = n()) %>% 
  mutate(total.counts = sum(counts)) %>% 
  mutate(percent = counts/total.counts*100) %>%
  filter(group == "Cachexia_Lung") %>%
  arrange(desc(percent))

muscle.pop = s.integrated@meta.data %>% 
  group_by(group, cellType) %>% 
  summarise(counts = n()) %>% 
  mutate(total.counts = sum(counts)) %>% 
  mutate(percent = counts/total.counts*100) %>%
  filter(group == "Cachexia_Muscle") %>%
  arrange(desc(percent))

ggplot(lung.pop, aes(x = reorder(cellType, percent), y=percent, fill=cellType)) + 
  coord_flip()+
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  geom_text(aes(label = paste0(round(lung.pop$percent, 1), "%")), position = position_dodge(1), size = 4, hjust = 0.5, vjust = 0.5) +
  labs(x = "cell type", y = "cell counts (%)", title = ' ',fill = "") +
  scale_fill_manual(values = mycol) + theme_bw() + theme_minimal() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.title.y = element_blank(),  
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "left",
        plot.title = element_text(size = 10, hjust = 0.5))

ggplot(muscle.pop, aes(x = reorder(cellType, percent), y=percent, fill=cellType)) + 
  coord_flip()+
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  geom_text(aes(label = paste0(round(muscle.pop$percent, 1), "%")), position = position_dodge(1), size = 4, hjust = 0.5, vjust = 0.5) +
  labs(x = "cell type", y = "cell counts (%)", title = ' ',fill = "") +
  scale_fill_manual(values = mycol) + theme_bw() + theme_minimal() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.title.y = element_blank(), 
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "left",
        plot.title = element_text(size = 10, hjust = 0.5))



############################
# Figure 2D-1
############################

library('Seurat')
library('tidyverse')
library('ggplot2')
library('patchwork')
library(stringr)
library(dplyr)
library(ggplot2)
#remotes::install_version("Matrix", "1.6.1")
library(Matrix)
library(celldex)
library(dbplyr)
library(dittoSeq)
library(colorRamp2)
library(viridis)
library(randomcoloR)
library(GPTCelltype)
library(openai)
library(cowplot)
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
library(pheatmap)
library(RColorBrewer)
library(GEOquery)
library(stringr)
library(tidyverse)
library(clusterProfiler)
library(org.Mm.eg.db)
library(gridExtra)
library(stringr)
library(ggvenn)
library(msigdbr)
library(fgsea)
library(GOfuncR)
library(biomaRt)



s.integrated@meta.data$cachexia = factor(s.integrated@meta.data$cachexia, levels = c("muscle", "lung"))

FeaturePlot(s.integrated, features = 'Itga1', order = TRUE, split.by = 'cachexia', min.cutoff = 0, max.cutoff = 'q90', by.col = F) & theme(legend.position = "right", legend.title = element_text(size = 12), legend.text = element_text(size = 10), axis.ticks = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), strip.text.y = element_blank()) 

FeaturePlot(s.integrated, features = 'Cd8a', order = TRUE, split.by = 'cachexia', min.cutoff = 0, max.cutoff = 0.08, by.col = F) & theme(legend.position = "right", legend.title = element_text(size = 12), legend.text = element_text(size = 10), axis.ticks = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())


genes = c('Cd3d','Cd3e','Cd3g')
FeaturePlot(s.integrated, features = genes, order = TRUE, split.by = 'cachexia', min.cutoff = 0, max.cutoff = 'q81', by.col = F) & theme(legend.position = "right", legend.title = element_text(size = 12), legend.text = element_text(size = 10), axis.ticks = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())




############################
# Figure 2D-2
############################

genes = c('Prf1', 'Gzma', 'Gzmb')
s.integrated@meta.data$cachexia = factor(s.integrated@meta.data$cachexia, levels = c("muscle", "lung"))

FeaturePlot(s.integrated, features = genes, order = TRUE, split.by = 'cachexia', min.cutoff = 0, max.cutoff = 'q95', by.col = F) & theme(legend.position = "right", legend.title = element_text(size = 12), legend.text = element_text(size = 10), axis.ticks = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())




############################
# Figure 2E
############################

library('Seurat')
library('tidyverse')
library('ggplot2')
library('patchwork')
library(stringr)
library(dplyr)
library(ggplot2)
#remotes::install_version("Matrix", "1.6.1")
library(Matrix)
library(celldex)
library(dbplyr)
library(dittoSeq)
library(colorRamp2)
library(viridis)
library(randomcoloR)
library(GPTCelltype)
library(openai)
library(cowplot)
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
library(pheatmap)
library(RColorBrewer)
library(GEOquery)
library(stringr)
library(tidyverse)
library(clusterProfiler)
library(org.Mm.eg.db)
library(gridExtra)
library(stringr)
library(ggvenn)
library(msigdbr)
library(fgsea)
library(GOfuncR)
library(biomaRt)
#library(devEMF)


# T cell vs Cd49a
tcell.vs.tlikel.updeg = tcell.vs.tlike[tcell.vs.tlike$avg_log2FC > 0 & tcell.vs.tlike$p_val < 0.01 & !is.na(tcell.vs.tlike$p_val), ]
tcell.vs.tlike.ora = enrichGO(gene = rownames(tcell.vs.tlikel.updeg), OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "BP", pvalueCutoff = 1, pAdjustMethod = "fdr", qvalueCutoff = 1, minGSSize = 10, maxGSSize = 500, readable = T)
tcell.vs.tlike.ora = as.data.frame(gofilter(tcell.vs.tlike.ora, level = 4))
tcell.top = head(tcell.vs.tlike.ora, 5)
tcell.top$p.adjust = -log10(tcell.top$p.adjust)
tcell.top = tcell.top[order(tcell.top$p.adjust),]

par(mai = c(0.5,0.1,0.1,1))
bp = barplot(tcell.top$p.adjust, xlim = c(0, 28), horiz = T, xaxt = 'n', yaxt = 'n', 
             xlab = "", names.arg = NA, 
             width = 0.7, border = NA, col = "lightsalmon", main = " ", cex.main = 1.5)
axis(1, at=seq(0,28,4), labels = seq(0,28,4), cex.axis=1, las=1)
abline(v=seq(4,28,4), lty=3, col= "dimgrey")
bp = barplot(tcell.top$p.adjust, xlim = c(0, 28), horiz = T, xaxt = 'n', yaxt = 'n', 
             xlab = "", names.arg = NA, 
             width = 0.7, border = NA, col = "lightsalmon", main = " ", cex.main = 1.5, add = T)
abline(v=0, lty=1)
text(x=0.2, y=bp ,labels=tcell.top$Description, col = "black", xpd=T, cex=1.2, adj=0)


# Cd49a vs T cell
tlike.vs.tcell.updeg = tlike.vs.tcell[tlike.vs.tcell$avg_log2FC > 0 & tlike.vs.tcell$p_val < 0.01 & !is.na(tlike.vs.tcell$p_val), ]
tlike.vs.tcell.ora = enrichGO(gene = rownames(tlike.vs.tcell.updeg), OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "BP", pvalueCutoff = 1, pAdjustMethod = "fdr", qvalueCutoff = 1, minGSSize = 10, maxGSSize = 500, readable = T)
tlike.vs.tcell.ora = as.data.frame(gofilter(tlike.vs.tcell.ora, level = 4))
tlike.top = head(tlike.vs.tcell.ora, 5)
tlike.top$p.adjust = -log10(tlike.top$p.adjust)
tlike.top = tlike.top[order(tlike.top$p.adjust),]

par(mai = c(0.5,0.1,0.1,1))
bp = barplot(tlike.top$p.adjust, xlim = c(0, 12), horiz = T, xaxt = 'n', yaxt = 'n', 
             xlab = "", names.arg = NA, 
             width = 0.7, border = NA, col = "lightsalmon", main = "", cex.main = 1.5)
axis(1, at=seq(0,12,4), labels = seq(0,12,4), cex.axis=1, las=1)
abline(v=seq(4,12,4), lty=3, col= "dimgrey")
bp = barplot(tlike.top$p.adjust, xlim = c(0, 12), horiz = T, xaxt = 'n', yaxt = 'n', 
             xlab = "", names.arg = NA, 
             width = 0.7, border = NA, col = "lightsalmon", main = "", cex.main = 1.5, add = T)
abline(v=0, lty=1)
text(x=0.2, y=bp ,labels=tlike.top$Description, col = "black", xpd=T, cex=1.2, adj=0)



############################
# Figure 2F
############################

cluster.updeg = list()
i=2
for (i in 1:length(markers.li)){
  deg = markers.li[[i]]  
  updeg = lapply(deg, function(j) j[j$avg_log2FC > 2 & j$p_val_adj < 0.01 & !is.na(j$p_val_adj), ])
  cluster.updeg[[names(markers.li)[i]]] = updeg
}

cluster.upKEGG = list()
i=2
j=1
for (i in 1:length(cluster.updeg)){
  updeg = cluster.updeg[[i]]
  kegg = lapply(updeg, function(j) {
    entrez = bitr(rownames(j), fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = "org.Mm.eg.db")
    kg = as.data.frame(enrichKEGG(gene = entrez$ENTREZID, organism = 'mmu', pAdjustMethod = "BH", qvalueCutoff = 0.05, minGSSize = 10, maxGSSize = 500))
  })
  cluster.upKEGG[[names(markers.li)[i]]] = kegg
}

muscle.tlike = cluster.upKEGG$muscle$`T-like`
muscle.tlike.top = head(muscle.tlike, 5)
muscle.tlike.top$p.adjust = -log10(muscle.tlike.top$p.adjust)
muscle.tlike.top = muscle.tlike.top[order(muscle.tlike.top$p.adjust),]
muscle.tlike.top$Description = gsub(" - Mus.*", "", muscle.tlike.top$Description)

par(mai = c(0.5,0.1,0.1,1))
bp = barplot(muscle.tlike.top$p.adjust, xlim = c(0, 10), horiz = T, xaxt = 'n', yaxt = 'n', 
             xlab = "", names.arg = NA, 
             width = 0.7, border = NA, col = "lightsalmon", main = "", cex.main = 1.5)
axis(1, at=seq(0,10,2), labels = seq(0,10,2), cex.axis=1, las=1)
abline(v=seq(2,10,2), lty=3, col= "dimgrey")
bp = barplot(muscle.tlike.top$p.adjust, xlim = c(0, 10), horiz = T, xaxt = 'n', yaxt = 'n', 
             xlab = "", names.arg = NA, 
             width = 0.7, border = NA, col = "lightsalmon", main = "", cex.main = 1.5, add = T)
abline(v=0, lty=1)
text(x=0.2, y=bp ,labels=muscle.tlike.top$Description, col = "black", xpd=T, cex=1.2, adj=0)


muscle.tcell = cluster.upKEGG$muscle$`T cell`
muscle.tcell.top = head(muscle.tcell, 5)
muscle.tcell.top$p.adjust = -log10(muscle.tcell.top$p.adjust)
muscle.tcell.top = muscle.tcell.top[order(muscle.tcell.top$p.adjust),]
muscle.tcell.top$Description = gsub(" - Mus.*", "", muscle.tcell.top$Description)

par(mai = c(0.5,0.1,0.1,1))
bp = barplot(muscle.tcell.top$p.adjust, xlim = c(0, 14), horiz = T, xaxt = 'n', yaxt = 'n', 
             xlab = "", names.arg = NA, 
             width = 0.7, border = NA, col = "lightsalmon", main = "", cex.main = 1.5)
axis(1, at=seq(0,14,2), labels = seq(0,14,2), cex.axis=1, las=1)
abline(v=seq(2,14,2), lty=3, col= "dimgrey")
bp = barplot(muscle.tcell.top$p.adjust, xlim = c(0, 14), horiz = T, xaxt = 'n', yaxt = 'n', 
             xlab = "", names.arg = NA, 
             width = 0.7, border = NA, col = "lightsalmon", main = "", cex.main = 1.5, add = T)
abline(v=0, lty=1)
text(x=0.2, y=bp ,labels=muscle.tcell.top$Description, col = "black", xpd=T, cex=1.2, adj=0)


