####################
# Figure 3A
####################

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

ta.deg = ta.degl$`IgG vs WT`
ta.updeg = ta.deg[ta.deg$log2FoldChange > 1 & ta.deg$padj < 0.01 & !is.na(ta.deg$padj),] #1988
ta.downdeg = ta.deg[ta.deg$log2FoldChange < -1 & ta.deg$padj < 0.01 & !is.na(ta.deg$padj),] #1953

cd8.deg = ta.degl$`Anti-CD8 vs IgG`
cd8.updeg = cd8.deg[cd8.deg$log2FoldChange > 1 & cd8.deg$padj<0.01 & !is.na(cd8.deg$padj),] #1210
cd8.downdeg = cd8.deg[cd8.deg$log2FoldChange < -1 & cd8.deg$padj<0.01 & !is.na(cd8.deg$padj),] #1408

ga.deg = ga.lung.degl$`Cachexia_GA vs Control`
ga.updeg = ga.deg[ga.deg$log2FoldChange > 1.5 & ga.deg$padj < 0.01& !is.na(ga.deg$padj), ] #711
ga.downdeg = ga.deg[ga.deg$log2FoldChange < -1.5 & ga.deg$padj < 0.01 & !is.na(ga.deg$padj), ] #542

lung.deg = ga.lung.degl$`Cachexia_Lung vs Control`
lung.updeg = lung.deg[lung.deg$log2FoldChange > 1.5 & lung.deg$padj < 0.01& !is.na(lung.deg$padj), ] #2872
lung.downdeg = lung.deg[lung.deg$log2FoldChange < -1.5 & lung.deg$padj < 0.01 & !is.na(lung.deg$padj), ] #4304

ta.inter = intersect(ta.updeg$Genes, cd8.downdeg$Genes) #1303
ga.ta.inter = intersect(ta.inter, ga.updeg$Genes) #511 (Ctsl)
lung.ta.inter = intersect(ta.inter, lung.updeg$Genes) #108 (Ctsl)
total.inter = intersect(ga.ta.inter, lung.updeg$Genes)

ga.ta.interli = list(TA_IgG=ta.updeg$Genes, 'TA_anti-CD8'=cd8.downdeg$Genes, GA=ga.updeg$Genes)
lung.ta.interli = list(TA_IgG=ta.updeg$Genes, 'TA_anti-CD8'=cd8.downdeg$Genes, Lung=lung.updeg$Genes)
inter.li =  list(GA=ga.updeg$Genes, TA_IgG=ta.updeg$Genes, Lung=lung.updeg$Genes, 'TA_anti-CD8'=cd8.downdeg$Genes)


ggvenn(inter.li, show_percentage = F, 
       fill_color = c("mediumturquoise","purple","lightgoldenrod1","dodgerblue"), fill_alpha = 0.6, stroke_color = "black", 
       stroke_size = 0.3, set_name_size = 0, text_color = "black", text_size = 7) 


total.ora = as.data.frame(enrichGO(gene = total.inter, OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "BP", pvalueCutoff = 1, pAdjustMethod = "fdr", qvalueCutoff = 1, minGSSize = 10, maxGSSize = 500, readable = T))

inter.top5 = head(total.ora, 5)
inter.top5$p.adjust = -log10(inter.top5$p.adjust)
colnames(inter.top5)[6] = "log10FDR"
inter.top5 = inter.top5[order(inter.top5$pvalue, decreasing = T),]

par(mai = c(0.5,0.1,1,1))
bp = barplot(tail(inter.top5$log10FDR), xlim = c(0, 3), horiz = T, xaxt = 'n', yaxt = 'n', 
             xlab = "", names.arg = NA, 
             width = 0.7, border = NA, col = "sandybrown", main = " ", cex.main = 1.5)
axis(1, at=seq(0,3,1), labels = c(0,1,2,3), cex.axis=0.9, las=1)
abline(v=0, lty=1)
abline(v=seq(0,3,1), lty=3, col= "dimgrey")
text(x=0.05, y=bp ,labels = tail(inter.top5$Description), col = "black", xpd=T, cex=1.5, adj=0)


####################
# Figure 3B
####################

library(igraph)
library(RCy3)
library(ggplot2)
library(gridExtra)
library(clusterProfiler)
library(org.Mm.eg.db)
library(gridExtra)
library(stringr)
library(biomaRt)
library(GOfuncR)
library(ggvenn)


ppi = read.csv(file = "yourPath/10090.protein.links.v12.0.txt", sep = " ", header = T, stringsAsFactors = F, quote = "")
pinfo = read.csv(file = "yourPath/10090.protein.info.v12.0.txt", sep = "\t", header = T, stringsAsFactors = F, quote = "")
ppi$protein1 = pinfo$preferred_name[match(ppi$protein1, pinfo$X.string_protein_id)]
ppi$protein2 = pinfo$preferred_name[match(ppi$protein2, pinfo$X.string_protein_id)]
ppi1 = ppi[ppi$combined_score>400 & ppi$protein1 %in% total.inter & ppi$protein2 %in% total.inter, 1:2]

inter.graph = graph_from_edgelist(as.matrix(ppi1), directed = F)
inter.graph = igraph::simplify(inter.graph)

cluster = cluster_walktrap(inter.graph)
cluster.size = sizes(cluster)
cluster.size
cluster$membership[grep('Ctsl', cluster$names)] #1

clusters.to.keep = which(cluster.size > 5)
nodes.to.remove = V(inter.graph)$name[!(membership(cluster) %in% clusters.to.keep)]
graph.filtered = delete.vertices(inter.graph, nodes.to.remove)
cluster.filtered = cluster_walktrap(graph.filtered)
sizes(cluster.filtered) #2개 남음.
cluster.filtered$membership[grep('Ctsl', cluster.filtered$names)] #1

######### cytoscape
V(graph.filtered)[cluster.filtered$membership == 1]$cluster = 'cluster1'
V(graph.filtered)[cluster.filtered$membership == 2]$cluster = 'cluster2'

createNetworkFromIgraph(graph.filtered, "Figure4B")

module.kegg = c()
i=2
for (i in 1:length(cluster.filtered)){
  g = cluster.filtered[[i]]
  entrez = bitr(g, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = "org.Mm.eg.db")
  kegg = as.data.frame(enrichKEGG(gene = entrez$ENTREZID, organism = 'mmu', pAdjustMethod = "BH", qvalueCutoff = 0.05, minGSSize = 10, maxGSSize = 500))
  kegg$Description = gsub(" - Mus musculus.*", "", kegg$Description)
  kegg = kegg[order(kegg$pvalue),]
  module.kegg[[i]] = kegg
}
names(module.kegg) = c('module1', 'module2')

# module1 (Ctsl) barplot
mo1 = module.kegg$module1
mo1$p.adjust = -log10(mo1$p.adjust)
colnames(mo1)[8] = "log10FDR"
mo1 = mo1[order(mo1$log10FDR),]
mo1$Description = gsub(" - animal.*", "", mo1$Description)

par(mai = c(0.5,0.1,1,1))
bp = barplot(tail(mo1$log10FDR), xlim = c(0, 2), horiz = T, xaxt = 'n', yaxt = 'n', 
             xlab = "", names.arg = NA, 
             width = 0.7, border = NA, col = "cyan", main = " ", cex.main = 1.5)
axis(1, at=seq(0,2,1), labels = c(0,1,2), cex.axis=0.9, las=1)
abline(v=0, lty=1)
#abline(v=seq(0,2,1), lty=3, col= "dimgrey")
abline(v=-log10(0.05), lty=3, col= "dimgrey")
text(x=0.05, y=bp ,labels = tail(mo1$Description), col = "black", xpd=T, cex=1.5, adj=0)


####################
# Figure 3C
####################

library(pheatmap)
library(biomaRt)
library(ComplexHeatmap)
library(colorRamp2)

mod2 = c('Ctsl', 'Bnip3', 'Fkbp5', 'Hif3a', 'Egln3', 'Psmd8', 'Tnfrsf12a', 'Psma7')

rownames(annot) = annot$ID
muscle.annot = annot[annot$Cell_type == "muscle" & !annot$Case == "Fasted",] 
ml = sml[match(muscle.annot$ID, names(sml))] #muscle인 geo만 가져오기.

resl = lapply(ml, function(m){
  mm = m[match(mod2, m$Symbol),] 
  list(logFC = mm[,colnames(mm) %in% c("log2FoldChange","logFC")], Qval = mm[,colnames(mm) %in% c("adj.P.Val","padj")])
})

logfc = do.call("cbind", lapply(resl, function(l) l$logFC))
padj = do.call("cbind", lapply(resl, function(l) l$Qval))

rownames(logfc) = mod2
rownames(padj) = mod2
logfc[is.na(logfc)] = 0
padj[is.na(padj)] = 1

sn = apply(is.na(logfc), 1, function(v) sum(v, na.rm=T)) 
use_idx = sn!=ncol(logfc)

logfc1 = logfc[use_idx,]
padj1 = padj[use_idx,]
colnames(logfc1) = annot$label[match(colnames(logfc1), annot[,1])]

val.data = sort(rowSums(padj1<0.05), decreasing = T)
logfc1.sort = logfc1[match(names(val.data), rownames(logfc1)),]
logfc1.sort
padj1.sort = padj1[match(names(val.data), rownames(padj1)),]
hm.val.sort = matrix(NA, nc = ncol(logfc1.sort), nr = nrow(logfc1.sort))
hm.val.sort[padj1.sort > 0.05] = ""
hm.val.sort[padj1.sort < 0.05] = "*"
hm.val.sort[padj1.sort < 0.005] = "**"
logfc1.sort[padj1.sort > 1.5] = 1.5
logfc1.sort[padj1.sort< -1.5] = -1.5

color.ht=colorRamp2(c(-1.5,-1.0,0,1.0,1.5), c('midnightblue', 'mediumblue','white', 'firebrick', 'darkred'))
ha = HeatmapAnnotation("# of Data" = anno_barplot(val.data, gp = gpar(fill = 'pink', col = "pink", lwd = 0),  border = T, bar_width = 0.7), height = unit(15, 'mm'))
lg.ht=Legend(title="Log2FC", at=c(-1.5,0,1.5), col_fun=color.ht, border='black', title_position="topcenter")
hm = Heatmap(t(logfc1.sort), show_heatmap_legend = F, heatmap_legend_param = list(title = 'log2FC'), col = color.ht, border = T, cluster_columns = F, cluster_rows = T,
             show_column_names = T, cluster_column_slices = T, column_names_rot = -90, column_names_gp = gpar(fontsize = 11),
             row_names_gp = gpar(fontsize = 10, fontface = "bold"),
             width=ncol(logfc1.sort)*unit(0.8,"cm"), height=nrow(logfc1.sort)*unit(2.7,"mm"),
             cell_fun = function(j, i, x, y, width, height, fill) {grid.text(hm.val.sort[j, i], x, y, gp = gpar(fontsize=10, col="black"))},
             rect_gp = gpar(col = "white", lwd = 1), 
             top_annotation = ha)

draw(hm, annotation_legend_list=packLegend(list=list(lg.ht)), merge_legend=T, heatmap_legend_side="left", annotation_legend_side='right', padding=unit(c(1,0,0,0), "cm"))

           

############
# Figure3E,H
############

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
library(ggcorrplot)
library(openxlsx)
library(corrplot)


# lung, CAC vs CON
lung.degre = ddfl$`Cachexia_Lung vs Control`
lung.degre$padj = -log10(lung.degre$padj)
colnames(lung.degre)[7] = "logFDR"
up.lung.degre = lung.degre[lung.degre$log2FoldChange > 1.5 & lung.degre$logFDR > 2 & !is.na(lung.degre$logFDR), ]
down.lung.degre = lung.degre[lung.degre$log2FoldChange < -1.5 & lung.degre$logFDR > 2 & !is.na(lung.degre$logFDR), ]

V(ctsl.graph)$lung.log2fc = lung.degre[match(names(V(ctsl.graph)), lung.degre$Genes), 'log2FoldChange']
V(ctsl.graph)$lung.log2fc[is.na(V(ctsl.graph)$lung.log2fc)] = 0

V(ctsl.graph)[names(V(ctsl.graph)) %in% up.lung.degre$Genes]$lung.deg = 'updeg' 
V(ctsl.graph)[names(V(ctsl.graph)) %in% down.lung.degre$Genes]$lung.deg = 'downdeg' 


# Correlation
tumor.tpm = ga.lung.tpm[rownames(ga.lung.tpm) %in% c('Ctsl','Mmp7','Mmp13'),]
tumor.tpm = data.frame(t(log2(tumor.tpm+1)))
colnames(tumor.tpm) = gsub("CO7", "Ctrl", colnames(tumor.tpm))
colnames(tumor.tpm) = gsub("-C-", "-CAC-", colnames(tumor.tpm))

tumor.lung.tpm = tumor.tpm[c(9:11,1:5),]
tumor.ga.tpm = tumor.tpm[c(6:8,12:14),]

corrplot(cor(tumor.lung.tpm), method = "circle", col = colorRampPalette(c('midnightblue', 'mediumblue','white', 'firebrick', 'darkred'))(200), addCoef.col = "white",  tl.col = "black", tl.srt = 45)




############
# Figure3F,I
############
###### CAC vs CON (TA muscle)
ta.deg = degl$`IgG vs WT`
ta.deg = ta.deg[!is.na(ta.deg$padj),]
ta.degre = ta.deg[ta.deg$padj < 0.01,]

ta.updeg = ta.deg[ta.deg$log2FoldChange > 1 & ta.deg$padj < 0.01 & !is.na(ta.deg$padj),]
ta.downdeg = ta.deg[ta.deg$log2FoldChange < -1 & ta.deg$padj < 0.01 & !is.na(ta.deg$padj),]

#### cytoscape version
V(ctsl.graph)$ta.log2fc = ta.deg[match(names(V(ctsl.graph)), ta.deg$Genes), 'log2FoldChange']
V(ctsl.graph)$ta.log2fc[is.na(V(ctsl.graph)$ta.log2fc)] = 0

V(ctsl.graph)[names(V(ctsl.graph)) %in% ta.updeg$Genes]$ta.deg = 'updeg' 
V(ctsl.graph)[names(V(ctsl.graph)) %in% ta.downdeg$Genes]$ta.deg = 'downdeg' 


# ctsl-DEG : ORA
ta.ctsl.deg = names(V(ctsl.graph))[names(V(ctsl.graph)) %in% ta.updeg$Genes]
ta.ctsl.deg.entrez = bitr(ta.ctsl.deg, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = "org.Mm.eg.db")
ta.ctsl.deg.kegg = as.data.frame(enrichKEGG(gene = ta.ctsl.deg.entrez$ENTREZID, organism = 'mmu', pAdjustMethod = "BH", qvalueCutoff = 0.05, minGSSize = 10, maxGSSize = 500))
ta.ctsl.deg.kegg$Description = gsub(" - Mus musculus.*", "", ta.ctsl.deg.kegg$Description)
ta.ctsl.deg.kegg = ta.ctsl.deg.kegg[order(ta.ctsl.deg.kegg$p.adjust),]

ta.ctsl.deg.kegg$p.adjust = -log10(ta.ctsl.deg.kegg$p.adjust)
colnames(ta.ctsl.deg.kegg)[8] = "log10FDR"
ta.ctsl.deg.kegg = ta.ctsl.deg.kegg[order(ta.ctsl.deg.kegg$log10FDR),]

par(mai = c(0.5,0.1,1,1))
bp = barplot(tail(ta.ctsl.deg.kegg$log10FDR, 5), xlim = c(0, 5), horiz = T, xaxt = 'n', yaxt = 'n', 
             xlab = "", names.arg = NA, 
             width = 0.7, border = NA, col = "sandybrown", main = " ", cex.main = 1.5)
axis(1, at=seq(0,5,1), labels = c(0,1,2,3,4,5), cex.axis=0.7, las=1)
abline(v=0, lty=1)
abline(v=seq(0,5,1), lty=3, col= "dimgrey")
text(x=0.05, y=bp ,labels = tail(ta.ctsl.deg.kegg$Description, 5), col = "black", xpd=T, cex=1, adj=0)



############
# Figure3G,J
############
###### anti-CD8a vs IgG (TA muscle)
cd8.deg = degl$`Anti-CD8 vs IgG`
cd8.deg = cd8.deg[!is.na(cd8.deg$padj),]
cd8.degre = cd8.deg[cd8.deg$padj < 0.01,]

cd8.updeg = cd8.deg[cd8.deg$log2FoldChange > 1 & cd8.deg$padj<0.01 & !is.na(cd8.deg$padj),] 
cd8.downdeg = cd8.deg[cd8.deg$log2FoldChange < -1 & cd8.deg$padj<0.01 & !is.na(cd8.deg$padj),] 

#### cytoscape version
V(ctsl.graph)$CD8.log2fc = cd8.deg[match(names(V(ctsl.graph)), cd8.deg$Genes), 'log2FoldChange']
V(ctsl.graph)$CD8.log2fc[is.na(V(ctsl.graph)$CD8.log2fc)] = 0

V(ctsl.graph)[names(V(ctsl.graph)) %in% cd8.updeg$Genes]$CD8.deg = 'updeg' 
V(ctsl.graph)[names(V(ctsl.graph)) %in% cd8.downdeg$Genes]$CD8.deg = 'downdeg' 


# ctsl-DEG : ORA
CD8.ctsl.deg = names(V(ctsl.graph))[names(V(ctsl.graph)) %in% cd8.downdeg$Genes]
CD8.ctsl.deg.entrez = bitr(CD8.ctsl.deg, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = "org.Mm.eg.db")
CD8.ctsl.deg.kegg = as.data.frame(enrichKEGG(gene = CD8.ctsl.deg.entrez$ENTREZID, organism = 'mmu', pAdjustMethod = "BH", qvalueCutoff = 0.05, minGSSize = 10, maxGSSize = 500))
CD8.ctsl.deg.kegg$Description = gsub(" - Mus musculus.*", "", CD8.ctsl.deg.kegg$Description)
CD8.ctsl.deg.kegg = CD8.ctsl.deg.kegg[order(CD8.ctsl.deg.kegg$p.adjust),]

CD8.ctsl.deg.kegg$p.adjust = -log10(CD8.ctsl.deg.kegg$p.adjust)
colnames(CD8.ctsl.deg.kegg)[8] = "log10FDR"
CD8.ctsl.deg.kegg = CD8.ctsl.deg.kegg[order(CD8.ctsl.deg.kegg$log10FDR),]

par(mai = c(0.5,0.1,1,1))
bp = barplot(CD8.ctsl.deg.kegg$log10FDR, xlim = c(0, 3), horiz = T, xaxt = 'n', yaxt = 'n', 
             xlab = "", names.arg = NA, 
             width = 0.7, border = NA, col = "lightsteelblue2", main = " ", cex.main = 1.5)
axis(1, at=seq(0,3,1), labels = c(0,1,2,3), cex.axis=0.7, las=1)
abline(v=0, lty=1)
abline(v=seq(0,3,1), lty=3, col= "dimgrey")
text(x=0.05, y=bp ,labels = CD8.ctsl.deg.kegg$Description, col = "black", xpd=T, cex=1, adj=0)



###############################
# Figure 3K
###############################

library(TCGAbiolinks)
library(SummarizedExperiment)
library(stringr)
library(lubridate)
library(ggpubr)
library(gridExtra)
#install.packages("survival")
library(survival)
#install.packages('survminer')
library(survminer)


################################
# LUAD
################################

# ENSG00000135047 CTSL

tcga.tpm = list(LUAD = tpml$LUAD)
tcga.sinfo = list(LUAD = clinical.sinfo$LUAD)

i=1
for(i in 1:length(tcga.tpm)){
  tpm = tcga.tpm[[i]]
  sinfo = tcga.sinfo[[i]]
  sinfo = sinfo[sinfo$sample_type_id %in% "01", ]
  
  ctsl = tpm[grep('CTSL', rownames(tpm)),]
  ctsl = ctsl[order(ctsl)]
  
  q = ceiling(quantile(1:length(ctsl), probs = c(0.5, 0.5))) # low, high 분류하기 위해 삼등분.
  low = ctsl[1:q[1]] 
  high = ctsl[q[2]:length(ctsl)] 
  low.id = names(low)
  high.id = names(high)
  
  low.idx = which(sinfo$bcr_patient_barcode %in% substr(low.id,1,16)) 
  high.idx = which(sinfo$bcr_patient_barcode %in% substr(high.id,1,16)) 
  
  sinfo$CTSL = NA
  sinfo[low.idx,]$CTSL = "low"
  sinfo[high.idx,]$CTSL = "high"
  
  sinfo$dead = sinfo$vital_status == "Dead"
  sinfo$dead = sinfo$vital_status == "Dead"
  sinfo$overall_survival = ifelse(sinfo$dead, sinfo$days_to_death, sinfo$days_to_last_follow_up) #dead이면 TRUE -> days_to_death, FALSE -> days_to_last_follow_up.
  sinfo = sinfo[sinfo$overall_survival <= 3600, ]
  
  fit = survfit(Surv(overall_survival, dead) ~ CTSL, data = sinfo)
  cox = coxph(Surv(overall_survival, dead) ~ CTSL, data = sinfo)
  hr = round(summary(cox)$conf.int[2], 2)
  
  sp = ggsurvplot(fit, pval=F, pval.coord = c(max(fit$time)-1500, 0.9), risk.table=T, risk.table.col="strata", risk.table.height=0.2, cex.title=8, title = names(tcga.tpm[i]), legend.title=element_text("CTSL expression level"),legend.labs=c("high", "low"), palette=c("red", "black"))
  
  s.plot = sp$plot #+ annotate("text", x = max(fit$time)-1000, y = 0.8, label = paste0("HR=", hr), size = 5, color = "black")
  
  ggsave(file = "yourPath/Figure3K.tiff", plot = s.plot, width = 12, height = 10, units = 'cm', dpi = 300)
}



###############################
# Figure 3L
###############################

library(TCGAbiolinks)
library(SummarizedExperiment)
library(stringr)
library(lubridate)
library(ggpubr)
library(gridExtra)


tpm = tpml$LUAD
annot = unlist(lapply(str_split(colnames(tpm), "-"), function(i) i[4]))
annot = str_sub(annot,1,2)
tumor.ann = annot[annot >= '01' & annot <= '09']
tumor.tpm = tpm[,annot %in% tumor.ann]
ctsl= tumor.tpm[match('CTSL', rownames(tumor.tpm)),]
pdl1= tumor.tpm[match('CD274', rownames(tumor.tpm)),]

cor = cor.test(x = ctsl, y = pdl1, method = "pearson")

par(mar=c(1.5,2,2,1), mgp = c(0,0.5,0))
plot(pdl1 ~ ctsl, xlab = "", ylab = "", pch=19, cex = 0.6, col = "plum", main = "", cex.main = 1.5, cex.axis=0.8, cex.lab=0.2)
abline(lm(pdl1 ~ ctsl), col = "red", lwd=2, lty = 3)



###############################
# Figure 3M
###############################

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
library(pheatmap)
library(RColorBrewer)
library(GEOquery)
library(stringr)
library(corrplot)
library(data.table)
library(msigdbr)


gse = getGEO("GSE283829", GSEMatrix = T,  AnnotGPL=TRUE)
gse = gse[[1]]
eset = exprs(gse)
pset = phenoData(gse)
pset.df = pset@data

tpm = read.csv("yourPath/Figure3M.txt", header=T, sep="\t", row.names = 1)
colnames(tpm) = gsub("S","",colnames(tpm))
colnames(tpm) = pset.df$geo_accession[match(colnames(tpm), pset.df$title)]

CompleteResponse = pset.df[grep('CR', pset.df$`characteristics_ch1.2`),]
StableDisease = pset.df[grep('SD', pset.df$`characteristics_ch1.2`),]
ProgressiveDisease = pset.df[grep('PD', pset.df$`characteristics_ch1.2`),]

ginfo[grep('CD274|CTSL', ginfo$gene_name),]
ctsl = data.frame(t(log2(tpm[grep('ENSG00000135047|ENSG00000120217', rownames(tpm)),]+1)))
colnames(ctsl) = c("CTSL", "PD-L1")
ctsl$sample = ifelse(rownames(ctsl) %in% CompleteResponse$geo_accession, '1_CR', 
                     ifelse(rownames(ctsl) %in% StableDisease$geo_accession, '2_SD', 
                            ifelse(rownames(ctsl) %in% ProgressiveDisease$geo_accession, '3_PD', '-')))

corrt = cor.test(x = ctsl[[1]], y = ctsl[[2]], method = "pearson")

par(mar=c(3,3,2,1), mgp = c(2.5,1,0))
plot(ctsl[[2]] ~ ctsl[[1]], xlab = "", ylab = "", pch=16, cex = 1.5, col = "indianred1", main = "", cex.axis=1, ylim=c(5,12), xlim=c(3,10))
abline(lm(ctsl[[2]] ~ ctsl[[1]]), col = "black", lwd=2, lty = 3)



###############################
# Figure 3N
###############################

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

                      
ddf[grep('CTSL', ddf$Genes),]
tpm[grep('CTSL', rownames(tpm)),]
ctsl = data.frame(log2(tpm[grep('CTSL', rownames(tpm)),]+1))
colnames(ctsl) = 'CTSL'
ctsl$sample = '1'
ctsl$sample[grep('Cachexia_', rownames(ctsl))] = '2'


par(mfrow=c(1,1), mar=c(0.2,1.5,1,1),  mgp=c(0,0.2,0), cex.axis=0.5, cex.lab=0.5, cex.main=0.5, tck=-0.02, las=1, bty="l")
boxplot(ctsl[,1] ~ ctsl$sample, data = ctsl, col = "white", xlab = "", ylab = "", ylim = c(4, 5), lwd = 1.0, cex.main = 0.9, boxwex = 0.8, xaxt = "n", horizontal = F, outline = F)
grid(nx = NA, ny = NULL, lty = 3, lwd = 0.8)
for (j in 1:length(unique(ctsl$sample))) {
  mycol = ifelse(unique(ctsl$sample) == "1", "dodgerblue", "mediumvioletred")
  stripchart(ctsl[,1][ctsl$sample == unique(ctsl$sample)[j]] ~ ctsl$sample[ctsl$sample == unique(ctsl$sample)[j]], method = "jitter", jitter = 0.23, vertical = T,  pch = 19, lwd = 1, cex = 0.5, col = mycol[j], at=j, add = T, bg = "black")
}


