#############################################
# Figure 7A GTEx
#############################################

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
library(homologene)
library(GSVA)
library(fgsea)
library(DESeq2)
library(edgeR)
library(dplyr)
library(org.Hs.eg.db)


colnames(muscle.rcm) = gsub("\\.", "-", colnames(muscle.rcm))
rownames(muscle.rcm) = muscle.rcm$Name
ginfo = muscle.rcm[,c(1:2)]
muscle.rcm = muscle.rcm[,-c(1:2)]

ann = read.table("E:/Dropbox/PNU/시스템생물학연구실/DB/GTEx/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt", header = T, row.names = 1, sep = "\t", quote = "\"", fill = T)
tissue = tapply(rownames(ann), ann$SMTSD, FUN = function(i) i) #55 tissues 별로 id 분리.
tissue = tissue$`Muscle - Skeletal`
rcm = muscle.rcm[,colnames(muscle.rcm) %in% tissue]
rcm = na.omit(rcm)


###########################################################################################

new_rownames = ginfo$Description[match(rownames(rcm), ginfo$Name)]
rcm = rcm[!duplicated(new_rownames), ]
rownames(rcm) = new_rownames[!duplicated(new_rownames)]

ctsl = as.numeric(rcm[match('CTSL', rownames(rcm)),])
names(ctsl) = colnames(rcm[match('CTSL', rownames(rcm)),])
ctsl = ctsl[order(ctsl)]
q.ctsl = ceiling(quantile(1:length(ctsl), probs = c(0.5, 0.5))) # low, high 분류하기 위해 삼등분.
low.ctsl = ctsl[1:q.ctsl[1]] 
high.ctsl = ctsl[q.ctsl[2]:length(ctsl)] 
low.id.ctsl = names(low.ctsl) #402
high.id.ctsl = names(high.ctsl) #402

bnip3 = as.numeric(rcm[match('BNIP3', rownames(rcm)),])
names(bnip3) = colnames(rcm[match('BNIP3', rownames(rcm)),])
bnip3 = bnip3[order(bnip3)]
q.bnip3 = ceiling(quantile(1:length(bnip3), probs = c(0.5, 0.5)))
low.bnip3 = bnip3[1:q.bnip3[1]] 
high.bnip3 = bnip3[q.bnip3[2]:length(bnip3)] 
low.id.bnip3 = names(low.bnip3) #402
high.id.bnip3 = names(high.bnip3) #402

low.id = intersect(low.id.ctsl, low.id.bnip3) #280
high.id = intersect(high.id.ctsl, high.id.bnip3)  #280

low.muscle.rcm = rcm[,match(low.id, colnames(rcm))]
high.muscle.rcm = rcm[,match(high.id, colnames(rcm))]

colnames(low.muscle.rcm) = paste0(colnames(low.muscle.rcm), "_low")
colnames(high.muscle.rcm) = paste0(colnames(high.muscle.rcm), "_high")

rcmForDeg = cbind(low.muscle.rcm, high.muscle.rcm) #560

grdf = data.frame(test = "high", control = "low")
sin = data.frame(ID = colnames(rcmForDeg), Group = ifelse(grepl("_low", colnames(rcmForDeg)), "low", "high"))

gr = grdf[1,]
cid = sin$ID[grep(gr[2], sin$Group)]
tid = sin$ID[grep(gr[1], sin$Group)]

dg = rbind(data.frame(sample = cid, fac = 1), data.frame(sample = tid, fac = 2))
ce.idx = match(dg$sample, colnames(rcmForDeg))
dce = round(rcmForDeg[,ce.idx], digits = 0)
dlv = as.factor(dg$fac)
des = model.matrix(~dlv)
keep = filterByExpr(dce, des, group = dlv)
col = data.frame(Condition = dlv)
dds = DESeqDataSetFromMatrix(countData = dce, colData = col, design = ~Condition)

ftdf = data.frame(gene = rownames(rcmForDeg))
mcols(dds) = DataFrame(mcols(dds), ftdf)
deg = DESeq(dds)

deg = results(deg, pAdjustMethod = "fdr", independentFiltering = F)
deg = data.frame(deg)
deg = data.frame(Genes = rownames(dce), deg[match(rownames(dce), rownames(deg)), ], keep, check.names = F)


###################################################################################

# Ctsl - Bnip3 network
ppi = read.csv(file = "E:/Dropbox/PNU/시스템생물학연구실/DB/string/mus_musculus/10090.protein.links.v12.0.txt", sep = " ", header = T, stringsAsFactors = F, quote = "")
pinfo = read.csv(file = "E:/Dropbox/PNU/시스템생물학연구실/DB/string/mus_musculus/10090.protein.info.v12.0.txt", sep = "\t", header = T, stringsAsFactors = F, quote = "")
ppi$protein1 = pinfo$preferred_name[match(ppi$protein1, pinfo$X.string_protein_id)]
ppi$protein2 = pinfo$preferred_name[match(ppi$protein2, pinfo$X.string_protein_id)]
ppi1 = ppi[ppi$combined_score>500 ,]
cb.ppi = ppi1[ppi1$protein1 %in% c('Ctsl','Bnip3') | ppi1$protein2 %in% c('Ctsl','Bnip3'), 1:2]

# IgG vs CON , aCD8 vs IgG (TA muscle)
load(file = sprintf("%s/Rdata/Figure6_RCM_TPM_DEGL.Rdata", dir)) #ta.rcm, ta.tpm, ta.ginfo, ta.sinfo, ta.degl, ga.lung.rcm, ga.lung.tpm, ga.lung.ginfo, ga.lung.sinfo, ga.lung.degl
ta.deg = ta.degl$`IgG vs WT`
ta.deg = ta.deg[!is.na(ta.deg$padj),]
ta.deg = ta.deg[ta.deg$padj < 0.01,]
ta.deg$padj = -log10(ta.deg$padj)
colnames(ta.deg)[7] = "logFDR"
ta.updeg = ta.deg[ta.deg$log2FoldChange > 1 & ta.deg$logFDR > 2 & !is.na(ta.deg$logFDR),]
ta.downdeg = ta.deg[ta.deg$log2FoldChange < -1 & ta.deg$logFDR > 2 & !is.na(ta.deg$logFDR),]

cd8.deg = ta.degl$`Anti-CD8 vs IgG`
cd8.deg = cd8.deg[!is.na(cd8.deg$padj),]
cd8.deg = cd8.deg[cd8.deg$padj < 0.01,]
cd8.deg$padj = -log10(cd8.deg$padj)
colnames(cd8.deg)[7] = "logFDR"
cd8.updeg = cd8.deg[cd8.deg$log2FoldChange > 1 & cd8.deg$logFDR > 2 & !is.na(cd8.deg$logFDR),] 
cd8.downdeg = cd8.deg[cd8.deg$log2FoldChange < -1 & cd8.deg$logFDR > 2 & !is.na(cd8.deg$logFDR),] 

ud = intersect(ta.updeg$Genes, cd8.downdeg$Genes)

ta = cb.ppi[cb.ppi$protein1 %in% ud & cb.ppi$protein2 %in% ud,]
ta = graph_from_edgelist(as.matrix(ta), directed = F)
ta = igraph::simplify(ta)

converted_genes = homologene(names(V(ta)), inTax = 10090, outTax = 9606)
length(names(V(ta)))
converted_genes[26, ] = c('Prkn', 'PRKN')
converted_genes[6,2] = 'CTSL'
V(ta)$name = converted_genes$'9606'[match(names(V(ta)), converted_genes$'10090')]
V(ta)$normalMuscle.log2fc = deg[match(names(V(ta)), deg$Genes), 'log2FoldChange']
V(ta)$normalMuscle.log2fc[is.na(V(ta)$normalMuscle.log2fc)] = 0

createNetworkFromIgraph(ta, "Figure7A")



#############################################
# Figure 7B GTEx
#############################################

library(fgsea)
library(dplyr)
library(openxlsx)

# gsea
deg = deg[match(protein_coding_genes$hgnc_symbol, deg$Genes),]
deg = deg[order(deg$padj),]
idx.use = !is.na(deg$stat)
use = deg[idx.use,]
rk = use$stat
names(rk) = use$Genes
rk = sort(rk, decreasing = T)


gsl = gmtPathways("E:/Dropbox/PNU/시스템생물학연구실/DB/msigdb/c2.all.v2024.1.Hs.symbols.gmt")
gsea = as.data.frame(fgsea(pathways = gsl, stats = rk, minSize = 10, maxSize = 500, nperm=100000))
gsea.df = gsea %>% filter(padj<=0.01) %>% arrange(desc(NES), desc(padj))
gsea.df = gsea.df[grep('WP_',gsea.df$pathway),]
gsea.df = gsea.df %>% slice(c(1:5,(n()-4):n()))
gsea.df$pathway = gsub("WP_","",gsea.df$pathway)
gsea.df$pathway = gsub("_", " ", gsea.df$pathway)
gsea.df$pathway = stringr::str_to_sentence(gsea.df$pathway)

ggplot(gsea.df, aes(x = reorder(pathway, NES), y = NES)) + 
  geom_bar(aes(fill = NES), stat = "identity") +
  geom_text(aes(label = pathway, hjust = -0.1), 
            position = position_stack(vjust = 0),  
            hjust = 0, size = 5)+
  geom_hline(yintercept = 0, color = "black", size = 0.2)+
  geom_vline(xintercept = 5.5, color = "black", size = 0.2) + 
  scale_fill_gradient2(low = "deepskyblue", high = "salmon", 
                       mid = "white", midpoint = 0) + 
  scale_y_continuous(breaks = seq(-3, 3, by = 1), limits = c(-3, 3)) + 
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(), 
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_line(color = "gainsboro", size = 0.3),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 0.5))+
  coord_flip()


#############################################
# Figure 7C GTEx
#############################################

# gsea
deg = deg[match(protein_coding_genes$hgnc_symbol, deg$Genes),]
deg = deg[order(deg$padj),]
idx.use = !is.na(deg$stat)
use = deg[idx.use,]
rk = use$stat
names(rk) = use$Genes
rk = sort(rk, decreasing = T)


gsl = gmtPathways("E:/Dropbox/PNU/시스템생물학연구실/DB/msigdb/c2.all.v2024.1.Hs.symbols.gmt")
gset = gsl[grep('PROTEASOME', names(gsl))]
gset = gset[15]
names(gset)
fname = sprintf("%s/0_reviewer/figure/3-8_metastasis.tif", dir)
labs = list(mt="Proteasome degradation", redgroup.lab="CTSL-BNIP3 high", bluegroup.lab="CTSL-BNIP3 low", mlab="")
length(rk)
xmax = 21000
hcol=c("blue", "white", "red")

gsl = gmtPathways("E:/Dropbox/PNU/시스템생물학연구실/DB/msigdb/c2.all.v2024.1.Hs.symbols.gmt")
gset = gsl[grep('ELECTRON_TRANSPORT_CHAIN', names(gsl))]
fname = sprintf("%s/0_reviewer/figure/3-8_metastasis.tif", dir)
labs = list(mt="Electron Transport Chain Oxphos System In Mitochondria", redgroup.lab="CTSL-BNIP3 high", bluegroup.lab="CTSL-BNIP3 low", mlab="")
length(rk)
xmax = 21000
hcol=c("blue", "white", "red")





