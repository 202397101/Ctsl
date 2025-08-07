##########################
# Figure 6A
##########################
dir = "E:/Dropbox/PNU/시스템생물학연구실/data/cachexia"
dir_TCGA = "E:/Dropbox/PNU/시스템생물학연구실/DB/TCGA"

load(file = sprintf("%s/Rdata/Figure6_RCM_TPM_DEGL.Rdata", dir)) #ta.rcm, ta.tpm, ta.ginfo, ta.sinfo, ta.degl, ga.lung.rcm, ga.lung.tpm, ga.lung.ginfo, ga.lung.sinfo, ga.lung.degl
load(file = sprintf("%s/Rdata/TCGA_DEG_tumor_vs_normal.Rdata", dir_TCGA)) # degl


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


luad = degl$LUAD

# Ctsl - Bnip3 network
ppi = read.csv(file = "E:/Dropbox/PNU/시스템생물학연구실/DB/string/mus_musculus/10090.protein.links.v12.0.txt", sep = " ", header = T, stringsAsFactors = F, quote = "")
pinfo = read.csv(file = "E:/Dropbox/PNU/시스템생물학연구실/DB/string/mus_musculus/10090.protein.info.v12.0.txt", sep = "\t", header = T, stringsAsFactors = F, quote = "")
ppi$protein1 = pinfo$preferred_name[match(ppi$protein1, pinfo$X.string_protein_id)]
ppi$protein2 = pinfo$preferred_name[match(ppi$protein2, pinfo$X.string_protein_id)]
ppi1 = ppi[ppi$combined_score>500 ,]
cb.ppi = ppi1[ppi1$protein1 %in% c('Ctsl','Bnip3') | ppi1$protein2 %in% c('Ctsl','Bnip3'), 1:2]


# lung, CAC vs CON
lung.degre = ga.lung.degl$`Cachexia_Lung vs Control`
lung.degre$padj = -log10(lung.degre$padj)
colnames(lung.degre)[7] = "logFDR"
up.lung.degre = lung.degre[lung.degre$log2FoldChange > 1.5 & lung.degre$logFDR > 2 & !is.na(lung.degre$logFDR), ]
down.lung.degre = lung.degre[lung.degre$log2FoldChange < -1.5 & lung.degre$logFDR > 2 & !is.na(lung.degre$logFDR), ]

cb.lung = cb.ppi[cb.ppi$protein1 %in% up.lung.degre$Genes & cb.ppi$protein2 %in% up.lung.degre$Genes,]
cb.lung = graph_from_edgelist(as.matrix(cb.lung), directed = F)
cb.lung = igraph::simplify(cb.lung)


# mouse -> human gene
converted_genes = homologene(names(V(cb.lung)), inTax = 10090, outTax = 9606)
converted_genes[4,2] = 'CTSL'
converted_genes[30, ] = c('Pmaip1', 'PMAIP1')
V(cb.lung)$name = converted_genes$'9606'[match(names(V(cb.lung)), converted_genes$'10090')]

V(cb.lung)$luad.log2fc = luad[match(names(V(cb.lung)), luad$Genes), 'log2FoldChange']
V(cb.lung)$luad.log2fc[is.na(V(cb.lung)$luad.log2fc)] = 0

createNetworkFromIgraph(cb.lung, "Figure6A")



##########################
# Figure 6C
##########################
dir = "E:/Dropbox/PNU/시스템생물학연구실/data/cachexia"
dir_TCGA = "E:/Dropbox/PNU/시스템생물학연구실/DB/TCGA"

load(file = sprintf("%s/Rdata/TCGA_tpml_pan-cancer_v39.Rdata", dir_TCGA)) # tpml, clinical.sinfo

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


# pan-cancer, LUAD, LUSC

tcga.tpm = list(pancancer = tpml$pan_cancer, LUAD = tpml$LUAD, LUSC = tpml$LUSC)
tcga.sinfo = list(pancancer = clinical.sinfo$pan_cancer, LUAD = clinical.sinfo$LUAD, LUSC = clinical.sinfo$LUSC)

i=1
for(i in 1:length(tcga.tpm)){
  tpm = tcga.tpm[[i]]
  sinfo = tcga.sinfo[[i]]
  sinfo = sinfo[sinfo$sample_type_id %in% "01",]
  
  ctsl = tpm[match('CTSL', rownames(tpm)),]
  ctsl = ctsl[order(ctsl)]
  q.ctsl = ceiling(quantile(1:length(ctsl), probs = c(0.5, 0.5))) # low, high 분류하기 위해 삼등분.
  low.ctsl = ctsl[1:q.ctsl[1]] 
  high.ctsl = ctsl[q.ctsl[2]:length(ctsl)] 
  low.id.ctsl = names(low.ctsl)
  high.id.ctsl = names(high.ctsl)
  
  bnip3 = tpm[match('BNIP3', rownames(tpm)),]
  bnip3 = bnip3[order(bnip3)]
  q.bnip3 = ceiling(quantile(1:length(bnip3), probs = c(0.5, 0.5)))
  low.bnip3 = bnip3[1:q.bnip3[1]] 
  high.bnip3 = bnip3[q.bnip3[2]:length(bnip3)] 
  low.id.bnip3 = names(low.bnip3)
  high.id.bnip3 = names(high.bnip3)
  
  low.id = intersect(low.id.ctsl, low.id.bnip3) #1736
  high.id = intersect(high.id.ctsl, high.id.bnip3) #1464
  
  low.idx = which(sinfo$bcr_patient_barcode %in% substr(low.id,1,16))
  high.idx = which(sinfo$bcr_patient_barcode %in% substr(high.id,1,16))
  
  sinfo$CTSL_BNIP3 = NA
  sinfo[low.idx,]$CTSL_BNIP3 = "low"
  sinfo[high.idx,]$CTSL_BNIP3 = "high"
  
  sinfo$dead = sinfo$vital_status == "Dead"
  sinfo$overall_survival = ifelse(sinfo$dead, sinfo$days_to_death, sinfo$days_to_last_follow_up) #dead이면 TRUE -> days_to_death, FALSE -> days_to_last_follow_up.
  
  sinfo = sinfo[sinfo$overall_survival <= 3600, ]
  
  fit = survfit(Surv(overall_survival, dead) ~ CTSL_BNIP3, data = sinfo)
  cox = coxph(Surv(overall_survival, dead) ~ CTSL_BNIP3, data = sinfo)
  hr = round(summary(cox)$conf.int[2], 2)
  
  sp = ggsurvplot(fit, pval=T, pval.coord = c(max(fit$time)-1500, 0.9), risk.table=T, risk.table.col="strata", risk.table.height=0.2, cex.title=8, title = names(tcga.tpm[i]), legend.title=element_text("CTSL-BNIP3 expression level"),legend.labs=c("high", "low"), palette=c("red", "black"))
  
  s.plot = sp$plot + annotate("text", x = max(fit$time)-1000, y = 0.8, label = paste0("HR=", hr), size = 5, color = "black")
  
  ggsave(file = paste0(dir, "/figure_revision/Figure6_", names(tcga.tpm)[i] ,".tiff"), plot = s.plot, width = 12, height = 10, units = 'cm', dpi = 300)
}


##########################
# Figure 6D TCGA-LUAD
##########################
dir = "C:/Dropbox/PNU/시스템생물학연구실/data/cachexia"
dir_TCGA = "C:/Dropbox/PNU/시스템생물학연구실/DB/TCGA"

load(file = sprintf("%s/Rdata/TCGA_rcml_pan-cancer_v39.Rdata", dir_TCGA)) #rcml, clinical.sinfo
load(file = sprintf("%s/Rdata/Figure6_TCGA_high-CTSL-BNIP3_vs_low-CTSL-BNIP3_DEGs.Rdata", dir)) #deg
source("C:/Dropbox/PNU/시스템생물학연구실/data/gluconeogenesis/bulkRNA/Rscripts/11-2_GSEAplot.R")

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
library(GSVA)
library(GOfuncR)
library(msigdbr)
library(ggplot2)
library(org.Hs.eg.db)
library(fgsea)
library(stringr)
library(DESeq2)
library(edgeR)
library(dplyr)


rcm = rcml$LUAD
sinfo = clinical.sinfo$LUAD

sinfo = sinfo[sinfo$sample_type_id %in% "01",]

ctsl = rcm[match('CTSL', rownames(rcm)),]
ctsl = ctsl[order(ctsl)]
q.ctsl = ceiling(quantile(1:length(ctsl), probs = c(0.5, 0.5))) # low, high 분류하기 위해 삼등분.
low.ctsl = ctsl[1:q.ctsl[1]] 
high.ctsl = ctsl[q.ctsl[2]:length(ctsl)] 
low.id.ctsl = names(low.ctsl)
high.id.ctsl = names(high.ctsl)

bnip3 = rcm[match('BNIP3', rownames(rcm)),]
bnip3 = bnip3[order(bnip3)]
q.bnip3 = ceiling(quantile(1:length(bnip3), probs = c(0.5, 0.5)))
low.bnip3 = bnip3[1:q.bnip3[1]] 
high.bnip3 = bnip3[q.bnip3[2]:length(bnip3)] 
low.id.bnip3 = names(low.bnip3)
high.id.bnip3 = names(high.bnip3)

low.id = intersect(low.id.ctsl, low.id.bnip3) #185
high.id = intersect(high.id.ctsl, high.id.bnip3)  #184

low.idx = sinfo[sinfo$bcr_patient_barcode %in% substr(low.id,1,16),] #156
high.idx = sinfo[sinfo$bcr_patient_barcode %in% substr(high.id,1,16),] #183

low.rcm = rcm[,match(rownames(low.idx), colnames(rcm))]
high.rcm = rcm[,match(rownames(high.idx), colnames(rcm))]

colnames(low.rcm) = paste0(colnames(low.rcm), "_low")
colnames(high.rcm) = paste0(colnames(high.rcm), "_high")

rcmForDeg = cbind(low.rcm, high.rcm) #339

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

#save(deg, file = sprintf("%s/Rdata/Figure6_TCGA_high-CTSL-BNIP3_vs_low-CTSL-BNIP3_DEGs.Rdata", dir))

########################### GSEA
deg = deg[order(deg$padj),]

idx.use = !is.na(deg$stat)
use = deg[idx.use,]
rk = use$stat
names(rk) = use$Genes
rk = rk[!is.na(rk)]
rk = sort(rk, decreasing = T)
sum(duplicated(names(rk)))
rk = rk[!duplicated(names(rk))]


gsl = gmtPathways("C:/Dropbox/PNU/시스템생물학연구실/DB/msigdb/c2.all.v2024.1.Hs.symbols.gmt")
gsea = as.data.frame(fgsea(pathways = gsl, stats = rk, minSize = 10, maxSize = 500, nperm=100000))
gsea.df = gsea %>% filter(padj<=0.01) %>% arrange(desc(NES))

gsea.alonso = gsea[grep('ALONSO_', gsea$pathway),]

gset = gsl[grep('ALONSO_METASTASIS_UP', names(gsl))]
fname = sprintf("%s/0_reviewer/figure/3-8_metastasis.tif", dir)
labs = list(mt="Metastasis", redgroup.lab="CTSL-BNIP high", bluegroup.lab="CTSL-BNIP3 low", mlab="")
length(rk)
xmax = 20000
gseaPlot(fname, rk, gset, labs, xmax)
hcol=c("blue", "white", "red")



##########################
# Figure 6D our data
##########################
dir = "E:/Dropbox/PNU/시스템생물학연구실/data/cachexia"

load(file = sprintf("%s/Rdata/Figure6_RCM_TPM_DEGL.Rdata", dir)) #ta.rcm, ta.tpm, ta.ginfo, ta.sinfo, ta.degl, ga.lung.rcm, ga.lung.tpm, ga.lung.ginfo, ga.lung.sinfo, ga.lung.degl
source("E:/Dropbox/PNU/시스템생물학연구실/data/gluconeogenesis/bulkRNA/Rscripts/11-2_GSEAplot.R")

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
library(homologene)
library(GSVA)
library(fgsea)


tumor.deg = ga.lung.degl$`Cachexia_Lung vs Control`
tumor.deg = tumor.deg[order(tumor.deg$padj),]

idx.use = !is.na(tumor.deg$stat)
use = tumor.deg[idx.use,]
rk = use$stat
names(rk) = use$Genes
rk = sort(rk, decreasing = T)

gsl = gmtPathways("E:/Dropbox/PNU/시스템생물학연구실/DB/msigdb/c2.all.v2024.1.Hs.symbols.gmt")
gset = gsl[grep('ALONSO_METASTASIS_UP', names(gsl))]
gset = lapply(gset, function(i) {
  gene = homologene(i, inTax = 9606, outTax = 10090)
  symbol = unique(gene$'10090')
})

fname = sprintf("%s/0_reviewer/figure/3-8_metastasis.tif", dir)
labs = list(mt="Metastasis", redgroup.lab="Tumor", bluegroup.lab="CON", mlab="")
length(rk)
xmax = 32500
hcol=c("blue", "white", "red")


