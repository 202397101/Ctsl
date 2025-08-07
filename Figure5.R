######################
# Figure 5A - lung
######################

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
library(pheatmap)
library(biomaRt)
library(ComplexHeatmap)
library(colorRamp2)



ctsl.neighbor = names(V(ctsl.graph))

cor.tpm = ga.lung.tpm[match(ctsl.neighbor, rownames(ga.lung.tpm)),]
cor.tpm = data.frame(t(log2(cor.tpm+1)))
colnames(cor.tpm) = gsub('\\.', '-', colnames(cor.tpm))

lung.cor = cor.tpm[c(9:11,1:5),]

lcor.res = c()
i=1
for (i in colnames(lung.cor)) {  
  ctsl = "Ctsl"  
  
  if(sd(lung.cor[[i]]) ==0 || sd(lung.cor[[ctsl]]) ==0) {
    lcor.res[i]=0
  } else {
    corrt = tryCatch({
      cor.test(x = lung.cor[[i]], y = lung.cor[[ctsl]], method = "pearson")
    }, error = function(e) return(NULL))
    
    if (!is.null(corrt) && !is.na(corrt$estimate)) {
      lcor.res[i] = corrt$estimate
    }
  }
}
lcor.res = lcor.res[!is.na(lcor.res)]  
lcor.res = sort(lcor.res, decreasing = TRUE) 


######### rank plot
lcor.df = data.frame(Gene = names(lcor.res), Correlation = lcor.res)
lcor.df$Rank = rank(-lcor.df$Correlation)  
lcor.df = lcor.df[-1,]
l.top10 = lcor.df[order(lcor.df$Correlation, decreasing = T),][1:10,]


ggplot(lcor.df, aes(x = Rank, y = Correlation)) +
  geom_point(color = "dimgrey", size = 0.5) +  
  geom_line(color = "dimgrey", size = 0.1) +  
  geom_point(data = l.top10, aes(x = Rank, y = Correlation), color = "red", size = 0.5) + 
  #geom_text_repel(data = l.top10, aes(label = Gene), color = "red", size = 3, max.overlaps = 20, box.padding = 0.5, force = 1, segment.size = 0.2, segment.alpha = 0.7,nudge_y = 0.2, nudge_x = 0.5) + 
  geom_vline(xintercept = 0, color = "black", size=1) +
  geom_hline(yintercept = 0, color = "black", size=0.5) +
  geom_hline(yintercept = -1, color = "black", size=0.5) +
  scale_x_continuous(expand = c(0, 0), limits = c(0,65)) +  
  scale_y_continuous(expand = c(0,0), limits = c(-1,1), breaks = seq(-1,1,0.5)) +  
  theme_classic(base_size = 14) +  
  theme_minimal() +
  labs(title = "",
       x = "",
       y = "") +
  theme(legend.position = "none",
        axis.text = element_text(size = 12, color = "black"),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_line(size = 0.2, color = "grey"),
        panel.grid.major.x = element_line(size = 0.2, color = "grey"))


##########################
# Figure 5A - muscle
##########################

################ 1. GA muscle
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
library(ggrepel)



ctsl.neighbor = names(V(ctsl.graph))

cor.tpm = ga.lung.tpm[match(ctsl.neighbor, rownames(ga.lung.tpm)),]
cor.tpm = data.frame(t(log2(cor.tpm+1)))
colnames(cor.tpm) = gsub('\\.', '-', colnames(cor.tpm))

ga.cor = cor.tpm[c(6:8,12:14),]


################ 2. GSE107470

for (i in seq_along(rdata)) {
  file_path = rdata[i]  
  file_name = fname[i]  
  
  temp_env = new.env()
  load(file_path, envir = temp_env)  #
  
  for (j in ls(temp_env)) {
    new_name = paste0(file_name, "_", j) 
    assign(new_name, get(j, envir = temp_env), envir = .GlobalEnv)  
  }
  
  cat("Loaded:", file_name, "\n")
}

GSE107470.annot = GSE107470_annotation[GSE107470_annotation$tissue == "Mouse gastrocnemius muscle" & !GSE107470_annotation$group == "fasted",] #muscle인 것만 골라오기.

GSE107470_em = GSE107470_em[,match(rownames(GSE107470.annot),colnames(GSE107470_em))]
ctsl.mhc1 = GSE107470_ginfo$Ensembl[GSE107470_ginfo$Symbol %in% ctsl.neighbor]

GSE107470_em = GSE107470_em[match(ctsl.mhc1, rownames(GSE107470_em)),]
colnames(GSE107470_em) = GSE107470_sinfo$`group:ch1`[match(colnames(GSE107470_em), rownames(GSE107470_sinfo))]
colnames(GSE107470_em) = paste0(colnames(GSE107470_em), seq_along(colnames(GSE107470_em)))
GSE107470_em = GSE107470_em[,c(1,3,7,9,10,2,4,5,6,8)]
rownames(GSE107470_em) = GSE107470_ginfo$Symbol[match(rownames(GSE107470_em), GSE107470_ginfo$Ensembl)]

GSE107470_em = data.frame(t(log2(GSE107470_em+1)))
d1.cor = GSE107470_em[1:10,]


############# 3. GSE144567

for (i in seq_along(rdata)) {
  file_path = rdata[i]  
  file_name = fname[i]  
  
  temp_env = new.env()
  load(file_path, envir = temp_env)  #
  
  for (j in ls(temp_env)) {
    new_name = paste0(file_name, "_", j) 
    assign(new_name, get(j, envir = temp_env), envir = .GlobalEnv)  
  }
  
  cat("Loaded:", file_name, "\n")
}

GSE144567.annot = GSE144567_annotation[GSE144567_annotation$tissue == "Tibialis anterior skeletal muscle" ,] #muscle인 것만 골라오기.

GSE144567_em = GSE144567_em[,match(rownames(GSE144567.annot),colnames(GSE144567_em))]
ctsl.mhc1 = GSE144567_ginfo[GSE144567_ginfo$Symbol %in% ctsl.neighbor,]

GSE144567_em = GSE144567_em[match(ctsl.mhc1$Ensembl, rownames(GSE144567_em)),]
colnames(GSE144567_em) = GSE144567_sinfo$`title`[match(colnames(GSE144567_em), rownames(GSE144567_sinfo))]
rownames(GSE144567_em) = GSE144567_ginfo$Symbol[match(rownames(GSE144567_em), GSE144567_ginfo$Ensembl)]

GSE144567_em = data.frame(t(log2(GSE144567_em+1)))
d3.cor = GSE144567_em[1:10,]


############# 5. GSE114820

for (i in seq_along(rdata)) {
  file_path = rdata[i]  
  file_name = fname[i]  
  
  temp_env = new.env()
  load(file_path, envir = temp_env)  #
  
  for (j in ls(temp_env)) {
    new_name = paste0(file_name, "_", j) 
    assign(new_name, get(j, envir = temp_env), envir = .GlobalEnv)  
  }
  
  cat("Loaded:", file_name, "\n")
}

GSE114820.annot = GSE114820_annotation[GSE114820_annotation$tissue == "gastrocnemius muscle" ,] 

GSE114820_em = GSE114820_em[,match(rownames(GSE114820.annot),colnames(GSE114820_em))]
ctsl.mhc1 = GSE114820_ginfo[GSE114820_ginfo$Symbol %in% ctsl.neighbor,]

GSE114820_em = GSE114820_em[match(ctsl.mhc1$Ensembl, rownames(GSE114820_em)),]
colnames(GSE114820_em) = GSE114820_sinfo$`title`[match(colnames(GSE114820_em), rownames(GSE114820_sinfo))]
rownames(GSE114820_em) = GSE114820_ginfo$Symbol[match(rownames(GSE114820_em), GSE114820_ginfo$Ensembl)]
GSE114820_em = GSE114820_em[,grep("PBS|4wks", colnames(GSE114820_em))]

GSE114820_em = data.frame(t(log2(GSE114820_em+1)))
d5.cor = GSE114820_em


############# TA 

cor.tpm = ta.tpm[match(ctsl.neighbor, rownames(ta.tpm)),]
cor.tpm = data.frame(t(log2(cor.tpm+1)))
ta.cor = cor.tpm[c(1:5,9:12),]

c.genes = Reduce(intersect, list(colnames(ta.cor), colnames(d1.cor), colnames(d3.cor), colnames(d5.cor), colnames(ga.cor))) #20016
merge.df = do.call(rbind, list(ta.cor[,c.genes], ga.cor[,c.genes], d1.cor[,c.genes], d3.cor[,c.genes], d5.cor[,c.genes]))

cor.res = c()
for (i in colnames(merge.df)) {  
  ctsl = "Ctsl"  
  corrt = tryCatch({
    cor.test(x = merge.df[[i]], y = merge.df[[ctsl]], method = "pearson")
  }, error = function(e) return(NULL))
  
  if (!is.null(corrt) && !is.na(corrt$estimate)) {
    cor.res[i] = corrt$estimate
  }
}
cor.res = cor.res[!is.na(cor.res)]  
cor.res = sort(cor.res, decreasing = TRUE)  


######### rank plot
mcor.df = data.frame(Gene = names(cor.res), Correlation = cor.res)
mcor.df$Rank = rank(-mcor.df$Correlation)  
mcor.df = mcor.df[-1,]
m.top10 = mcor.df[order(mcor.df$Correlation, decreasing = T),][1:10,]


ggplot(mcor.df, aes(x = Rank, y = Correlation)) +
  geom_point(color = "dimgrey", size = 0.5) +  
  geom_line(color = "dimgrey", size = 0.1) +  
  geom_point(data = m.top10, aes(x = Rank, y = Correlation), color = "red", size = 0.5) + 
  #geom_text_repel(data = m.top10, aes(label = Gene), color = "red", size = 3, max.overlaps = 20, box.padding = 0.5, force = 1, segment.size = 0.2, segment.alpha = 0.7,nudge_y = 0.2, nudge_x = 0.5) + 
  geom_vline(xintercept = 0, color = "black", size=1) +
  geom_hline(yintercept = 0, color = "black", size=1) +
  scale_x_continuous(expand = c(0, 0), limits = c(0,65)) +  
  scale_y_continuous(expand = c(0, 0), limits = c(0,1), breaks = seq(0.25,1,0.25)) +  
  theme_classic(base_size = 14) +  
  theme_minimal() +
  labs(title = "",
       x = "",
       y = "") +
  theme(legend.position = "none",
        axis.text = element_text(size = 12, color = "black"),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_line(size = 0.2, color = "grey"),
        panel.grid.major.x = element_line(size = 0.2, color = "grey"))



##########################
# Figure 5B
##########################

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
library(pheatmap)
library(biomaRt)
library(ComplexHeatmap)
library(colorRamp2)
library(ggrepel)



mcor.df$Rank = mcor.df$Rank-1
lcor.df$Rank = lcor.df$Rank-1
rdf = data.frame(Muscle=mcor.df, Tumor=lcor.df[match(mcor.df$Gene, lcor.df$Gene), ])
identical(rdf$Muscle.Gene, rdf$Tumor.Gene)
rdf.sorted = rdf %>% arrange(Muscle.Rank + Tumor.Rank)
top10 = head(rdf.sorted, 10)

ggplot(rdf, aes(x = Muscle.Rank, y = Tumor.Rank, label = Muscle.Gene)) +
  geom_point(aes(color = ifelse(Muscle.Gene == "Bnip3", "red", "darkgrey")), alpha = 0.8, size = 1) +
  geom_point(data = subset(rdf, Muscle.Gene == "Bnip3"), aes(x = Muscle.Rank, y = Tumor.Rank), color = "red", size = 3, alpha = 0.3) +
  #geom_text(data = top10, aes(x = Muscle.Rank, y = Tumor.Rank, label = Muscle.Gene), vjust = -0.5, hjust = 0.5, size = 2, color = "blue") + 
  scale_color_identity() +
  scale_x_continuous(expand = c(0,0), limits = c(0,62), breaks = seq(1, 62, by = 10)) +  
  scale_y_reverse(expand = c(0,0), limits = c(62,0), breaks = seq(1, 62, by = 10)) +  
  theme_minimal() +
  labs(title = "", x = "", y = "") +
  theme(legend.position = "none",
        axis.text = element_text(size = 12, color = "black"),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_line(size = 0.2, color = "grey"),
        panel.grid.major.x = element_line(size = 0.2, color = "grey"),
        plot.margin = margin(5,15,5,5),
        panel.border = element_rect(color = "black", fill = NA, size = 0.7))


                   
##########################
# Figure 5C - tumor
##########################

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


cor.tpm = ga.lung.tpm[match(c('Ctsl','Bnip3'), rownames(ga.lung.tpm)),]
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

par(mar=c(3,3,2,1), mgp = c(2,1,0))
plot(cor.tpm[[2]] ~ cor.tpm[[1]], xlab = "", ylab = "", pch=pch_values, cex = 1, col = cor.tpm$color, main = "", cex.main = 1.5, cex.axis=1, cex.lab=1, xlim=c(6,10), ylim = c(4,11), axes = F)
axis(1, lwd=1, cex.axis=1.2)
axis(2, lwd=1, cex.axis=1.2, at=seq(4,12,2), labels=seq(4,12,2))
box(lwd=1)
abline(lm(cor.tpm[[2]] ~ cor.tpm[[1]]), col = "black", lwd=2, lty = 3)




##########################
# Figure 5C - muscle
##########################

################ 1. GA muscle

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


cor.tpm = ga.lung.tpm[match(c('Ctsl','Bnip3'), rownames(ga.lung.tpm)),]
cor.tpm = data.frame(t(log2(cor.tpm+1)))
colnames(cor.tpm) = gsub('\\.', '-', colnames(cor.tpm))

ga.cor = cor.tpm[c(6:8,12:14),]


################ 2. GSE107470

for (i in seq_along(rdata)) {
  file_path = rdata[i]  
  file_name = fname[i]  
  
  temp_env = new.env()
  load(file_path, envir = temp_env)  #
  
  for (j in ls(temp_env)) {
    new_name = paste0(file_name, "_", j) 
    assign(new_name, get(j, envir = temp_env), envir = .GlobalEnv)  
  }
  
  cat("Loaded:", file_name, "\n")
}

GSE107470.annot = GSE107470_annotation[GSE107470_annotation$tissue == "Mouse gastrocnemius muscle" & !GSE107470_annotation$group == "fasted",] #muscle인 것만 골라오기.

GSE107470_em = GSE107470_em[,match(rownames(GSE107470.annot),colnames(GSE107470_em))]
ctsl.mhc1 = GSE107470_ginfo$Ensembl[GSE107470_ginfo$Symbol %in% c('Ctsl','Bnip3')]

GSE107470_em = GSE107470_em[match(ctsl.mhc1, rownames(GSE107470_em)),]
colnames(GSE107470_em) = GSE107470_sinfo$`group:ch1`[match(colnames(GSE107470_em), rownames(GSE107470_sinfo))]
colnames(GSE107470_em) = paste0(colnames(GSE107470_em), seq_along(colnames(GSE107470_em)))
GSE107470_em = GSE107470_em[,c(1,3,7,9,10,2,4,5,6,8)]
rownames(GSE107470_em) = GSE107470_ginfo$Symbol[match(rownames(GSE107470_em), GSE107470_ginfo$Ensembl)]

GSE107470_em = data.frame(t(log2(GSE107470_em+1)))
d1.cor = GSE107470_em[1:10,]


############# 3. GSE144567

for (i in seq_along(rdata)) {
  file_path = rdata[i]  
  file_name = fname[i]  
  
  temp_env = new.env()
  load(file_path, envir = temp_env)  #
  
  for (j in ls(temp_env)) {
    new_name = paste0(file_name, "_", j) 
    assign(new_name, get(j, envir = temp_env), envir = .GlobalEnv)  
  }
  
  cat("Loaded:", file_name, "\n")
}

GSE144567.annot = GSE144567_annotation[GSE144567_annotation$tissue == "Tibialis anterior skeletal muscle" ,] #muscle인 것만 골라오기.

GSE144567_em = GSE144567_em[,match(rownames(GSE144567.annot),colnames(GSE144567_em))]
ctsl.mhc1 = GSE144567_ginfo[GSE144567_ginfo$Symbol %in% c('Ctsl','Bnip3'),]

GSE144567_em = GSE144567_em[match(ctsl.mhc1$Ensembl, rownames(GSE144567_em)),]
colnames(GSE144567_em) = GSE144567_sinfo$`title`[match(colnames(GSE144567_em), rownames(GSE144567_sinfo))]
rownames(GSE144567_em) = GSE144567_ginfo$Symbol[match(rownames(GSE144567_em), GSE144567_ginfo$Ensembl)]

GSE144567_em = data.frame(t(log2(GSE144567_em+1)))
d3.cor = GSE144567_em[1:10,]


############# 5. GSE114820

for (i in seq_along(rdata)) {
  file_path = rdata[i]  
  file_name = fname[i]  
  
  temp_env = new.env()
  load(file_path, envir = temp_env)  #
  
  for (j in ls(temp_env)) {
    new_name = paste0(file_name, "_", j) 
    assign(new_name, get(j, envir = temp_env), envir = .GlobalEnv)  
  }
  
  cat("Loaded:", file_name, "\n")
}

GSE114820.annot = GSE114820_annotation[GSE114820_annotation$tissue == "gastrocnemius muscle" ,] #muscle인 것만 골라오기.

GSE114820_em = GSE114820_em[,match(rownames(GSE114820.annot),colnames(GSE114820_em))]
ctsl.mhc1 = GSE114820_ginfo[GSE114820_ginfo$Symbol %in% c('Ctsl','Bnip3'),]

GSE114820_em = GSE114820_em[match(ctsl.mhc1$Ensembl, rownames(GSE114820_em)),]
colnames(GSE114820_em) = GSE114820_sinfo$`title`[match(colnames(GSE114820_em), rownames(GSE114820_sinfo))]
rownames(GSE114820_em) = GSE114820_ginfo$Symbol[match(rownames(GSE114820_em), GSE114820_ginfo$Ensembl)]
GSE114820_em = GSE114820_em[,grep("PBS|4wks", colnames(GSE114820_em))]

GSE114820_em = data.frame(t(log2(GSE114820_em+1)))
d5.cor = GSE114820_em


########### plot 

ga.cor
colnames(ga.cor) = gsub('-', '\\.', colnames(ga.cor))
rownames(ga.cor) = paste0(rownames(ga.cor), "_ga")

d1.cor
rownames(d1.cor) = paste0(rownames(d1.cor), "_d1")
rownames(d1.cor)[grep('non-CACS', rownames(d1.cor))] = gsub("non-CACS", "healthy",rownames(d1.cor)[grep('non-CACS', rownames(d1.cor))] )
d1.cor = d1.cor[,c(2,1)]

d3.cor
rownames(d3.cor) = paste0(rownames(d3.cor), "_d3")

d5.cor
rownames(d5.cor) = paste0(rownames(d5.cor), "_d5")

colnames(d5.cor)
colnames(d3.cor)
colnames(d1.cor)
colnames(ga.cor)

total.cor = rbind(ga.cor, d1.cor, d3.cor, d5.cor)
rownames_vector = rownames(total.cor)
total.cor$sample[grepl("CO7", rownames_vector)] = "ga_con"
total.cor$sample[grepl("KR", rownames_vector)] = "ga_cac"
total.cor$sample[grepl("healthy", rownames_vector)] = "d1_con"
total.cor$sample[grepl("CACS", rownames_vector)] = "d1_cac"
total.cor$sample[grepl("Control", rownames_vector)] = "d3_con"
total.cor$sample[grepl("LLC", rownames_vector)] = "d3_cac"
total.cor$sample[grepl("PBS", rownames_vector)] = "d5_con"
total.cor$sample[grepl("4wks", rownames_vector)] = "d5_cac"
total.cor$sample = factor(total.cor$sample)

s.col = c("ga_con" = "gold",
          "ga_cac" = "goldenrod",
          "d1_con" = "skyblue1",
          "d1_cac" = "dodgerblue1",
          "d3_con" = "orchid1",
          "d3_cac" = "violetred1",
          "d5_con" = "green1",
          "d5_cac" = "forestgreen")
total.cor$color = s.col[as.character(total.cor$sample)]


############# TA 

cor.tpm = ta.tpm[match(c('Ctsl','Bnip3'), rownames(ta.tpm)),]
cor.tpm = data.frame(t(log2(cor.tpm+1)))

ta.cor = cor.tpm[c(1:5,9:12),]
ta.cor$sample[grepl("CON", rownames(ta.cor))] = "TA_con"
ta.cor$sample[grepl("IgG", rownames(ta.cor))] = "TA_IgG"
ta.cor$sample = factor(ta.cor$sample)

s.col = c("TA_con" = "gold",
          "TA_IgG" = "goldenrod")
ta.cor$color = s.col[as.character(ta.cor$sample)]

colnames(total.cor)
colnames(ta.cor)
all.cor = rbind(total.cor, ta.cor)
pch_values = ifelse(grepl("con", all.cor$sample), 16, 17)


par(mar=c(5,0,0,0))
plot(1, type="n", axes=FALSE, xlab="", ylab="")
legend("bottom", legend=c("muscle_CON", "muscle_TB", "D1_CON", "D1_TB", "D3_CON", "D3_TB", "D5_CON", "D5_TB"), col= c("gold", "goldenrod","skyblue1", "dodgerblue1","orchid1","violetred1","green1","forestgreen"), pch=rep(c(16,17), times = 5), pt.cex=1.5, horiz=F, x.intersp=0.7, ncol = 1) 


corrt = cor.test(x = cor.tpm[[2]], y = cor.tpm[[1]], method = "pearson")

par(mar=c(3,3,2,1), mgp = c(2,1,0))
plot(all.cor[[2]] ~ all.cor[[1]], xlab = "", ylab = "", pch=pch_values, cex = 1, col = all.cor$color, main = "", cex.main = 1.5, cex.axis=1, cex.lab=1, xlim=c(1,15), ylim = c(3,15), axes = F)
axis(1, lwd=1, cex.axis=1.2, at=seq(1,15,3), labels=seq(1,15,3))
axis(2, lwd=1, cex.axis=1.2, at=seq(3,15,3), labels=seq(3,15,3))
box(lwd=1)
abline(lm(all.cor[[2]] ~ all.cor[[1]]), col = "black", lwd=2, lty = 3)
#x.pos = par("usr")[1] + diff(par("usr")[1:2])*0.8
#y.pos = par("usr")[4] - diff(par("usr")[3:4])*0.8
#text(x.pos, y.pos, paste0("p=", ifelse(corrt$p.value<=1e-6, "1e-6", round(corrt$p.value, 5))), cex=1.2, col = "black")
#x.pos1 = par("usr")[1] + diff(par("usr")[1:2])*0.8
#y.pos1 = par("usr")[4] - diff(par("usr")[3:4])*0.9
#text(x.pos1, y.pos1, paste0("r =", round(corrt$estimate, 3)), cex=1.2, col = "black")  



#################
# Figure 5D
#################

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
library(pheatmap)
library(biomaRt)
library(ComplexHeatmap)
library(colorRamp2)


# lung, CAC vs CON
lung.degre = ga.lung.degl$`Cachexia_Lung vs Control`
lung.degre$padj = -log10(lung.degre$padj)
colnames(lung.degre)[7] = "logFDR"
up.lung.degre = lung.degre[lung.degre$log2FoldChange > 1.5 & lung.degre$logFDR > 2 & !is.na(lung.degre$logFDR), ]
down.lung.degre = lung.degre[lung.degre$log2FoldChange < -1.5 & lung.degre$logFDR > 2 & !is.na(lung.degre$logFDR), ]

V(ctsl.graph)$lung.log2fc = lung.degre[match(names(V(ctsl.graph)), lung.degre$Genes), 'log2FoldChange']
V(ctsl.graph)$lung.log2fc[is.na(V(ctsl.graph)$lung.log2fc)] = 0

V(ctsl.graph)[names(V(ctsl.graph)) %in% up.lung.degre$Genes]$lung.deg = 'updeg' 
V(ctsl.graph)[names(V(ctsl.graph)) %in% down.lung.degre$Genes]$lung.deg = 'downdeg' 

up.lung=up.lung.degre[match(names(V(ctsl.graph)[names(V(ctsl.graph)) %in% up.lung.degre$Genes]), up.lung.degre$Genes),] 
up.lung[order(up.lung$log2FoldChange, decreasing=T),]
up.lung[order(up.lung$logFDR, decreasing=T),]
down.lung=down.lung.degre[match(names(V(ctsl.graph)[names(V(ctsl.graph)) %in% down.lung.degre$Genes]), down.lung.degre$Genes),] 
down.lung[order(down.lung$log2FoldChange, decreasing=T),]
down.lung[order(down.lung$logFDR, decreasing=T),]
up.down.lung = rbind(up.lung, down.lung)



# CAC vs CON (TA muscle)
ta.deg = ta.degl$`IgG vs WT`
ta.deg = ta.deg[!is.na(ta.deg$padj),]
ta.deg = ta.deg[ta.deg$padj < 0.01,]
ta.deg$padj = -log10(ta.deg$padj)
colnames(ta.deg)[7] = "logFDR"

ta.updeg = ta.deg[ta.deg$log2FoldChange > 1 & ta.deg$logFDR > 2 & !is.na(ta.deg$logFDR),]
ta.downdeg = ta.deg[ta.deg$log2FoldChange < -1 & ta.deg$logFDR > 2 & !is.na(ta.deg$logFDR),]

V(ctsl.graph)$ta.log2fc = ta.deg[match(names(V(ctsl.graph)), ta.deg$Genes), 'log2FoldChange']
V(ctsl.graph)$ta.log2fc[is.na(V(ctsl.graph)$ta.log2fc)] = 0

V(ctsl.graph)[names(V(ctsl.graph)) %in% ta.updeg$Genes]$ta.deg = 'updeg' 
V(ctsl.graph)[names(V(ctsl.graph)) %in% ta.downdeg$Genes]$ta.deg = 'downdeg' 

up.igg=ta.updeg[match(names(V(ctsl.graph)[names(V(ctsl.graph)) %in% ta.updeg$Genes]), ta.updeg$Genes),] 
up.igg[order(up.igg$log2FoldChange, decreasing=T),]
up.igg[order(up.igg$logFDR, decreasing=T),]

down.igg=ta.downdeg[match(names(V(ctsl.graph)[names(V(ctsl.graph)) %in% ta.downdeg$Genes]), ta.downdeg$Genes),] 
down.igg[order(down.igg$log2FoldChange, decreasing=T),]
down.igg[order(down.igg$logFDR, decreasing=T),]
up.down.igg = rbind(up.igg,down.igg)



# anti-CD8a vs IgG (TA muscle)
cd8.deg = ta.degl$`Anti-CD8 vs IgG`
cd8.deg = cd8.deg[!is.na(cd8.deg$padj),]
cd8.deg = cd8.deg[cd8.deg$padj < 0.01,]
cd8.deg$padj = -log10(cd8.deg$padj)
colnames(cd8.deg)[7] = "logFDR"

cd8.updeg = cd8.deg[cd8.deg$log2FoldChange > 1 & cd8.deg$logFDR > 2 & !is.na(cd8.deg$logFDR),] 
cd8.downdeg = cd8.deg[cd8.deg$log2FoldChange < -1 & cd8.deg$logFDR > 2 & !is.na(cd8.deg$logFDR),] 

V(ctsl.graph)$CD8.log2fc = cd8.deg[match(names(V(ctsl.graph)), cd8.deg$Genes), 'log2FoldChange']
V(ctsl.graph)$CD8.log2fc[is.na(V(ctsl.graph)$CD8.log2fc)] = 0

V(ctsl.graph)[names(V(ctsl.graph)) %in% cd8.updeg$Genes]$CD8.deg = 'updeg' 
V(ctsl.graph)[names(V(ctsl.graph)) %in% cd8.downdeg$Genes]$CD8.deg = 'downdeg' 

up.cd8=cd8.updeg[match(names(V(ctsl.graph)[names(V(ctsl.graph)) %in% cd8.updeg$Genes]), cd8.updeg$Genes),] 
up.cd8[order(up.cd8$log2FoldChange, decreasing=T),]
up.cd8[order(up.cd8$logFDR, decreasing=T),]

down.cd8=cd8.downdeg[match(names(V(ctsl.graph)[names(V(ctsl.graph)) %in% cd8.downdeg$Genes]), cd8.downdeg$Genes),] 
down.cd8[order(down.cd8$log2FoldChange, decreasing=T),]
down.cd8[order(down.cd8$logFDR, decreasing=T),]
up.down.cd8 = rbind(up.cd8,down.cd8)

ctsl.list = list(lung=up.lung, IgG.TA=up.igg, aCD8.TA=down.cd8)


###############################################
lung.gene = ctsl.list$lung$Genes[-1]
igg.gene = ctsl.list$IgG.TA$Genes[-1]
cd8.gene = ctsl.list$aCD8.TA$Genes[-1]
commonGeneli = list(Lung=lung.gene, TA=igg.gene, 'anti-CD8a'=cd8.gene)

commonGenes = Reduce(intersect, list(lung.gene, igg.gene, cd8.gene))

ggvenn(commonGeneli, show_percentage = F, 
       fill_color = c("lightgoldenrod1","purple","dodgerblue"), fill_alpha = 0.5, stroke_color = "black", 
       stroke_size = 0.3, set_name_size = 0, text_color = "black", text_size = 7) 





#############################################
# Figure 5E,F,G,H
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


# Ctsl - Bnip3 network
ppi = read.csv(file = "yourPath/10090.protein.links.v12.0.txt", sep = " ", header = T, stringsAsFactors = F, quote = "")
pinfo = read.csv(file = "yourPath/10090.protein.info.v12.0.txt", sep = "\t", header = T, stringsAsFactors = F, quote = "")
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

V(cb.lung)$lung.log2fc = lung.degre[match(names(V(cb.lung)), lung.degre$Genes), 'log2FoldChange']
V(cb.lung)$lung.log2fc[is.na(V(cb.lung)$lung.log2fc)] = 0

V(cb.lung)[names(V(cb.lung)) %in% up.lung.degre$Genes]$lung.deg = 'updeg' 
V(cb.lung)[names(V(cb.lung)) %in% down.lung.degre$Genes]$lung.deg = 'downdeg' 

createNetworkFromIgraph(cb.lung, "Figure5E")



# IgG vs CON , aCD8 vs IgG (TA muscle)
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

V(ta)$ta.log2fc = ta.deg[match(names(V(ta)), ta.deg$Genes), 'log2FoldChange']
V(ta)$ta.log2fc[is.na(V(ta)$ta.log2fc)] = 0
V(ta)[names(V(ta)) %in% ta.updeg$Genes]$ta.deg = 'updeg' 
V(ta)[names(V(ta)) %in% ta.downdeg$Genes]$ta.deg = 'downdeg' 

V(ta)$cd8.log2fc = cd8.deg[match(names(V(ta)), cd8.deg$Genes), 'log2FoldChange']
V(ta)$cd8.log2fc[is.na(V(ta)$cd8.log2fc)] = 0
V(ta)[names(V(ta)) %in% cd8.updeg$Genes]$cd8.deg = 'updeg' 
V(ta)[names(V(ta)) %in% cd8.downdeg$Genes]$cd8.deg = 'downdeg' 

createNetworkFromIgraph(ta, "Figure5F")


net.gl = list(Lung = names(V(cb.lung)),
              TA = names(V(ta)))

gofilter.li = list()
go.li = list()
i=1
for (i in 1:length(net.gl)){
  gene = net.gl[[i]]
  go = enrichGO(gene = gene, OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "BP", pvalueCutoff = 1, pAdjustMethod = "fdr", qvalueCutoff = 1, minGSSize = 10, maxGSSize = 500, readable = T)
  gofilter = as.data.frame(gofilter(go, level = 5))
  gofilter$p.adjust = -log10(gofilter$p.adjust)
  colnames(gofilter)[6] = "logFDR"
  gofilter = as.data.frame(gofilter[order(gofilter$logFDR),])
  gofilter.li[[names(net.gl)[i]]] = gofilter
  go.li[[names(net.gl)[i]]] = go
}

lung.go = go.li$Lung@result
lung.go$p.adjust = -log10(lung.go$p.adjust)
colnames(lung.go)[6] = "logFDR"
lung.go = lung.go[order(lung.go$logFDR), ]
lung.go = lung.go[!grepl('positive|negative', lung.go$Description),]
lung.go$Description = str_replace(lung.go$Description, "^(\\w)", toupper)

par(mai = c(0.5,0.1,0.1,1))
bp = barplot(tail(lung.go$logFDR,7), xlim = c(0, 7), horiz = T, xaxt = 'n', yaxt = 'n', 
             xlab = "", names.arg = NA, 
             width = 0.7, border = NA, col = "lightsalmon", main = "", cex.main = 1.5)
axis(1, at=seq(0,7,1), labels = seq(0,7,1), cex.axis=1, las=1)
abline(v=seq(1,7,1), lty=3, col= "dimgrey")
bp = barplot(tail(lung.go$logFDR,7), xlim = c(0, 7), horiz = T, xaxt = 'n', yaxt = 'n', 
             xlab = "", names.arg = NA, 
             width = 0.7, border = NA, col = "lightsalmon", main = "", cex.main = 1.5, add = T)
abline(v=0, lty=1)
text(x=0.2, y=bp ,labels=tail(lung.go$Description,7), col = "black", xpd=T, cex=1.2, adj=0)



ta.go = go.li$TA@result
ta.go$p.adjust = -log10(ta.go$p.adjust)
colnames(ta.go)[6] = "logFDR"
ta.go = ta.go[order(ta.go$logFDR), ]
ta.go = ta.go[!grepl('positive|negative', ta.go$Description),]
ta.go$Description = str_replace(ta.go$Description, "^(\\w)", toupper)


par(mai = c(0.5,0.1,0.1,1))
bp = barplot(tail(ta.go$logFDR,7), xlim = c(0, 16), horiz = T, xaxt = 'n', yaxt = 'n', 
             xlab = "", names.arg = NA, 
             width = 0.7, border = NA, col = "moccasin", main = "", cex.main = 1.5)
axis(1, at=seq(0,16,4), labels = seq(0,16,4), cex.axis=1, las=1)
abline(v=seq(4,16,4), lty=3, col= "dimgrey")
bp = barplot(tail(ta.go$logFDR,7), xlim = c(0, 16), horiz = T, xaxt = 'n', yaxt = 'n', 
             xlab = "", names.arg = NA, 
             width = 0.7, border = NA, col = "moccasin", main = "", cex.main = 1.5, add = T)
abline(v=0, lty=1)
text(x=0.2, y=bp ,labels=tail(ta.go$Description,7), col = "black", xpd=T, cex=1.2, adj=0)




