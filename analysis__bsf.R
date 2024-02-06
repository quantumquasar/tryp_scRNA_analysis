#analysing the BSF data against the original STF data
#this script is intended to be executed following the original script

library(ggplot2)
library(ComplexHeatmap)
library(EnhancedVolcano)
library(circlize)#colorramp2
library(viridis)#for viridis
library(pheatmap)


#colData_filt, coldata_tmp objects inherited from the original script
colData_bsf <- colData_filt
colData_bsf <- subset(colData_bsf, select = -c(nCell, n.gene, lib.size))
coldata_tmp <- subset(coldata_tmp, select = -c(nCell, n.gene, lib.size))
coldata_tmp$time <- "bsf"
colData_bsf <- rbind(coldata_tmp, colData_bsf)
dim(colData_bsf)
#203 2

h0 <- factor(gsub("bsf", "rest", gsub("7d","rest",gsub("4h","rest",
                                                       gsub("12h","rest",gsub("24h", "rest", gsub("0h", "aim", colData_bsf$time)))))))
h4 <- factor(gsub("bsf", "rest", gsub("7d","rest",gsub("4h","aim",
                                                       gsub("12h","rest",gsub("24h", "rest", gsub("0h", "rest", colData_bsf$time)))))))
h12 <- factor(gsub("bsf", "rest", gsub("7d","rest",gsub("4h","rest",
                                                        gsub("12h","aim",gsub("24h", "rest", gsub("0h", "rest", colData_bsf$time)))))))
h24 <- factor(gsub("bsf", "rest", gsub("7d","rest",gsub("4h","rest",
                                                        gsub("12h","rest",gsub("24h", "aim", gsub("0h", "rest", colData_bsf$time)))))))
d7 <- factor(gsub("bsf", "rest", gsub("7d","aim",gsub("4h","rest",
                                                      gsub("12h","rest",gsub("24h", "rest", gsub("0h", "rest", colData_bsf$time)))))))
hbsf <- factor(gsub("bsf", "aim", gsub("7d","rest",gsub("4h","rest",
                                                      gsub("12h","rest",gsub("24h", "rest", gsub("0h", "rest", colData_bsf$time)))))))

colData_bsf <- data.frame(colData_bsf, h0, h4, h12, h24, d7, hbsf)
colData_bsf$col <- as.character(colData_bsf$col)
rownames(colData_bsf) <- colData_bsf$col
#colData_bsf$batch <- ifelse(startsWith(colData_bsf$col, "bsf"), 1, 2)
#colData_bsf$batch <- as.character(colData_bsf$batch)
#colData_bsf$batch <- factor(colData_bsf$batch)
colData_bsf$time <- factor(colData_bsf$time)


#normalize
sfSingle_bsf <- estimateSizeFactorsForMatrix(combined_cm)
norm_bsf <- t(t(combined_cm/sfSingle_bsf))
norm_bsf <- as.data.frame(norm_bsf)

#####paiwise non-batch#####
#####BSF vs MCFs#####

colData_bsf_pair <- colData_bsf[colData_bsf$time == "bsf"| colData_bsf$time == "0h",]
colData_bsf_pair <- subset(colData_bsf_pair, select = c(col, time))
bsf_cm_pair <- combined_cm[,rownames(colData_bsf_pair)]

#normalise
sfSingle_bsf <- estimateSizeFactorsForMatrix(bsf_cm_pair)
norm_bsf_pair <- t(t(bsf_cm_pair/sfSingle_bsf))
norm_bsf_pair <- as.data.frame(norm_bsf_pair)

dds_bsf_pair <- DESeqDataSetFromMatrix(bsf_cm_pair,
                                       colData = colData_bsf_pair,
                                       design = ~time)
dds_bsf_pair <- DESeq(dds_bsf_pair)

resultsNames(dds_bsf_pair)
res_bsf_pair <- results(dds_bsf_pair, name = "time_bsf_vs_0h")
plotMA(res_bsf_pair)
res_bsf_pair <- data.frame(res_bsf_pair)
res_bsf_pair <- res_bsf_pair[complete.cases(res_bsf_pair),]



#sanity check for signs after DE testing
gene <- norm_bsf_pair["Tb927.10.16560",]
gene<-reshape2::melt(gene)
gene$time <- colData_bsf_pair$time
gene


ggplot((gene), aes(x=time, y=log10(value+1), colour=time))+
  geom_boxplot(outlier.shape = NA) +
  geom_jitter()+
  theme_minimal()+
  #ylim(0, 400)+
  #scale_color_ipsum()+
  ggtitle(paste("Tb927.10.16560"," lfc < 0"))+ #or >
  ylab("normalised expression")



#heatmap
res_bsf_0h_vs_bsf_hm <- res_bsf_pair[res_bsf_pair$padj < 0.01,]
res_bsf_0h_vs_bsf_hm <- res_bsf_0h_vs_bsf_hm[res_bsf_0h_vs_bsf_hm$log2FoldChange < -2 | res_bsf_0h_vs_bsf_hm$log2FoldChange > 2,]

#upregulated
res_bsf_0h_vs_bsf_hm <- res_bsf_0h_vs_bsf_hm[res_bsf_0h_vs_bsf_hm$log2FoldChange < -2, ]
res_bsf_0h_vs_bsf_hm <- res_bsf_0h_vs_bsf_hm[res_bsf_0h_vs_bsf_hm$log2FoldChange > 2, ]


norm <- norm_bsf_pair[rownames(res_bsf_0h_vs_bsf_hm), rownames(colData_bsf_pair)]
scaled <- t(apply(norm, 1, scale))
colnames(scaled) <- colnames(norm)
write.csv(scaled, file = "results/reproduction_mini_bulk/lfc_both_above_below_2_genes.csv")
dim(scaled)

max(scaled)
min(scaled)

p <- Heatmap(scaled, cluster_columns = F, cluster_rows = T,
             col = colorRamp2(seq(-2,2,0.1), viridis(41)),
             row_names_gp = gpar(fontsize = 2),
             show_column_names = TRUE, column_title = "BSF vs MCF: upregulated genes, lfc > 2; padj<0.01")
print(p)


###volcano plots

pdf(file = "smartseq_new/bsf_mapped/plots_de_common_genes/5_deseq_norm/7/pairwise_BSF_vs_MCF.pdf")#,
    #height = 5)

EnhancedVolcano(res_bsf_pair,
                lab = rownames(res_bsf_pair),
                x = 'log2FoldChange',
                y = 'padj',
                title = 'padj<0.01, lfc>2, 0h',
                pCutoff = 0.01,
                FCcutoff = 2)

dev.off()

########BSF vs 7d#####

colData_bsf_pair <- colData_bsf[colData_bsf$time == "bsf"| colData_bsf$time == "7d",]
colData_bsf_pair <- subset(colData_bsf_pair, select = c(col, time))
bsf_cm_pair <- combined_cm[,rownames(colData_bsf_pair)]

#normalise
sfSingle_bsf <- estimateSizeFactorsForMatrix(bsf_cm_pair)
norm_bsf_pair <- t(t(bsf_cm_pair/sfSingle_bsf))
norm_bsf_pair <- as.data.frame(norm_bsf_pair)

dds_bsf_pair <- DESeqDataSetFromMatrix(bsf_cm_pair,
                                       colData = colData_bsf_pair,
                                       design = ~time)
dds_bsf_pair <- DESeq(dds_bsf_pair)

resultsNames(dds_bsf_pair)
res_bsf_pair <- results(dds_bsf_pair, name = "time_bsf_vs_7d")
plotMA(res_bsf_pair)
res_bsf_pair <- data.frame(res_bsf_pair)
res_bsf_pair <- res_bsf_pair[complete.cases(res_bsf_pair),]

#heatmap
res_bsf_7d_vs_bsf_hm <- res_bsf_pair[res_bsf_pair$padj < 0.01,]
res_bsf_7d_vs_bsf_hm <- res_bsf_7d_vs_bsf_hm[res_bsf_7d_vs_bsf_hm$log2FoldChange < -2 | res_bsf_7d_vs_bsf_hm$log2FoldChange > 2,]



norm <- norm_bsf_pair[rownames(res_bsf_7d_vs_bsf_hm), rownames(colData_bsf_pair)]
scaled <- t(apply(norm, 1, scale))
colnames(scaled) <- colnames(norm)
write.csv(scaled, file = "results/reproduction_mini_bulk/lfc_both_above_below_2_genes.csv")
dim(scaled)

max(scaled)
min(scaled)

p <- Heatmap(scaled, cluster_columns = F, cluster_rows = T,
             col = colorRamp2(seq(-2,2,0.1), viridis(41)),
             row_names_gp = gpar(fontsize = 2),
             show_column_names = TRUE, column_title = "BSF vs 7d: absolute lfc > 2; padj<0.01")
print(p)


pdf(file = "smartseq_new/bsf_mapped/plots_de_common_genes/5_deseq_norm/7/pairwise_BSF_vs_7d.pdf")#,
#height = 5)

EnhancedVolcano(res_bsf_pair,
                lab = rownames(res_bsf_pair),
                x = 'log2FoldChange',
                y = 'padj',
                title = 'padj<0.01, lfc>2, BSF vs 7d',
                pCutoff = 0.01,
                FCcutoff = 2)

dev.off()


######BSF vs 24h#####

colData_bsf_pair <- colData_bsf[colData_bsf$time == "bsf"| colData_bsf$time == "24h",]
colData_bsf_pair <- subset(colData_bsf_pair, select = c(col, time))
bsf_cm_pair <- combined_cm[,rownames(colData_bsf_pair)]

#normalise
sfSingle_bsf <- estimateSizeFactorsForMatrix(bsf_cm_pair)
norm_bsf_pair <- t(t(bsf_cm_pair/sfSingle_bsf))
norm_bsf_pair <- as.data.frame(norm_bsf_pair)

dds_bsf_pair <- DESeqDataSetFromMatrix(bsf_cm_pair,
                                       colData = colData_bsf_pair,
                                       design = ~time)
dds_bsf_pair <- DESeq(dds_bsf_pair)

resultsNames(dds_bsf_pair)
res_bsf_pair <- results(dds_bsf_pair, name = "time_bsf_vs_24h")
plotMA(res_bsf_pair)
res_bsf_pair <- data.frame(res_bsf_pair)
res_bsf_pair <- res_bsf_pair[complete.cases(res_bsf_pair),]

#heatmap
res_bsf_24h_vs_bsf_hm <- res_bsf_pair[res_bsf_pair$padj < 0.01,]
res_bsf_24h_vs_bsf_hm <- res_bsf_24h_vs_bsf_hm[res_bsf_24h_vs_bsf_hm$log2FoldChange < -2 | res_bsf_24h_vs_bsf_hm$log2FoldChange > 2,]



norm <- norm_bsf_pair[rownames(res_bsf_24h_vs_bsf_hm), rownames(colData_bsf_pair)]
scaled <- t(apply(norm, 1, scale))
colnames(scaled) <- colnames(norm)
write.csv(scaled, file = "results/reproduction_mini_bulk/lfc_both_above_below_2_genes.csv")
dim(scaled)

max(scaled)
min(scaled)

p <- Heatmap(scaled, cluster_columns = T, cluster_rows = F,
             col = colorRamp2(seq(-2,2,0.1), viridis(41)),
             row_names_gp = gpar(fontsize = 2),
             show_column_names = TRUE, column_title = "BSF vs 24h: absolute lfc > 2; padj<0.01")
print(p)


pdf(file = "smartseq_new/bsf_mapped/plots_de_common_genes/5_deseq_norm/7/pairwise_BSF_vs_24h.pdf")#,
#height = 5)

EnhancedVolcano(res_bsf_pair,
                lab = rownames(res_bsf_pair),
                x = 'log2FoldChange',
                y = 'padj',
                title = 'padj<0.01, lfc>2, 0h',
                pCutoff = 0.01,
                FCcutoff = 2)

dev.off()




#w p value and log fold cutoffs
EnhancedVolcano(res_h4_bsf,
                lab = rownames(res_h4_bsf),
                x = 'log2FoldChange',#
                y = 'padj',
                title = 'padj<0.01, lfc>2, 4h',
                pCutoff = 0.01,
                FCcutoff = 2)







#####DE with common genes between BSF and STFs#####
##no batch correction done

#single_replaced object is taken from the original analysis for the STFs

cm_common_genes <- merge(cm_bsf_unnormalised, single_replaced, by = 0)
dim(cm_common_genes)
#[1] 6875  204
rownames(cm_common_genes) <- cm_common_genes$Row.names
cm_common_genes <- subset(cm_common_genes, select = -c(Row.names))
colnames(cm_common_genes)[1:length(colnames(cm_bsf))] <- paste0("bsf_",1:length(colnames(cm_bsf)))

dim(cm_common_genes)
#[1] 6875  203
dim(cm_common_genes > 0)

sfSingle_bsf <- estimateSizeFactorsForMatrix(cm_common_genes)
norm_bsf_common <- t(t(cm_common_genes/sfSingle_bsf))
norm_bsf_common <- as.data.frame(norm_bsf_common)


#normalising STFs and BSFs separately and merge
#sfSingle <- estimateSizeFactorsForMatrix(single_replaced)
#nSingle <- t(t(single_replaced) / sfSingle)

#norm_bsf_common <- merge(cm_bsf, nSingle, by = 0)
#dim(norm_bsf_common)
#6875 204
#rownames(norm_bsf_common) <- norm_bsf_common$Row.names
#norm_bsf_common <- subset(norm_bsf_common, select = -c(Row.names))
#dim(norm_bsf_common)
#6875 203



#for precomputed size factors using scran::computesumfactors
#norm_diff_common <- assay(norm_common_preclust, "logcounts")
#norm_common_nopreclust
#dim(norm_diff_common)

#####h0####
dds_bsf <- DESeqDataSetFromMatrix(countData = cm_common_genes,
                                  colData = colData_bsf,
                                  design = ~ h0)

#apply size factors already computed using computesumfactors
sizeFactors(dds_bsf) <- sizeFactors(sce2)
head(counts(dds_bsf, normalized = T))


dds_bsf <- DESeq(dds_bsf)

model.matrix(~h0)
resultsNames(dds_bsf)
res_h0_bsf <- results(dds_bsf, name = "h0_rest_vs_aim")
plotMA(res_h0_bsf)
res_h0_bsf <- data.frame(res_h0_bsf)
res_h0_bsf <- res_h0_bsf[complete.cases(res_h0_bsf),]


#aim is the target. if lfc < 0, aim is overexpressed, underexpressed in everything else
#lfc > 0 overexpressed in rest
head(res_h0_bsf[order(res_h0_bsf$log2FoldChange, decreasing = T),])
#lfc > 0
head(res_h0_bsf[order(res_h0_bsf$log2FoldChange, decreasing = F),])
#lfc < 0

norm_bsf_common["AnTat1.1",]
norm_bsf["Tb927.10.11980",]

#####
#sanity check for signs after DE testing
gene <- norm_bsf_common["AnTat1.1",]
gene<-reshape2::melt(gene)
gene$time <- colData_bsf$time
gene


p<-ggplot((gene), aes(x=time, y=log10(value+1), colour=time))+
  geom_boxplot(outlier.shape = NA) +
  geom_jitter()+
  theme_minimal()+
  #ylim(0, 400)+
  #scale_color_ipsum()+
  ggtitle(paste("AnTat1.1"," lfc > 0, h0"))+ #or >
  ylab("normalised expression")
print(p)





#####h4####
dds_bsf <- DESeqDataSetFromMatrix(countData = cm_common_genes,
                                  colData = colData_bsf,
                                  design = ~ h4)

sizeFactors(dds_bsf) <- sizeFactors(sce2)

dds_bsf <- DESeq(dds_bsf)

model.matrix(~h4)
resultsNames(dds_bsf)
res_h4_bsf <- results(dds_bsf, name = "h4_rest_vs_aim")
res_h4_bsf <- data.frame(res_h4_bsf)
res_h4_bsf <- res_h4_bsf[complete.cases(res_h4_bsf),]


results(dds_bsf, contrast = c("time", "aim", "rest"))

#aim is the target. if lfc < 0, aim is overexpressed, underexpressed in everything else
#lfc > 0 overexpressed in rest
head(res_h4_bsf[order(res_h4_bsf$log2FoldChange, decreasing = T),])
#lfc > 0
head(res_h4_bsf[order(res_h4_bsf$log2FoldChange, decreasing = F),])
#lfc < 0

norm_bsf["Tb927.10.1710",]
norm_bsf["Tb927.3.1280",]

#####
#sanity check for signs after DE testing
gene <- norm_bsf["Tb927.3.1280",]
gene<-reshape2::melt(gene)
gene$time <- colData_bsf$time
gene


p<-ggplot((gene), aes(x=time, y=log10(value+1), colour=time))+
  geom_boxplot(outlier.shape = NA) +
  geom_jitter()+
  theme_minimal()+
  #ylim(0, 400)+
  #scale_color_ipsum()+
  ggtitle(paste("Tb927.3.1280"," lfc < 0, h4"))+ #or >
  ylab("normalised expression")
print(p)






#####h12####
dds_bsf <- DESeqDataSetFromMatrix(countData = cm_common_genes,
                                  colData = colData_bsf,
                                  design = ~h12)
sizeFactors(dds_bsf) <- sizeFactors(sce2)

dds_bsf <- DESeq(dds_bsf)

model.matrix(~h12)
resultsNames(dds_bsf)
res_h12_bsf <- results(dds_bsf, name = "h12_rest_vs_aim")
res_h12_bsf <- data.frame(res_h12_bsf)
res_h12_bsf <- res_h12_bsf[complete.cases(res_h12_bsf),]


head(res_h12_bsf[order(res_h12_bsf$log2FoldChange, decreasing = T),])
#lfc > 0, underexpressed
head(res_h12_bsf[order(res_h12_bsf$log2FoldChange, decreasing = F),])
#lfc < 0, overexpressed

norm_bsf["AnTat1.1",]
norm_bsf["Tb927.9.8650",]


#####
#sanity check for signs after DE testing
gene <- norm_bsf["AnTat1.1",]
gene<-reshape2::melt(gene)
gene$time <- colData_bsf$time
gene


p<-ggplot((gene), aes(x=time, y=log10(value+1), colour=time))+
  geom_boxplot(outlier.shape = NA) +
  geom_jitter()+
  theme_minimal()+
  #ylim(0, 400)+
  #scale_color_ipsum()+
  ggtitle(paste("AnTat1.1"," lfc > 0, h12"))+ #or >
  ylab("normalised expression")
print(p)








#####h24####
dds_bsf <- DESeqDataSetFromMatrix(countData = cm_common_genes,
                                  colData = colData_bsf,
                                  design = ~h24)

sizeFactors(dds_bsf) <- sizeFactors(sce2)

dds_bsf <- DESeq(dds_bsf)

model.matrix(~h24)
resultsNames(dds_bsf)
res_h24_bsf <- results(dds_bsf, name = "h24_rest_vs_aim")
res_h24_bsf <- data.frame(res_h24_bsf)
res_h24_bsf <- res_h24_bsf[complete.cases(res_h24_bsf),]

head(res_h24_bsf[order(res_h24_bsf$log2FoldChange, decreasing = T),])
#lfc > 0, under
head(res_h24_bsf[order(res_h24_bsf$log2FoldChange, decreasing = F),])
#lfc < 0, over

norm_bsf["AnTat1.1",]
norm_bsf["Tb927.10.5570",]

#####
#sanity check for signs after DE testing
gene <- norm_bsf["Tb927.10.5570",]

gene<-reshape2::melt(gene)
gene$time <- colData_bsf$time
gene


p<-ggplot((gene), aes(x=time, y=log10(value+1), colour=time))+
  geom_boxplot(outlier.shape = NA) +
  geom_jitter()+
  theme_minimal()+
  #ylim(0, 400)+
  #scale_color_ipsum()+
  ggtitle(paste("Tb927.10.5570"," lfc < 0, h24"))+ #or >
  ylab("normalised expression")
print(p)











#####7d####
dds_bsf <- DESeqDataSetFromMatrix(countData = cm_common_genes,
                                  colData = colData_bsf,
                                  design = ~d7)

sizeFactors(dds_bsf) <- sizeFactors(sce2)

dds_bsf <- DESeq(dds_bsf)

model.matrix(~d7)
resultsNames(dds_bsf)
res_d7_bsf <- results(dds_bsf, name = "d7_rest_vs_aim")
res_d7_bsf <- data.frame(res_d7_bsf)
res_d7_bsf <- res_d7_bsf[complete.cases(res_d7_bsf),]

head(res_d7_bsf[order(res_d7_bsf$log2FoldChange, decreasing = T),])
#lfc > 0
head(res_d7_bsf[order(res_d7_bsf$log2FoldChange, decreasing = F),])
#lfc < 0

norm_bsf["Tb927.10.13240",]
norm_bsf["Tb927.6.1420",]

#####
#sanity check for signs after DE testing
gene <- norm_bsf["Tb927.7.5580",]
gene

gene<-reshape2::melt(gene)
gene$time <- colData_bsf$time
gene


p<-ggplot((gene), aes(x=time, y=log10(value+1), colour=time))+
  geom_boxplot(outlier.shape = NA) +
  geom_jitter()+
  theme_minimal()+
  #ylim(0, 400)+
  #scale_color_ipsum()+
  ggtitle(paste("Tb927.7.5580"," lfc < 0, 7d"))+ #or >
  ylab("normalised expression")
print(p)








#####bsf####

dds_bsf <- DESeqDataSetFromMatrix(countData = cm_common_genes,
                                  colData = colData_bsf,
                                  design = ~hbsf)
sizeFactors(dds_bsf) <- sizeFactors(sce2)

dds_bsf <- DESeq(dds_bsf)

model.matrix(~hbsf)
resultsNames(dds_bsf)
res_hbsf_bsf <- results(dds_bsf, name = "hbsf_rest_vs_aim")
res_hbsf_bsf <- data.frame(res_hbsf_bsf)
res_hbsf_bsf <- res_hbsf_bsf[complete.cases(res_hbsf_bsf),]

head(res_hbsf_bsf[order(res_hbsf_bsf$log2FoldChange, decreasing = T),])
#lfc > 0
head(res_hbsf_bsf[order(res_hbsf_bsf$log2FoldChange, decreasing = F),])
#lfc < 0

norm_bsf["Tb927.3.1150",]
norm_bsf["Tb927.5.950",]

#####
#sanity check for signs after DE testing
gene <- norm_bsf["Tb927.3.1150",]
gene<-reshape2::melt(gene)
gene$time <- colData_bsf$time
gene


p<-ggplot((gene), aes(x=time, y=log10(value+1), colour=time))+
  geom_boxplot(outlier.shape = NA) +
  geom_jitter()+
  theme_minimal()+
  #ylim(0, 400)+
  #scale_color_ipsum()+
  ggtitle(paste("Tb927.3.1150"," lfc > 0, bsf"))+ #or >
  ylab("normalised expression")
print(p)






#####heatmaps####
#####absolute logfc 2 OR 1####
res_h0_hm_bsf <- res_h0_bsf[res_h0_bsf$padj < 0.01,]

#write.csv(res_h0_hm_bsf[res_h0_hm_bsf$log2FoldChange < -2,], "plots_de_common_genes/bsf_vs_all_p0.05_overexpressed_h0.csv")
#write.csv(res_h0_hm_bsf[res_h0_hm_bsf$log2FoldChange > 2,], "plots_de_common_genes/1_new_norm/bsf_vs_all_p0.05_underexpressed_h0.csv")

res_h0_hm_bsf <- res_h0_hm_bsf[res_h0_hm_bsf$log2FoldChange < -2 | res_h0_hm_bsf$log2FoldChange > 2,]

res_h0_hm_bsf <- res_h0_hm_bsf[res_h0_hm_bsf$log2FoldChange < -2,]


res_h24_hm_bsf <- res_h24_bsf[res_h24_bsf$padj < 0.01,]

#write.csv(res_h24_hm_bsf[res_h24_hm_bsf$log2FoldChange < -1,], "plots_de_common_genes/1_new_norm/bsf_vs_all_p0.05_overexpressed_h24.csv")
#write.csv(res_h24_hm_bsf[res_h24_hm_bsf$log2FoldChange > 1,], "plots_de_common_genes/1_new_norm/bsf_vs_all_p0.05_underexpressed_h24.csv")


res_h24_hm_bsf <- res_h24_hm_bsf[res_h24_hm_bsf$log2FoldChange < -2 | res_h24_hm_bsf$log2FoldChange > 2,]
res_h24_hm_bsf <- res_h24_hm_bsf[res_h24_hm_bsf$log2FoldChange < -2,]


res_d7_hm_bsf <- res_d7_bsf[res_d7_bsf$padj < 0.01,]
#write.csv(res_d7_hm_bsf[res_d7_hm_bsf$log2FoldChange < -1,],
#          "plots_de_common_genes/1_new_norm/expressed_genes_lfc_2/bsf_vs_all_p0.05_overexpressed_d7.csv")
#write.csv(res_d7_hm_bsf[res_d7_hm_bsf$log2FoldChange > 1,],
#          "plots_de_common_genes/1_new_norm/expressed_genes_lfc_2/bsf_vs_all_p0.05_underexpressed_d7.csv")


res_d7_hm_bsf <- res_d7_hm_bsf[res_d7_hm_bsf$log2FoldChange < -2 | res_d7_hm_bsf$log2FoldChange > 2,]
res_d7_hm_bsf <- res_d7_hm_bsf[res_d7_hm_bsf$log2FoldChange < -2,]


res_h4_hm_bsf <- res_h4_bsf[res_h4_bsf$padj < 0.01,]
#write.csv(res_h4_hm_bsf[res_h4_hm_bsf$log2FoldChange < -1,],
#          "plots_de_common_genes/1_new_norm/expressed_genes_lfc_2/bsf_vs_all_p0.05_overexpressed_h4.csv")
#write.csv(res_h4_hm_bsf[res_h4_hm_bsf$log2FoldChange > 1,],
#          "plots_de_common_genes/1_new_norm/expressed_genes_lfc_2/bsf_vs_all_p0.05_underexpressed_h4.csv")


res_h4_hm_bsf <- res_h4_hm_bsf[res_h4_hm_bsf$log2FoldChange < -2 | res_h4_hm_bsf$log2FoldChange > 2,]
res_h4_hm_bsf <- res_h4_hm_bsf[res_h4_hm_bsf$log2FoldChange < -2,]


res_h12_hm_bsf <- res_h12_bsf[res_h12_bsf$padj < 0.01,]
#write.csv(res_h12_hm_bsf[res_h12_hm_bsf$log2FoldChange < -1,],
#          "plots_de_common_genes/1_new_norm/expressed_genes_lfc_2/bsf_vs_all_p0.05_overexpressed_h12.csv")
#write.csv(res_h12_hm_bsf[res_h12_hm_bsf$log2FoldChange > 1,],
#          "plots_de_common_genes/1_new_norm/expressed_genes_lfc_2/bsf_vs_all_p0.05_underexpressed_h12.csv")


res_h12_hm_bsf <- res_h12_hm_bsf[res_h12_hm_bsf$log2FoldChange < -2 | res_h12_hm_bsf$log2FoldChange > 2,]
res_h12_hm_bsf <- res_h12_hm_bsf[res_h12_hm_bsf$log2FoldChange < -2,]


res_hbsf_hm_bsf <- res_hbsf_bsf[res_hbsf_bsf$padj < 0.01,]
#write.csv(res_hbsf_hm_bsf[res_hbsf_hm_bsf$log2FoldChange < -1,],
#          "plots_de_common_genes/1_new_norm/expressed_genes_lfc_2/bsf_vs_all_p0.05_overexpressed_hbsf.csv")
#write.csv(res_hbsf_hm_bsf[res_hbsf_hm_bsf$log2FoldChange > 1,],
#          "plots_de_common_genes/1_new_norm/expressed_genes_lfc_2/bsf_vs_all_p0.05_underexpressed_hbsf.csv")


res_hbsf_hm_bsf <- res_hbsf_hm_bsf[res_hbsf_hm_bsf$log2FoldChange < -2 | res_hbsf_hm_bsf$log2FoldChange > 2,]
res_hbsf_hm_bsf <- res_hbsf_hm_bsf[res_hbsf_hm_bsf$log2FoldChange < -2,]


hm_gene_bsf <- unique(c(rownames(res_h0_hm_bsf), rownames(res_h24_hm_bsf), rownames(res_d7_hm_bsf),
                        rownames(res_h4_hm_bsf), rownames(res_h12_hm_bsf), rownames(res_hbsf_hm_bsf)))

length(unique(c(rownames(res_hbsf_hm_bsf), rownames(res_h0_hm_bsf), rownames(res_h24_hm_bsf), rownames(res_d7_hm_bsf),
         rownames(res_h4_hm_bsf), rownames(res_h12_hm_bsf))))


#####p 0.05 expressed genes####
dim(cm_bsf_unnormalised)
gene_exp <- (cm_common_genes[hm_gene_bsf,colData_bsf[colData_bsf$time == "0h", ]$col])
gene_exp <- gene_exp[rowSums(gene_exp) > 0,]


gene_exp <- (cm_common_genes[hm_gene_bsf,colData_bsf[colData_bsf$time == "bsf", ]$col])
gene_exp <- gene_exp[rowSums(gene_exp) > 0,]


gene_exp <- (cm_common_genes[hm_gene_bsf,colData_bsf[colData_bsf$time == "12h", ]$col])
gene_exp <- gene_exp[rowSums(gene_exp) > 0,]


gene_exp <- (cm_common_genes[hm_gene_bsf,colData_bsf[colData_bsf$time == "24h", ]$col])
gene_exp <- gene_exp[rowSums(gene_exp) > 0,]


gene_exp <- (cm_common_genes[hm_gene_bsf,colData_bsf[colData_bsf$time == "7d", ]$col])
gene_exp <- gene_exp[rowSums(gene_exp) > 0,]




#####p 0.01 expressed genes####
gene_exp <- (cm_common_genes[hm_gene_bsf,colData_bsf[colData_bsf$time == "0h", ]$col])
gene_exp <- gene_exp[rowSums(gene_exp) > 0,]

gene_exp <- (cm_common_genes[hm_gene_bsf,colData_bsf[colData_bsf$time == "bsf", ]$col])
gene_exp <- gene_exp[rowSums(gene_exp) > 0,]

gene_exp <- (cm_common_genes[hm_gene_bsf,colData_bsf[colData_bsf$time == "12h", ]$col])
gene_exp <- gene_exp[rowSums(gene_exp) > 0,]

gene_exp <- (cm_common_genes[hm_gene_bsf,colData_bsf[colData_bsf$time == "24h", ]$col])
gene_exp <- gene_exp[rowSums(gene_exp) > 0,]


gene_exp <- (cm_common_genes[hm_gene_bsf,colData_bsf[colData_bsf$time == "7d", ]$col])
gene_exp <- gene_exp[rowSums(gene_exp) > 0,]


#####
### creating the heatmap
h24_bsf <- as.character(colData_bsf[colData_bsf$time == "24h", ]$col)
h0_bsf <- as.character(colData_bsf[colData_bsf$time == "0h", ]$col)
d7_bsf <- as.character(colData_bsf[colData_bsf$time == "7d", ]$col)
h4_bsf <- as.character(colData_bsf[colData_bsf$time == "4h",]$col)
h12_bsf <- as.character(colData_bsf[colData_bsf$time == "12h",]$col)
hbsf_bsf <- as.character(colData_bsf[colData_bsf$time == "bsf",]$col)
ht_order <- c(h0_bsf, h4_bsf, h12_bsf,  h24_bsf, d7_bsf, hbsf_bsf)

norm <- norm_bsf_common[hm_gene_bsf, ht_order]

scaled <- t(apply(norm, 1, scale))
colnames(scaled) <- colnames(norm)

dim(scaled)

max(scaled)
min(scaled)

p <- Heatmap(scaled, cluster_columns = F, cluster_rows = T,
             col = colorRamp2(seq(-2,2,0.1), viridis(41)),
             row_names_gp = gpar(fontsize = 2),
             show_column_names = TRUE, column_title = "BSF vs all timepoints, common genes:
             abs lfc > 2; padj<0.01")
print(p)

Heatmap(scaled, cluster_columns = F, cluster_rows = T,
             col = colorRamp2(seq(-2,2,0.1), viridis(41)),
             row_names_gp = gpar(fontsize = 2),
             show_column_names = TRUE, column_title = "BSF vs all timepoints, common genes:
              lfc < -2; padj<0.01")


##finding common upregulated genes between stf+bsf data and old STF data
length(intersect(hm_gene, hm_gene_bsf))

norm <- norm_bsf_common[intersect(hm_gene, hm_gene_bsf), ht_order]
scaled <- t(apply(norm, 1, scale))
colnames(scaled) <- colnames(norm)

dim(scaled)

max(scaled)
min(scaled)
Heatmap(scaled, cluster_columns = F, cluster_rows = T,
        col = colorRamp2(seq(-2,2,0.1), viridis(41)),
        row_names_gp = gpar(fontsize = 2),
        show_column_names = TRUE, column_title = "common genes between orig STF and incl. BSF data:
              abs lfc > 2; padj<0.01")

########volcano plots non-batch##########
####+++++++++++++###
####

EnhancedVolcano(res_h0_bsf,
                lab = rownames(res_h0_bsf),
                x = 'log2FoldChange',
                y = 'padj',
                title = 'padj<0.01, lfc>2, 0h',
                pCutoff = 0.01,
                FCcutoff = 2)


#w p value and log fold cutoffs
EnhancedVolcano(res_h4_bsf,
                lab = rownames(res_h4_bsf),
                x = 'log2FoldChange',#
                y = 'padj',
                title = 'padj<0.01, lfc>2, 4h',
                pCutoff = 0.01,
                FCcutoff = 2)
#pointSize = 3.0,
#labSize = 6.0)

EnhancedVolcano(res_h24_bsf,
                lab = rownames(res_h24_bsf),
                x = 'log2FoldChange',#
                y = 'padj',
                title = 'padj<0.01, lfc>2, 24h',
                pCutoff = 0.01,
                FCcutoff = 2)

EnhancedVolcano(res_d7_bsf,
                lab = rownames(res_d7_bsf),
                x = 'log2FoldChange',#
                y = 'padj',
                title = 'padj<0.01, lfc>2, 7d',
                pCutoff = 0.01,
                FCcutoff = 2)

EnhancedVolcano(res_h12_bsf,
                lab = rownames(res_h12_bsf),
                x = 'log2FoldChange',#
                y = 'padj',
                title = 'padj<0.01, lfc>2, 12h',
                pCutoff = 0.01,
                FCcutoff = 2)

EnhancedVolcano(res_hbsf_bsf,
                lab = rownames(res_hbsf_bsf),
                x = 'log2FoldChange',#
                y = 'padj',
                title = 'padj<0.01, lfc>2, bsf',
                pCutoff = 0.01,
                FCcutoff = 2)


#####PCA + cor####

##this is normalised using computesumfactors (norm method for sc specifically)
#norm_diff_common <- assay(norm_common_preclust, "logcounts")
#dim(norm_diff_common)
#1] 6875  203
#no preclustering
#norm_diff_common <- assay(norm_common_nopreclust, "logcounts")

#topVarGenesSingle <- head(order(rowVars(as.matrix(norm_diff_common)), decreasing = TRUE),1000) #1000 in manuscript
topVarGenesSingle <- head(order(rowVars(as.matrix(norm_bsf_common)), decreasing = TRUE),1000) #1000 in manuscript
#nSingle_fin <- log(norm_diff_common[topVarGenesSingle,] + 1)
#nSingle_fin <- (norm_diff_common[topVarGenesSingle,])
nSingle_fin <- norm_bsf_common[topVarGenesSingle,]
### Performing PCA (single cell)
res.pca <- prcomp(t(nSingle_fin), scale = FALSE)

variance_prop <- round(as.data.frame(summary(res.pca)[[6]])[2,],3)
variance_prop <- c(variance_prop[1], variance_prop[2])
variance_prop[[1]]#0.253
variance_prop[[2]]#0.089

#pdf(file = "smartseq_new/bsf_mapped/plots_de_common_genes/5_deseq_norm/pca_w_BSF_corrected.pdf")

ggplot(data.frame(PC1 = res.pca$x[,"PC1"], PC2 = res.pca$x[,"PC2"]),
       #time = factor(colData_filt$time,levels = c("0h", "4h", "12h", "24h", "7d"))),
       aes(x = PC1, y = PC2, color = colData_bsf$time)) +
  labs(
    x = paste("PC1 (",variance_prop[[1]]*100, "% )", sep = ""),#c("PC1",variance_prop[[1]]),
    y = paste("PC2 (",variance_prop[[2]]*100, "% )", sep = ""),
    color="time"
  )+
  coord_cartesian(xlim = c(-80, 80))+
  geom_jitter(size = 2) +
  theme_classic() +
  theme() +
  geom_hline(yintercept=0, linetype="dashed") +
  geom_vline(xintercept=0, linetype="dashed")

#dev.off()

##correlation


#calculating the correlation
#pearson on log-transformed read counts
cor_bsf_common <- cor(log(norm_bsf_common + 1), method = "pearson")
cor_bsf_common <- cor((norm_bsf_common), method = "pearson")

Heatmap(cor_bsf_common, cluster_columns = F, cluster_rows = F, #row_dend_side = "right",
        #col = colorRamp2(seq(0.97,1,0.005), viridis(41)),
        col = colorRamp2(seq(min(cor_bsf), seq(max(cor_bsf)), length = 3), c("light blue", "dark blue", "yellow")),
        #col = rev(rainbow(10)),
        row_names_gp = gpar(fontsize = 4),
        show_column_names = TRUE, column_title = "Pearson correlation on BSF + STFs (common genes)",
        clustering_method_rows = "single")

library(corrplot)

corrplot(cor_rld_SS, method = "number", type = "full")


###z count##
top_z <- head(order(as.matrix(norm_bsf_log), decreasing = TRUE),1000)
z_count <- t(scale(t(top_z), center=T, scale=T))#for the genes with highest expression

z_count <- t(scale(t(as.matrix(log(norm_diff_common + 1))), center=T, scale=T))
z_count <- t(scale(t(as.matrix((norm_diff_common))), center=T, scale=T))


Heatmap(z_count, cluster_columns = T, cluster_rows = T, #row_dend_side = "right",
        #col = colorRamp2(seq(-2,2,0.1), viridis(41)),
        col = rev(rainbow(10)),
        row_names_gp = gpar(fontsize = 2),
        show_column_names = TRUE, column_title = "z-score BSF vs  STFs (common genes)",
        clustering_method_rows = "single")





#####pairwise common genes#####
#####BSF vs MCFs#####

colData_bsf_pair <- colData_bsf[colData_bsf$time == "bsf"| colData_bsf$time == "0h",]
colData_bsf_pair <- subset(colData_bsf_pair, select = c(col, time))
bsf_cm_pair <- cm_common_genes[,rownames(colData_bsf_pair)]

#normalise
sfSingle_bsf <- estimateSizeFactorsForMatrix(bsf_cm_pair)
norm_bsf_pair <- t(t(bsf_cm_pair/sfSingle_bsf))
norm_bsf_pair <- as.data.frame(norm_bsf_pair)


dds_bsf_pair <- DESeqDataSetFromMatrix(bsf_cm_pair,
                                       colData = colData_bsf_pair,
                                       design = ~time)

sizeFactors(dds_bsf_pair) <- sizeFactors(sce_pairs)

dds_bsf_pair <- DESeq(dds_bsf_pair)

resultsNames(dds_bsf_pair)
res_bsf_pair <- results(dds_bsf_pair, name = "time_bsf_vs_0h")
plotMA(res_bsf_pair)
res_bsf_pair <- data.frame(res_bsf_pair)
res_bsf_pair <- res_bsf_pair[complete.cases(res_bsf_pair),]


head(res_bsf_pair[order(res_bsf_pair$log2FoldChange, decreasing = T),])
#lfc > 0
head(res_bsf_pair[order(res_bsf_pair$log2FoldChange, decreasing = F),])
#lfc < 0


#####
#sanity check for signs after DE testing
gene <- norm_bsf_pair["AnTat1.1",]
gene<-reshape2::melt(gene)
gene$time <- colData_bsf_pair$time
gene


p<-ggplot((gene), aes(x=time, y=log10(value+1), colour=time))+
  geom_boxplot(outlier.shape = NA) +
  geom_jitter()+
  theme_minimal()+
  #ylim(0, 400)+
  #scale_color_ipsum()+
  ggtitle(paste("Tb927.11.4750"," lfc > 0, bsf vs mcf"))+ #or >
  ylab("normalised expression")
print(p)



#heatmap
res_bsf_0h_vs_bsf_hm <- res_bsf_pair[res_bsf_pair$padj < 0.01,]
#write.csv(res_bsf_0h_vs_bsf_hm[res_bsf_0h_vs_bsf_hm$log2FoldChange < -2,],
#          "smartseq_new/bsf_mapped/plots_de_common_genes/5_deseq_norm/7/pairwise_bsf_vs_mcf_up_mcf_down_bsf.csv")
          #"plots_de_common_genes/1_new_norm/expressed_genes_lfc_2/pairwise/bsf_vs_mcf_overexpressed_mcf.csv")
#write.csv(res_bsf_0h_vs_bsf_hm[res_bsf_0h_vs_bsf_hm$log2FoldChange > 2,],
#          "smartseq_new/bsf_mapped/plots_de_common_genes/5_deseq_norm/7/pairwise_bsf_vs_mcf_down_mcf_up_bs.csv")

bsf_vs_0h_overxp <- rownames(res_bsf_0h_vs_bsf_hm[res_bsf_0h_vs_bsf_hm$log2FoldChange < -2,])
bsf_vs_0h_underxp <- rownames(res_bsf_0h_vs_bsf_hm[res_bsf_0h_vs_bsf_hm$log2FoldChange > 2,])
length(intersect(bsf_vs_0h_overxp, bsf_vs_0h_underxp))


res_bsf_0h_vs_bsf_hm <- res_bsf_0h_vs_bsf_hm[res_bsf_0h_vs_bsf_hm$log2FoldChange < -2 | res_bsf_0h_vs_bsf_hm$log2FoldChange > 2,]
res_bsf_0h_vs_bsf_hm <- res_bsf_0h_vs_bsf_hm[res_bsf_0h_vs_bsf_hm$log2FoldChange < -2,]
res_bsf_0h_vs_bsf_hm <- res_bsf_0h_vs_bsf_hm[res_bsf_0h_vs_bsf_hm$log2FoldChange > 2,]



norm <- norm_bsf_pair[rownames(res_bsf_0h_vs_bsf_hm), rownames(colData_bsf_pair)]
scaled <- t(apply(norm, 1, scale))
colnames(scaled) <- colnames(norm)
#write.csv(scaled, file = "results/reproduction_mini_bulk/lfc_both_above_below_2_genes.csv")
dim(scaled)

max(scaled)
min(scaled)

p <- Heatmap(scaled, cluster_columns = F, cluster_rows = F,
             col = colorRamp2(seq(-2,2,0.1), viridis(41)),
             row_names_gp = gpar(fontsize = 2),
             show_column_names = TRUE, column_title = "BSF vs MCF, common genes (diff norm):
             abs lfc > 2; padj<0.01")
print(p)


###volcano plots

#pdf(file = "smartseq_new/bsf_mapped/plots_de_common_genes/5_deseq_norm/7/pairwise_BSF_vs_MCF.pdf")

EnhancedVolcano(res_bsf_pair,
                lab = rownames(res_bsf_pair),
                x = 'log2FoldChange',
                y = 'padj',
                title = 'padj<0.01, lfc>2, bsf_vs_mcf',
                pCutoff = 0.01,
                FCcutoff = 2)

#dev.off()

########BSF vs 7d#####

colData_bsf_pair <- colData_bsf[colData_bsf$time == "bsf"| colData_bsf$time == "7d",]
colData_bsf_pair <- subset(colData_bsf_pair, select = c(col, time))
bsf_cm_pair <- cm_common_genes[,rownames(colData_bsf_pair)]

#normalise
sfSingle_bsf <- estimateSizeFactorsForMatrix(bsf_cm_pair)
norm_bsf_pair <- t(t(bsf_cm_pair/sfSingle_bsf))
norm_bsf_pair <- as.data.frame(norm_bsf_pair)


dds_bsf_pair <- DESeqDataSetFromMatrix(bsf_cm_pair,
                                       colData = colData_bsf_pair,
                                       design = ~time)
sizeFactors(dds_bsf_pair) <- sizeFactors(sce_pairs)

dds_bsf_pair <- DESeq(dds_bsf_pair)

resultsNames(dds_bsf_pair)
res_bsf_pair <- results(dds_bsf_pair, name = "time_bsf_vs_7d")
plotMA(res_bsf_pair)
res_bsf_pair <- data.frame(res_bsf_pair)
res_bsf_pair <- res_bsf_pair[complete.cases(res_bsf_pair),]

#heatmap
res_bsf_7d_vs_bsf_hm <- res_bsf_pair[res_bsf_pair$padj < 0.01,]#2089
#write.csv(res_bsf_7d_vs_bsf_hm[res_bsf_7d_vs_bsf_hm$log2FoldChange < -2,],
#          "smartseq_new/bsf_mapped/plots_de_common_genes/5_deseq_norm/7/pairwise_bsf_vs_7d_up_7d_down_bsf.csv")
#write.csv(res_bsf_7d_vs_bsf_hm[res_bsf_7d_vs_bsf_hm$log2FoldChange > 2,],
#          "smartseq_new/bsf_mapped/plots_de_common_genes/5_deseq_norm/7/pairwise_bsf_vs_7d_down_7d_up_bsf.csv")


bsf_vs_0h_overxp <- (rownames(res_bsf_7d_vs_bsf_hm[res_bsf_7d_vs_bsf_hm$log2FoldChange < -2,]))
bsf_vs_0h_underxp <- rownames(res_bsf_7d_vs_bsf_hm[res_bsf_7d_vs_bsf_hm$log2FoldChange > 2,])
length(intersect(bsf_vs_0h_overxp, bsf_vs_0h_underxp))

res_bsf_7d_vs_bsf_hm <- res_bsf_7d_vs_bsf_hm[res_bsf_7d_vs_bsf_hm$log2FoldChange > 2 | res_bsf_7d_vs_bsf_hm$log2FoldChange < -2,]

res_bsf_7d_vs_bsf_hm <- res_bsf_7d_vs_bsf_hm[res_bsf_7d_vs_bsf_hm$log2FoldChange > 2,]
res_bsf_7d_vs_bsf_hm <- res_bsf_7d_vs_bsf_hm[res_bsf_7d_vs_bsf_hm$log2FoldChange < -2,]



norm <- norm_bsf_pair[rownames(res_bsf_7d_vs_bsf_hm), rownames(colData_bsf_pair)]
scaled <- t(apply(norm, 1, scale))
colnames(scaled) <- colnames(norm)
#write.csv(scaled, file = "results/reproduction_mini_bulk/lfc_both_above_below_2_genes.csv")
dim(scaled)

max(scaled)
min(scaled)

p <- Heatmap(scaled, cluster_columns = F, cluster_rows = F,
             col = colorRamp2(seq(-2,2,0.1), viridis(41)),
             row_names_gp = gpar(fontsize = 2),
             show_column_names = TRUE, column_title = "BSF vs 7d (common genes), diff norm:
             abs lfc > 2; padj<0.01")
print(p)


#pdf(file = "smartseq_new/bsf_mapped/plots_de_common_genes/5_deseq_norm/7/pairwise_BSF_vs_7d.pdf")

EnhancedVolcano(res_bsf_pair,
                lab = rownames(res_bsf_pair),
                x = 'log2FoldChange',
                y = 'padj',
                title = 'padj<0.01, lfc>2, BSF vs 7d',
                pCutoff = 0.01,
                FCcutoff = 2)

#dev.off()


`######BSF vs 24h#####

colData_bsf_pair <- colData_bsf[colData_bsf$time == "bsf"| colData_bsf$time == "24h",]
colData_bsf_pair <- subset(colData_bsf_pair, select = c(col, time))
bsf_cm_pair <- cm_common_genes[,rownames(colData_bsf_pair)]

#normalise
sfSingle_bsf <- estimateSizeFactorsForMatrix(bsf_cm_pair)
norm_bsf_pair <- t(t(bsf_cm_pair/sfSingle_bsf))
norm_bsf_pair <- as.data.frame(norm_bsf_pair)



dds_bsf_pair <- DESeqDataSetFromMatrix(bsf_cm_pair,
                                       colData = colData_bsf_pair,
                                       design = ~time)
sizeFactors(dds_bsf_pair) <- sizeFactors(sce_pairs)

dds_bsf_pair <- DESeq(dds_bsf_pair)

resultsNames(dds_bsf_pair)
res_bsf_pair <- results(dds_bsf_pair, name = "time_bsf_vs_24h")
plotMA(res_bsf_pair)
res_bsf_pair <- data.frame(res_bsf_pair)
res_bsf_pair <- res_bsf_pair[complete.cases(res_bsf_pair),]

#heatmap
res_bsf_24h_vs_bsf_hm <- res_bsf_pair[res_bsf_pair$padj < 0.01,]
#write.csv(res_bsf_24h_vs_bsf_hm[res_bsf_24h_vs_bsf_hm$log2FoldChange < -2,],
#          "smartseq_new/bsf_mapped/plots_de_common_genes/5_deseq_norm/7/pairwise_bsf_vs_24h_up_24h_down_bsf.csv")
#write.csv(res_bsf_24h_vs_bsf_hm[res_bsf_24h_vs_bsf_hm$log2FoldChange > 2,],
#          "smartseq_new/bsf_mapped/plots_de_common_genes/5_deseq_norm/7/pairwise_bsf_vs_24h_down_24h_up_bsf.csv")



bsf_vs_0h_overxp <- (rownames(res_bsf_24h_vs_bsf_hm[res_bsf_24h_vs_bsf_hm$log2FoldChange < -2,]))
bsf_vs_0h_underxp <- rownames(res_bsf_24h_vs_bsf_hm[res_bsf_24h_vs_bsf_hm$log2FoldChange > 2,])
length(intersect(bsf_vs_0h_overxp, bsf_vs_0h_underxp))



res_bsf_24h_vs_bsf_hm <- res_bsf_24h_vs_bsf_hm[res_bsf_24h_vs_bsf_hm$log2FoldChange < -2 | res_bsf_24h_vs_bsf_hm$log2FoldChange > 2,]



norm <- norm_bsf_pair[rownames(res_bsf_24h_vs_bsf_hm), rownames(colData_bsf_pair)]
scaled <- t(apply(norm, 1, scale))
colnames(scaled) <- colnames(norm)
#write.csv(scaled, file = "results/reproduction_mini_bulk/lfc_both_above_below_2_genes.csv")
dim(scaled)

max(scaled)
min(scaled)

p <- Heatmap(scaled, cluster_columns = F, cluster_rows = F,
             col = colorRamp2(seq(-2,2,0.1), viridis(41)),
             row_names_gp = gpar(fontsize = 2),
             show_column_names = TRUE, column_title = "BSF vs 24h (common genes), diff norm:
             absolute lfc > 2; padj<0.01")
print(p)


#pdf(file = "smartseq_new/bsf_mapped/plots_de_common_genes/5_deseq_norm/7/pairwise_BSF_vs_24h.pdf")

EnhancedVolcano(res_bsf_pair,
                lab = rownames(res_bsf_pair),
                x = 'log2FoldChange',
                y = 'padj',
                title = 'padj<0.01, lfc>2, BSF vs 24h',
                pCutoff = 0.01,
                FCcutoff = 2)

#dev.off()




#w p value and log fold cutoffs
EnhancedVolcano(res_bsf_pair,
                lab = rownames(res_bsf_pair),
                x = 'log2FoldChange',#
                y = 'padj',
                title = 'padj<0.01, lfc>2, bsf_vs_24h',
                pCutoff = 0.01,
                FCcutoff = 2)




#####BSF vs MCF vs 7d####

colData_bsf_pair <- colData_bsf[colData_bsf$time == "bsf"| colData_bsf$time == "0h" | colData_bsf$time == "7d",]
colData_bsf_pair <- subset(colData_bsf_pair, select = c(col, time))
bsf_cm_pair <- cm_common_genes[,rownames(colData_bsf_pair)]

#normalise
sfSingle_bsf <- estimateSizeFactorsForMatrix(bsf_cm_pair)
norm_bsf_pair <- t(t(bsf_cm_pair/sfSingle_bsf))
norm_bsf_pair <- as.data.frame(norm_bsf_pair)

h0 <- factor(gsub("bsf", "rest", gsub("7d","rest", gsub("0h", "aim", colData_bsf_pair$time))))
d7 <- factor(gsub("bsf", "rest", gsub("7d","aim", gsub("0h", "rest", colData_bsf_pair$time))))
hbsf <- factor(gsub("bsf", "aim", gsub("7d","rest", gsub("0h", "rest", colData_bsf_pair$time))))

colData_bsf_pair <- data.frame(colData_bsf_pair, h0, d7, hbsf)
colData_bsf_pair$col <- as.character(colData_bsf_pair$col)
rownames(colData_bsf_pair) <- colData_bsf_pair$col
colData_bsf$time <- factor(colData_bsf$time)


#mcf
dds_bsf_pair <- DESeqDataSetFromMatrix(bsf_cm_pair,
                                       colData = colData_bsf_pair,
                                       design = ~h0)
sizeFactors(dds_bsf_pair) <- sizeFactors(sce_pairs)
dds_bsf_pair <- DESeq(dds_bsf_pair)

resultsNames(dds_bsf_pair)
res_bsf_vs_0h <- results(dds_bsf_pair, name = "h0_rest_vs_aim")
plotMA(res_bsf_vs_0h)
res_bsf_vs_0h <- data.frame(res_bsf_vs_0h)
res_bsf_vs_0h <- res_bsf_vs_0h[complete.cases(res_bsf_vs_0h),]


head(res_bsf_vs_0h[order(res_bsf_vs_0h$log2FoldChange, decreasing = T),])
#lfc > 0
head(res_bsf_vs_0h[order(res_bsf_vs_0h$log2FoldChange, decreasing = F),])
#lfc < 0


#####
#sanity check for signs after DE testing
gene <- norm_bsf_pair["Tb927.11.4750",]
gene<-reshape2::melt(gene)
gene$time <- colData_bsf_pair$time
gene


p<-ggplot((gene), aes(x=time, y=log10(value+1), colour=time))+
  geom_boxplot(outlier.shape = NA) +
  geom_jitter()+
  theme_minimal()+
  #ylim(0, 400)+
  #scale_color_ipsum()+
  ggtitle(paste("Tb927.10.4080"," lfc > 0, bsf vs mcf"))+ #or >
  ylab("normalised expression")
print(p)

#lfc >0 underexpressed in tested time point; 0h
#lfc <0 overexpressed in tested time point; 0h


#7d
dds_bsf_pair <- DESeqDataSetFromMatrix(bsf_cm_pair,
                                       colData = colData_bsf_pair,
                                       design = ~d7)
sizeFactors(dds_bsf_pair) <- sizeFactors(sce_pairs)
dds_bsf_pair <- DESeq(dds_bsf_pair)

resultsNames(dds_bsf_pair)
res_bsf_vs_7d <- results(dds_bsf_pair, name = "d7_rest_vs_aim")

res_bsf_vs_7d <- data.frame(res_bsf_vs_7d)
res_bsf_vs_7d <- res_bsf_vs_7d[complete.cases(res_bsf_vs_7d),]



#bsf
dds_bsf_pair <- DESeqDataSetFromMatrix(bsf_cm_pair,
                                       colData = colData_bsf_pair,
                                       design = ~hbsf)
sizeFactors(dds_bsf_pair) <- sizeFactors(sce_pairs)
dds_bsf_pair <- DESeq(dds_bsf_pair)

resultsNames(dds_bsf_pair)
res_bsf_vs_all <- results(dds_bsf_pair, name = "hbsf_rest_vs_aim")
plotMA(res_bsf_vs_all)
res_bsf_vs_all <- data.frame(res_bsf_vs_all)
res_bsf_vs_all <- res_bsf_vs_all[complete.cases(res_bsf_vs_all),]


res_bsf_0h_hm <- res_bsf_vs_0h[res_bsf_vs_0h$padj<0.01,]
#write.csv(res_bsf_0h_hm[res_bsf_0h_hm$log2FoldChange < -2,],
#          "smartseq_new/bsf_mapped/plots_de_common_genes/5_deseq_norm/7/pairwise_bsf_vs_all_up_mcf_down_rest.csv")
#write.csv(res_bsf_0h_hm[res_bsf_0h_hm$log2FoldChange > 2,],
#          "smartseq_new/bsf_mapped/plots_de_common_genes/5_deseq_norm/7/pairwise_bsf_vs_all_down_mcf_up_rest.csv")
res_bsf_0h_hm <- res_bsf_0h_hm[res_bsf_0h_hm$log2FoldChange < -2 | res_bsf_0h_hm$log2FoldChange > 2,]

res_bsf_7d_hm <- res_bsf_vs_7d[res_bsf_vs_7d$padj<0.01,]
#write.csv(res_bsf_7d_hm[res_bsf_7d_hm$log2FoldChange < -2,],
#          "smartseq_new/bsf_mapped/plots_de_common_genes/5_deseq_norm/7/pairwise_bsf_vs_all_up_7d_down_rest.csv")
#write.csv(res_bsf_7d_hm[res_bsf_7d_hm$log2FoldChange > 2,],
#          "smartseq_new/bsf_mapped/plots_de_common_genes/5_deseq_norm/7/pairwise_bsf_vs_all_down_7d_up_rest.csv")
res_bsf_7d_hm <- res_bsf_7d_hm[res_bsf_7d_hm$log2FoldChange < -2 | res_bsf_7d_hm$log2FoldChange > 2,]

res_bsf_all_hm <- res_bsf_vs_all[res_bsf_vs_all$padj<0.01,]
#write.csv(res_bsf_all_hm[res_bsf_all_hm$log2FoldChange < -2,],
#          "smartseq_new/bsf_mapped/plots_de_common_genes/5_deseq_norm/7/pairwise_bsf_vs_all_up_bsf_down_rest.csv")
#write.csv(res_bsf_all_hm[res_bsf_all_hm$log2FoldChange > 2,],
#          "smartseq_new/bsf_mapped/plots_de_common_genes/5_deseq_norm/7/pairwise_bsf_vs_all_down_bsf_up_rest.csv")
res_bsf_all_hm <- res_bsf_all_hm[res_bsf_all_hm$log2FoldChange < -2 | res_bsf_all_hm$log2FoldChange > 2,]





#heatmap
hm_gene_bsf <- unique(c(rownames(res_bsf_0h_hm), rownames(res_bsf_7d_hm), rownames(res_bsf_all_hm)))

#####p 0.01 expressed genes####

gene_exp <- (cm_common_genes[hm_gene_bsf,colData_bsf[colData_bsf$time == "0h", ]$col])
gene_exp <- gene_exp[rowSums(gene_exp) > 0,]




#### creating the heatmap
h0_bsf <- as.character(colData_bsf_pair[colData_bsf_pair$time == "0h", ]$col)
d7_bsf <- as.character(colData_bsf_pair[colData_bsf_pair$time == "7d", ]$col)
hbsf_bsf <- as.character(colData_bsf_pair[colData_bsf_pair$time == "bsf",]$col)
ht_order <- c(h0_bsf, d7_bsf, hbsf_bsf)

norm <- norm_bsf_pair[hm_gene_bsf, ht_order]
#norm <- df[hm_gene, rownames(colData_DE)]
scaled <- t(apply(norm, 1, scale))
colnames(scaled) <- colnames(norm)

dim(scaled)

max(scaled)
min(scaled)

p <- Heatmap(scaled, cluster_columns = F, cluster_rows = T,
             col = colorRamp2(seq(-2,2,0.1), viridis(41)),
             row_names_gp = gpar(fontsize = 2),
             show_column_names = TRUE, column_title = "BSF vs MCF + 7d, common genes, diff norm:
             absolute lfc > 2; padj<0.01")
print(p)




res_bsf_24h_vs_bsf_hm <- res_bsf_pair[res_bsf_pair$padj < 0.01,]
#write.csv(res_bsf_24h_vs_bsf_hm[res_bsf_24h_vs_bsf_hm$log2FoldChange < -2,],
#          "plots_de_common_genes/1_new_norm/expressed_genes_lfc_2/pairwise/bsf_vs_24h_overexpressed_24h.csv")
#write.csv(res_bsf_24h_vs_bsf_hm[res_bsf_24h_vs_bsf_hm$log2FoldChange > 2,],
#          "plots_de_common_genes/1_new_norm/expressed_genes_lfc_2/pairwise/bsf_vs_24h_underexpressed_24h.csv")

bsf_vs_0h_overxp <- (rownames(res_bsf_24h_vs_bsf_hm[res_bsf_24h_vs_bsf_hm$log2FoldChange < -2,]))
bsf_vs_0h_underxp <- rownames(res_bsf_24h_vs_bsf_hm[res_bsf_24h_vs_bsf_hm$log2FoldChange > 2,])
length(intersect(bsf_vs_0h_overxp, bsf_vs_0h_underxp))



#w p value and log fold cutoffs
EnhancedVolcano(res_bsf_pair,
                lab = rownames(res_bsf_pair),
                x = 'log2FoldChange',#
                y = 'padj',
                title = 'padj<0.01, lfc>2, bsf_vs_24h',
                pCutoff = 0.01,
                FCcutoff = 2)
