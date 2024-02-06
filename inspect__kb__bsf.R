#preprocessing of the count matrix
#obtained post-mapping with bustools

library(Matrix)
library(irlba)
library(ggplot2)
library(dplyr)
library(scico)
library(tidyverse)
library(Seurat)
theme_set(theme_bw())


setwd("bsf_mapped/")
#read in count matrices
cm <- readMM("counts_unfiltered/cells_x_genes.mtx")
cm <- t(cm)
dim(cm)
#[1] 11706 64

object.size(cm)#sparse matrix - zeroes are dots
#723704 bytes

counts_per_gene <- Matrix::rowSums(cm) #abundance of a gene across cells
cat("counts per gene: ", counts_per_gene[1:5], "\n")  ## counts for first 5 genes

genes_per_cell <- Matrix::colSums(cm > 0) # count gene only if it has non-zero reads mapped.
cat("counts for non-zero genes: ", genes_per_cell[1:5])  ## counts for first 5 genes

# only count cells where the gene is expressed
cells_per_gene <- rowSums(cm>0)
cat("count of cells with expressed genes: ", cells_per_gene[1:5])

hist(log10(genes_per_cell+1), main='genes per cell', col='wheat')
ggplot(as.data.frame(genes_per_cell), aes(x=log10(genes_per_cell+1)))+
  geom_histogram()+
  xlim(c(0,NA))


median(genes_per_cell)
mean(genes_per_cell)

hist(log10(counts_per_gene+1), main='counts per gene', col='wheat')
plot(sort(genes_per_cell), xlab='cell', log='y', main='genes per cell (ordered)')
plot((genes_per_cell), xlab='cell', log='y', main='genes per cell (not ordered)')



#####
# Convert to dgCMatrix, which is a compressed, sparse matrix format

cm_matrix <- Matrix::as.matrix(cm)
cm <- as(cm, "dgCMatrix")
dim(cm)
#11704(genes) 64 (cells)


genes <- read.delim("counts_unfiltered/cells_x_genes.genes.txt", sep="\n",
                    header = F, stringsAsFactors = F)
barcodes <- read.delim("counts_unfiltered/cells_x_genes.barcodes.txt", sep="\n",
                       header = F, stringsAsFactors = F)

colnames(cm) <- barcodes$V1
rownames(cm) <- genes$V1
cm <- as(cm, "dgCMatrix")


vsn::meanSdPlot(as.matrix(cm[rowVars(cm),]))
(matrixStats::rowVars(as.matrix(cm[rownames(var_cm),])))

cm_mean <- cm[rowMeans(as.matrix(cm))<10000,]
cm_var <- matrixStats::rowVars(as.matrix(cm_mean))
cm_mv <- data.frame(mean=rowMeans(as.matrix(cm_mean)), var=cm_var)


ggplot(data.frame(cm_mv), aes(x=mean, y=var)) +
  geom_point()+
  geom_abline(slope = 1, color="red")
  #geom_text(aes(label=mean))


##seurat workflow
tryp <- CreateSeuratObject(counts = cm, project = "tryp_new_kb", min.cells = 3, min.features = 200)
tryp


# Visualize QC metrics as a violin plot
options(repr.plot.width=12, repr.plot.height=6)
VlnPlot(tryp, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(tryp, feature1 = "nCount_RNA", feature2 = "rrna.percent")
plot2 <- FeatureScatter(tryp, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot2
plot1 + plot2

#subsetting
tryp <- subset(tryp, subset = nFeature_RNA>1500)
tryp

#normalising
tryp <- NormalizeData(tryp)#, normalization.method = "LogNormalize", scale.factor = 10000)

#normalising differently to align with normalisation method of orig dataset
cm_bsf_unnormalised <- as.matrix(GetAssayData(object = tryp, slot = "counts"))
cm_bsf <- as.matrix(GetAssayData(object = tryp, slot = "counts"))
dim(cm_bsf)

sf_bsf <- DESeq2::estimateSizeFactorsForMatrix(cm_bsf)
norm_bsf <- t(t(cm_bsf) / sf_bsf)


#recreate Seurat object with normalised counts
tryp <- CreateSeuratObject(counts = norm_bsf, project = "tryp_new_kb")#, min.cells = 3, min.features = 200)
tryp

#highly variable genes 
tryp <- FindVariableFeatures(tryp, selection.method = "vst", nfeatures = 1000)


# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(tryp), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(tryp)
plot1 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1

tryp@assays
#Top 10 variable features:
# Tb927.10.8590, Tb11.v5.0178, Tb927.4.880, Tb927.1.2980, Tb927.8.1210, Tb927.11.9190, Tb927.6.2960, Tb927.2.3020,
# Tb927.5.2600, Tb927.7.3090

# Top 10 variable features: using seurat's lognormalize
#  Tb927.10.8590, Tb11.v5.0178, Tb927.4.880, Tb927.1.2980, Tb927.8.1210, Tb927.11.9190, Tb927.6.2960,
# Tb927.2.3020, Tb927.5.2600, Tb927.7.3090

cm_bsf <- as.matrix(GetAssayData(object = tryp, slot = "data"))#for seurat normalised count matrix
cm_bsf <- norm_bsf #for normalised count matrix using deseq

