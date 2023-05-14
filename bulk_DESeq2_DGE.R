setwd("/home/sozaez/Desktop/TFM/bulkRNAseq/Analyses/Uroth_reg_bulkRNAseq/")

# very handy!
vignette("DESeq2")

# Libraries
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(ComplexHeatmap)
library(tibble)
library(tidyverse)
library(DESeq2)
library(pheatmap)
library(ggrepel)
library(gplots)
library(ashr)
library(EnhancedVolcano)
library(eulerr)
library(plotly)
library(wesanderson)
library(biomaRt)
library(DEGreport)
library(maSigPro)
library(MASS)
library(gridExtra)
library(genefilter)

############## Data loading ##############
data <- read.table("data/data_uro.txt", header = TRUE, row.names = 1)
meta <- read.table("meta/uroth_reg_meta.txt", header = TRUE, row.names = 1)

############## Merge count matrices and rename columns ##############
# Count matrix generation from HTSeq files
setwd("/home/sozaez/Desktop/TFM/bulkRNAseq/Uroth_RNAseq/HTSeq")

# Get list of all files in directory
files = list.files(path = ".", pattern = "trimmed")

# Merge the first two files and store
file1 = read.table(files[1], col.names = c("gene_id", files[1]))
file2 = read.table(files[2], col.names = c("gene_id", files[2]))
out.file = merge (file1, file2, by = c("gene_id"))

# For loop to merge contents of remaining files
for(i in 3:length(files))
{
  file = read.table(files[i], col.names = c("gene_id",files[i]))
  out.file <- merge(out.file, file, by = c("gene_id"))
}

# write.table(out.file, file = "htseq_total_count.tsv", sep="\t", row.names = FALSE)
# Rename columns
final <- out.file %>% 
  rename(
    day_0_S40 = X40_trimmedmerged_read_counts.txt,
    day_0_S41 = X41_trimmedmerged_read_counts.txt,
    day_0_S42 = X42_trimmedmerged_read_counts.txt,
    day_0_S43 = X43_trimmedmerged_read_counts.txt,
    day_0.5_S44 = X44_trimmedmerged_read_counts.txt,
    day_0.5_S45 = X45_trimmedmerged_read_counts.txt,
    day_0.5_S46 = X46_trimmedmerged_read_counts.txt,
    day_0.5_S47 = X47_trimmedmerged_read_counts.txt,
    day_1_S48 = X48_trimmedmerged_read_counts.txt,
    day_1_S49 = X49_trimmedmerged_read_counts.txt,
    day_1_S50 = X50_trimmedmerged_read_counts.txt,
    day_3_S52 = X52_trimmedmerged_read_counts.txt,
    day_3_S53 = X53_trimmedmerged_read_counts.txt,
    day_3_S54 = X54_trimmedmerged_read_counts.txt,
    day_3_S55 = X55_trimmedmerged_read_counts.txt,
    day_7_S56 = X56_trimmedmerged_read_counts.txt,
    day_7_S57 = X57_trimmedmerged_read_counts.txt,
    day_7_S58 = X58_trimmedmerged_read_counts.txt,
    day_7_S59 = X59_trimmedmerged_read_counts.txt,
    day_30_S60 = X60_trimmedmerged_read_counts.txt,
    day_30_S61 = X61_trimmedmerged_read_counts.txt,
    day_30_S62 = X62_trimmedmerged_read_counts.txt,
    day_30_S63 = X63_trimmedmerged_read_counts.txt,
  )
# write.table(final, file = "htseq_total_count.tsv", sep="\t", row.names = FALSE)


############## Conversion of ensembl gene IDs to gene symbols ##########
gene_ids <- data[1]
onlygene_ids <- tail(data[1], -5) # como gene_ids, pero sin __alignment_not_unique, __ambiguous...

# biomaRt
listEnsembl()
listEnsembl(version=108)

ensembl <- useEnsembl(biomart = "genes")
datasets <- listDatasets(ensembl)
ensembl.con <- useMart("ensembl", dataset = 'mmusculus_gene_ensembl')
attr <- listAttributes(ensembl.con)
filters <- listFilters(ensembl.con)

biomart_ids <- getBM(attributes = c('ensembl_gene_id','external_gene_name', 'mgi_symbol'),
                     filter = "ensembl_gene_id",
                     values = onlygene_ids$gene_id,
                     mart = ensembl.con)

write.table(biomart_ids, file = "biomart_ids.txt", sep="\t", row.names = TRUE)

# Evaluate match between input list of gene ID and the retrieved gene-gene name known pairs:
gene_id_match = match(onlygene_ids$gene_id, biomart_ids$ensembl_gene_id) # match for every element in gene_id entry in ensembl_gene_id (row no.)
# All elements match!

isSorted(gene_id_match)

onlygene_ids$symbol <- biomart_ids$mgi_symbol[gene_id_match] # addsuro gene symbol column

gaps <- which(onlygene_ids$symbol == "") # identify gaps in symbol column
ens_gap_genes <- onlygene_ids$gene_id[gaps] # get the ensembl id corresponding to those gaps

onlygene_ids$symbol[onlygene_ids$symbol == ""] <- ens_gap_genes # substitute gaps with ensembl ids

onlygene_ids[17541,] # check if the gene_id matches the symbol where there were gaps
onlygene_ids <- onlygene_ids$symbol
check <- which(complete_ids$gene_id == complete_ids$symbol) # check if ensembl ids match ensembl ids filled in gaps
matches <- match(biomart_ids$external_gene_name, biomart_ids$mgi_symbol)

data_symbols = merge(data, complete_ids, by = c("gene_id"))

write.table(data_symbols, file = "./data/data_symbols.txt", sep="\t", row.names = TRUE)

data_final <- data_symbols[, c("symbol", "day_0_S40", "day_0_S41", "day_0_S42", "day_0_S43", "day_0.5_S44", "day_0.5_S45", "day_0.5_S46", "day_0.5_S47", "day_1_S48", "day_1_S49", "day_1_S50", "day_3_S52", "day_3_S53", "day_3_S54", "day_3_S55", "day_7_S56", "day_7_S57", "day_7_S58", "day_7_S59", "day_30_S60", "day_30_S61", "day_30_S62", "day_30_S63")] # reorder columns

write.table(data_final, file = "./data/data_final.txt", sep="\t", row.names = TRUE)

# When checking the uniqueness of the mgi symbols, we observed that there are 88 nomenclatures of repeated genes (two different ensembl IDs corresponded to the same mgi symbol). Therefore, we have sorted the mgi symbols by their ensembl ID:

length(Data$symbol)-length(unique(Data$symbol))# then 88, now 0
dupl_id <- which(duplicated(Data$symbol))# duplication positions
Data$symbol[dupl_id] <- onlygene_ids$gene_id[dupl_id]# replace ensembl_id in duplicate positions

data2 <- Data[,-1]
write.table(data2, file = "data/data_uro.txt", sep="\t", row.names = TRUE)


############## 1. COUNT NORMALIZATION ##########
# 1. Check sample names match in both files: column names count matrix == row names metadata
all(colnames(data) %in% rownames(meta)) # TRUE
all(colnames(data) == rownames(meta)) # TRUE

## Creation DESeq2Dataset object
ddsuro <- DESeqDataSetFromMatrix(countData = data, 
                                 colData = meta, 
                                 design = ~timepoint)

levels(ddsuro$timepoint) # check the current levels
# reorder the timepoint levels (for PCA legend order)
ddsuro$timepoint <- factor(ddsuro$timepoint, levels = levels(ddsuro$timepoint)[c(1:3, 5:6, 4)])
ddsuro <- estimateSizeFactors(ddsuro)
levels(ddsuro$timepoint) # check levels after reasignation

# Retrieve the normalized counts matrix from ddsuro
normalized_counts <- counts(ddsuro, normalized = TRUE)
vsd <- vst(ddsuro, blind = FALSE)
vst_counts <- assay(vsd)
vst_counts_df <- data.frame(vst_counts)

# vst_counts_df <- tibble::rownames_to_column(vst_counts_df, "genes") # convert row.names into first column 
# write.table(vst_counts_df, file = "data/vst_counts_df.xls", sep="\t", row.names = TRUE)

############## 2. QUALITY CONTROL ANALYSIS ##########
# Transform counts just for data visualization:
rld <- rlog(ddsuro, blind = TRUE)
rld_mat <- assay(rld)# function from the "SummarizedExperiment" package
rld_cor <- cor(rld_mat) # compute the pairwise correlation values for samples
head(rld_cor)   ## check the output of cor() rownames and colnames

# write.table(rld_mat, file = "data/rld_mat.xls", sep="\t", row.names = TRUE)

###     2.1. PCA ##########
plotPCA(rld, intgroup = "timepoint")
pcaData <- plotPCA(rld, intgroup = "timepoint", returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

# label samples
ggplot(pcaData, aes(PC1, PC2, color=timepoint)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  ggtitle("Urothelium samples") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_brewer(palette = "Set2") +
  geom_text(aes(label = rownames(pcaData)), size = 3.5) +
  coord_fixed()

# Retrieve PC3 and PC4 information
pca <- prcomp(t(rld_mat)) # prcomp() principal components info
mat_pca <- t(rld_mat)
df <- cbind(meta, pca$x) # data frame with metadata and PCs values
ggplot(df) + geom_point(aes(x = PC3, y = PC4, color = timepoint))

# Total variance explained by each principal component
var_explained = pca$sdev^2 / sum(pca$sdev^2)

# Scree plot
qplot(c(1:23), var_explained) + 
  geom_line() + 
  xlab("Principal Component") + 
  ylab("Variance Explained") +
  ggtitle("Scree Plot") +
  ylim(0, 0.3)

# 3D PCA:
ntop <- 500
rv <- rowVars(assay(rld))
# pick the top 500 genes by row variance:
select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
# subset data to those top 500 genes, transpose and get PCs info
meta$timepoint <- factor(meta$timepoint, levels=unique(meta$timepoint))
pcs_top500 <- prcomp(t(assay(rld)[select, ]))
pca_top500_df <- cbind(meta, pcs_top500$x)
colors <- brewer.pal(n = 6, name = "Set2")

meta$timepoint <- factor(meta$timepoint, levels=unique(meta$timepoint))
legend_group <- factor(meta$timepoint, levels=unique(meta$timepoint))
legend_group <- unique(meta$timepoint)

fig3D <- plot_ly(pca_top500_df, x = ~PC1, y = ~PC2, z = ~PC3, type = "scatter3d", mode = "markers", color = ~timepoint, size = 3, colors = colors)
fig3D <-  fig3D %>% add_markers() %>% layout(title = "3D PCA Urothelium samples", legend = list(traceorder = "reversed")) %>% add_trace(showlegend = FALSE)
fig3D


###     2.2. HIERARCHICAL CLUSTERING ##########
## Plot heatmap
# Sample distance heatmap
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( vsd$timepoint )
colnames(sampleDistMatrix) <- paste( vsd$timepoint )

colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,cluster_rows = FALSE, cluster_cols = FALSE,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors,
         main = "QC sample heatmap")


############## 3. DESeq2 ANALYSIS ##########
# DE analysis - timepoint as factor (desing ~ timepoint)
ddsuro <- DESeq(ddsuro)
res <- results(ddsuro)
# write.table(res, file="results/wald_res.xls", sep="\t", quote=F, col.names=NA)

resultsNames(ddsuro) # Contrasts

## Plot dispersion estimates
plotDispEsts(ddsuro)

# Top Gene Plot
topGene <- rownames(res)[which.min(res$padj)] # finds gene w/ lowest padj
plots_counts <- plotCounts(ddsuro, gene = "Ttr", intgroup=c("timepoint"), returnData = TRUE)

# Turn timepoint into a factor with the levels in the correct order
plots_counts$timepoint <- factor(plots_counts$timepoint, levels=unique(plots_counts$timepoint))

ggplot(plots_counts, aes(timepoint, count, color=timepoint)) +
  geom_point(size=3) +
  xlab("timepoint") +
  ylab("counts") + 
  ggtitle(paste0("Gene: ", topGene[1], " over time")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_brewer(palette = "Set2")

# Label samples
geneCounts <- plotCounts(ddsuro, gene = topGene, intgroup = c("timepoint"), returnData = TRUE)
ggplot(geneCounts, aes(x = timepoint, y = count, color = timepoint, group = timepoint)) + scale_y_log10() +  geom_beeswarm(cex = 3) + geom_text(
  label=rownames(meta),
  nudge_x = 0.25, nudge_y = 0.25,
  check_overlap = T)


###     3.1. Hypothesis testing: WALD TEST ##########
# We have 6 sample cases, we can make 15 possible pairwise comparisons.
# Define contrasts. The name provided in the second element is the level that is used as baseline

contrast_12_0 <- c("timepoint", "12h", "0h")
contrast_24_0 <- c("timepoint", "24h", "0h")
contrast_72_0 <- c("timepoint", "72h", "0h")
contrast_7d_0 <- c("timepoint", "7d", "0h")
contrast_30d_0 <- c("timepoint", "30d", "0h")
contrast_24_12 <- c("timepoint", "24h", "12h")
contrast_72_24 <- c("timepoint", "72h", "24h")
contrast_7d_72 <- c("timepoint", "7d", "72h")
contrast_30d_7d <- c("timepoint", "30d", "7d")

# Extract results table and shrink the log2 fold changes
# - 0h - 12h
res_table12h <- results(ddsuro, contrast = contrast_12_0, alpha = 0.05, lfcThreshold = 0.58)
res_table12_Shr <- lfcShrink(ddsuro, contrast = contrast_12_0, res = res_table12h, type = "ashr")
# - 0h - 24h
res_table24h <- results(ddsuro, contrast = contrast_24_0, alpha = 0.05, lfcThreshold = 0.58)
res_table24_Shr <- lfcShrink(ddsuro, contrast = contrast_24_0, res = res_table24h, type = "ashr")
# - 0h - 72h
res_table72h <- results(ddsuro, contrast = contrast_72_0, alpha = 0.05, lfcThreshold = 0.58)
res_table72_Shr <- lfcShrink(ddsuro, contrast = contrast_72_0, res = res_table72h, type = "ashr")
# - 0h - 7d
res_table7d <- results(ddsuro, contrast = contrast_7d_0, alpha = 0.05, lfcThreshold = 0.58)
res_table7d_Shr <- lfcShrink(ddsuro, contrast = contrast_7d_0, res = res_table7d, type = "ashr")
# - 0h - 30d
res_table30d <- results(ddsuro, contrast = contrast_30d_0, alpha = 0.05, lfcThreshold = 0.58)
res_table30d_Shr <- lfcShrink(ddsuro, contrast = contrast_30d_0, res = res_table30d, type = "ashr")

# - 24h - 12h
res_table24h_12h <- results(ddsuro, contrast = contrast_24_12, alpha = 0.05, lfcThreshold = 0.58)
res_table24h_12h_Shr <- lfcShrink(ddsuro, contrast = contrast_24_12, res = res_table24h_12h, type = "ashr")
# - 72h - 24h
res_table72h_24h <- results(ddsuro, contrast = contrast_72_24, alpha = 0.05, lfcThreshold = 0.58)
res_table72h_24h_Shr <- lfcShrink(ddsuro, contrast = contrast_72_24, res = res_table72h_24h, type = "ashr")
# - 7d - 72h
res_table7d_72h <- results(ddsuro, contrast = contrast_7d_72, alpha = 0.05, lfcThreshold = 0.58)
res_table7d_72h_Shr <- lfcShrink(ddsuro, contrast = contrast_7d_72, res = res_table7d_72h, type = "ashr")
# - 30d - 7d
res_table30d_7d <- results(ddsuro, contrast = contrast_30d_7d, alpha = 0.05, lfcThreshold = 0.58)
res_table30d_7d_Shr <- lfcShrink(ddsuro, contrast = contrast_30d_7d, res = res_table30d_7d, type = "ashr")

# If a row contains a sample with an extreme count outlier then the p-value and adjusted p-value will be set to NA.
# If a row is filtered by automatic independent filtering, for having a low mean normalized count, then only the adjusted p-value will be set to NA.

# Convert the results table into a tibble.
res_12h_tb <- res_table12_Shr %>% data.frame() %>% rownames_to_column(var="gene") %>% as_tibble()
res_24h_tb <- res_table24_Shr %>% data.frame() %>% rownames_to_column(var="gene") %>% as_tibble()
res_72h_tb <- res_table72_Shr %>% data.frame() %>% rownames_to_column(var="gene") %>% as_tibble()
res_7d_tb <- res_table7d_Shr %>% data.frame() %>% rownames_to_column(var="gene") %>%  as_tibble()
res_30d_tb <- res_table30d_Shr %>% data.frame() %>% rownames_to_column(var="gene") %>% as_tibble()

res_24h_12h_tb <- res_table24h_12h_Shr %>% data.frame() %>% rownames_to_column(var="gene") %>% as_tibble()
res_72h_24h_tb <- res_table72h_24h_Shr %>% data.frame() %>% rownames_to_column(var="gene") %>% as_tibble()
res_7d_72h_tb <- res_table7d_72h_Shr %>% data.frame() %>% rownames_to_column(var="gene") %>% as_tibble()
res_30d_7d_tb <- res_table30d_7d_Shr %>% data.frame() %>% rownames_to_column(var="gene") %>% as_tibble()

# Alternatively, set results as data frame table:
res_12h_df <- data.frame(res_table12_Shr) %>% arrange(padj)
# write.table(res_12h_df, file="results/res_12h_df.xls", sep="\t", quote=F, col.names=NA)
res_24h_df <- data.frame(res_table24_Shr) %>% arrange(padj)
# write.table(res_24h_df, file="results/res_24h_df.xls", sep="\t", quote=F, col.names=NA)
res_72h_df <- data.frame(res_table72_Shr) %>% arrange(padj)
# write.table(res_72h_df, file="results/res_72h_df.xls", sep="\t", quote=F, col.names=NA)
res_7d_df <- data.frame(res_table7d_Shr) %>% arrange(padj)
# write.table(res_7d_df, file="results/res_7d_df.xls", sep="\t", quote=F, col.names=NA)
res_30d_df <- data.frame(res_table30d_Shr) %>% arrange(padj)
# write.table(res_30d_df, file="results/res_30d_df.xls", sep="\t", quote=F, col.names=NA)

res_24h_12h_df <- data.frame(res_table24h_12h_Shr) %>% arrange(padj)
# write.table(res_24h_12h_df, file="results/res_24h_12h_df.xls", sep="\t", quote=F, col.names=NA)
res_72h_24h_df <- data.frame(res_table72h_24h_Shr) %>% arrange(padj)
# write.table(res_72h_24h_df, file="results/res_72h_24h_df.xls", sep="\t", quote=F, col.names=NA)
res_7d_72h_df <- data.frame(res_table7d_72h_Shr) %>% arrange(padj)
# write.table(res_7d_72h_df, file="results/res_7d_72h_df.xls", sep="\t", quote=F, col.names=NA)
res_30d_7d_df <- data.frame(res_table30d_7d_Shr) %>% arrange(padj)
# write.table(res_30d_7d_df, file="results/res_30d_7d_df.xls", sep="\t", quote=F, col.names=NA)


# Filter significant genes using the filter() function
padj.cutoff <- 0.05 # already done; FDR = 5% of positives after BH correction
lfc.cutoff <- 0.58 # log2 fold change > 1.5 cutoff
lfc.cutoff_1.3 <- 0.38 # log2 fold change > 1.3 cutoff

# Subsetting the tibble (or data frame) table:
sigDE12h <- res_12h_tb %>% 
  dplyr::filter(padj < padj.cutoff & abs(log2FoldChange) > lfc.cutoff) %>% 
  arrange(padj)
sigDE24h <- res_24h_tb %>% 
  dplyr::filter(padj < padj.cutoff & abs(log2FoldChange) > lfc.cutoff) %>% 
  arrange(padj)
sigDE72h <- res_72h_tb %>% 
  dplyr::filter(padj < padj.cutoff & abs(log2FoldChange) > lfc.cutoff) %>% 
  arrange(padj)
sigDE7d <- res_7d_tb %>% 
  dplyr::filter(padj < padj.cutoff & abs(log2FoldChange) > lfc.cutoff) %>% 
  arrange(padj)
sigDE30d <- res_30d_tb %>% 
  dplyr::filter(padj < padj.cutoff & abs(log2FoldChange) > lfc.cutoff) %>% 
  arrange(padj)

sigDE24h12h <- res_24h_12h_tb %>% 
  dplyr::filter(padj < padj.cutoff & abs(log2FoldChange) > lfc.cutoff) %>% 
  arrange(padj)
sigDE72h24h <- res_72h_24h_tb %>% 
  dplyr::filter(padj < padj.cutoff & abs(log2FoldChange) > lfc.cutoff) %>% 
  arrange(padj)
sigDE7d72h <- res_7d_72h_tb %>% 
  dplyr::filter(padj < padj.cutoff & abs(log2FoldChange) > lfc.cutoff) %>% 
  arrange(padj)
sigDE30d7d <- res_30d_7d_tb %>% 
  dplyr::filter(padj < padj.cutoff & abs(log2FoldChange) > lfc.cutoff) %>% 
  arrange(padj)

# write.table(sigDE12h, file="results/sigDE12_.xls", sep="\t", quote=F, col.names=NA)


# Select only the 20 first DEGs (for volcano plot):
gene_labels12h <- head(sigDE12h$gene, 20)
gene_labels24h <- head(sigDE24h$gene, 20)
gene_labels72h <- head(sigDE72h$gene, 20)
gene_labels7d <- head(sigDE7d$gene, 20)
gene_labels30d <- head(sigDE30d$gene, 20)
gene_labels24h12h <- head(sigDE24h12h$gene, 20)
gene_labels72h24h <- head(sigDE72h24h$gene, 20)
gene_labels7d72h <- head(sigDE7d72h$gene, 20)
gene_labels30d7d <- head(sigDE30d7d$gene, 20)


# Volcano plots:
EnhancedVolcano(res_table12_Shr, lab = rownames(res_table12_Shr), x = 'log2FoldChange', y = 'padj', title = '12h vs 0h', pCutoff = 0.05, FCcutoff = 1.5, selectLab = gene_labels12h, drawConnectors = TRUE, max.overlaps = 23, legendPosition = 'right', legendLabels=c('not sig.','log2FC','padj', 'padj & log2FC')) #pCutoff relates to whatever value we supply as the y-axis. Padj here.

EnhancedVolcano(res_table24_Shr, lab = rownames(res_table24_Shr), x = 'log2FoldChange', y = 'padj', title = '24h vs 0h', pCutoff = 0.05, FCcutoff = 1.5, selectLab = gene_labels24h, drawConnectors = TRUE, max.overlaps = 23, legendPosition = 'right', legendLabels=c('not sig.','log2FC','padj', 'padj & log2FC'))

EnhancedVolcano(res_table72_Shr, lab = rownames(res_table72_Shr), x = 'log2FoldChange', y = 'padj', title = '72h vs 0h', pCutoff = 0.05, FCcutoff = 1.5, selectLab = gene_labels72h, drawConnectors = TRUE, max.overlaps = 23, legendPosition = 'right', legendLabels=c('not sig.','log2FC','padj', 'padj & log2FC'))

EnhancedVolcano(res_table7d_Shr, lab = rownames(res_table7d_Shr), x = 'log2FoldChange', y = 'padj', title = '7d vs 0h', pCutoff = 0.05, FCcutoff = 1.5, selectLab = gene_labels7d, drawConnectors = TRUE, max.overlaps = 23, legendPosition = 'right', legendLabels=c('not sig.','log2FC','padj', 'padj & log2FC'))

EnhancedVolcano(res_table30d_Shr, lab = rownames(res_table30d_Shr), x = 'log2FoldChange', y = 'padj', title = '30d vs 0h', pCutoff = 0.05, FCcutoff = 1.5, drawConnectors = TRUE, legendPosition = 'right', legendLabels=c('not sig.','log2FC','padj', 'padj & log2FC'))

EnhancedVolcano(res_table24h_12h_Shr, lab = rownames(res_table24h_12h_Shr), x = 'log2FoldChange', y = 'padj', title = '24h vs 12h', pCutoff = 0.05, FCcutoff = 1.5, selectLab = gene_labels24h12h, drawConnectors = TRUE, max.overlaps = 23, legendPosition = 'right', legendLabels=c('not sig.','log2FC','padj', 'padj & log2FC'))

EnhancedVolcano(res_table72h_24h_Shr, lab = rownames(res_table72h_24h_Shr), x = 'log2FoldChange', y = 'padj', title = '72h vs 24h', pCutoff = 0.05, FCcutoff = 1.5, selectLab = gene_labels72h24h, drawConnectors = TRUE, max.overlaps = 23, legendPosition = 'right', legendLabels=c('not sig.','log2FC','padj', 'padj & log2FC'))
                
EnhancedVolcano(res_table7d_72h_Shr, lab = rownames(res_table7d_72h_Shr), x = 'log2FoldChange', y = 'padj', title = '7d vs 72h', pCutoff = 0.05, FCcutoff = 1.5, selectLab = gene_labels7d72h, drawConnectors = TRUE, max.overlaps = 23, legendPosition = 'right', legendLabels=c('not sig.','log2FC','padj', 'padj & log2FC'))

EnhancedVolcano(res_table30d_7d_Shr, lab = rownames(res_table30d_7d_Shr), x = 'log2FoldChange', y = 'padj', title = '30d vs 7d', pCutoff = 0.05, FCcutoff = 1.5, selectLab = gene_labels30d7d, drawConnectors = TRUE, max.overlaps = 23, legendPosition = 'right', legendLabels=c('not sig.','log2FC','padj', 'padj & log2FC'))

# (file="plots/Volcano/volcano_30d7d.", width = 12)
# EnhancedVolcano(res_table30d_7d_Shr, lab = rownames(res_table30d_7d_Shr), x = 'log2FoldChange', y = 'padj', title = '30d vs 7d', pCutoff = 0.05, FCcutoff = 1.5, selectLab = gene_labels30d7d, drawConnectors = TRUE, max.overlaps = 23, legendPosition = 'right', legendLabels=c('not sig.','log2FC','padj', 'padj & log2FC'))
# dev.off()


# Shrunken results:
plotMA(res_table12_Shr, ylim=c(-2,2))
plotMA(res_table24_Shr, ylim=c(-2,2))
plotMA(res_table72_Shr, ylim=c(-2,2))
plotMA(res_table7d_Shr, ylim=c(-2,2))
plotMA(res_table30d_Shr, ylim=c(-2,2))
plotMA(res_table24h_12h_Shr, ylim=c(-2,2))
plotMA(res_table72h_24h_Shr, ylim=c(-2,2))
plotMA(res_table7d_72h_Shr, ylim=c(-2,2))
plotMA(res_table30d_7d_Shr, ylim=c(-2,2))


# Venn diagrams:
# Retrieve only significant genes up and down regulated by L2FC
position_up12h <- which(sigDE12h$log2FoldChange >= 0)
DEGup_12h <- sigDE12h$gene[position_up12h]
position_down12h <- which(sigDE12h$log2FoldChange < 0)
DEGdown_12h <- sigDE12h$gene[position_down12h]

position_up24h <- which(sigDE24h$log2FoldChange >= 0)
DEGup_24h <- sigDE24h$gene[position_up24h]
position_down24h <- which(sigDE24h$log2FoldChange < 0)
DEGdown_24h <- sigDE24h$gene[position_down24h]

position_up72h <- which(sigDE72h$log2FoldChange >= 0)
DEGup_72h <- sigDE72h$gene[position_up72h]
position_down72h <- which(sigDE72h$log2FoldChange < 0)
DEGdown_72h <- sigDE72h$gene[position_down72h]

position_up7d <- which(sigDE7d$log2FoldChange >= 0)
DEGup_7d <- sigDE7d$gene[position_up7d]
position_down7d <- which(sigDE7d$log2FoldChange < 0)
DEGdown_7d <- sigDE7d$gene[position_down7d]

position_up2412 <- which(sigDE24h12h$log2FoldChange >= 0)
DEGup_2412 <- sigDE24h12h$gene[position_up2412]
position_down2412 <- which(sigDE24h12h$log2FoldChange < 0)
DEGdown_2412 <- sigDE24h12h$gene[position_down2412]

position_up7224 <- which(sigDE72h24h$log2FoldChange >= 0)
DEGup_7224 <- sigDE72h24h$gene[position_up7224]
position_down7224 <- which(sigDE72h24h$log2FoldChange < 0)
DEGdown_7224 <- sigDE72h24h$gene[position_down7224]

position_up7d72 <- which(sigDE7d72h$log2FoldChange >= 0)
DEGup_7d72 <- sigDE7d72h$gene[position_up7d72]
position_down7d72 <- which(sigDE7d72h$log2FoldChange < 0)
DEGdown_7d72 <- sigDE7d72h$gene[position_down7d72]

position_up30d7d <- which(sigDE30d7d$log2FoldChange >= 0)
DEGup_30d7d <- sigDE30d7d$gene[position_up30d7d]
position_down30d7d <- which(sigDE30d7d$log2FoldChange < 0)
DEGdown_30d7d <- sigDE30d7d$gene[position_down30d7d]

# Make list
list_all_DEup <- list("12h (472)" = DEGup_12h, "24h (316)" = DEGup_24h, "72h (42)" = DEGup_72h, "7d (19)" = DEGup_7d)
list_all_DEdown <- list("12h (744)" = DEGdown_12h, "24h (144)" = DEGdown_24h, "72h (8)" = DEGdown_72h, "7d (4)" = DEGdown_7d)

plot(venn(list_all_DEup), quantities = TRUE, legend = TRUE, labels = TRUE, adjust_labels = TRUE, main = "upreg DEGs Wald test")
plot(venn(list_all_DEdown), quantities = TRUE, legend = TRUE, labels = TRUE, adjust_labels = TRUE, main = "downreg DEGs Wald test")

plot(euler(list_all_DEup), quantities = TRUE, legend = TRUE,  adjust_labels = TRUE, main = "upregulated DEGs Wald test")
plot(euler(list_all_DEdown), quantities = TRUE, legend = TRUE, adjust_labels = TRUE, main = "downregulated DEGs Wald test")


# Heatmap:
# DEGs names extraction from Wald test pairwise comparisons:
DEGs12h <- sigDE12h[,1]
DEGs12h <- DEGs12h$gene

DEGs24h <- sigDE24h[,1]
DEGs24h <- DEGs24h$gene

DEGs72h <- sigDE72h[,1]
DEGs72h <- DEGs72h$gene

DEGs7d <- sigDE7d[,1]
DEGs7d <- DEGs7d$gene

DEGs30d <- sigDE30d[,1]
DEGs30d <- DEGs30d$gene

DEGs24h12h <- sigDE24h12h[,1]
DEGs24h12h <- DEGs24h12h$gene

DEGs72h24h <- sigDE72h24h[,1]
DEGs72h24h <- DEGs72h24h$gene

DEGs7d72h <- sigDE7d72h[,1]
DEGs7d72h <- DEGs7d72h$gene

DEGs30d7d <- sigDE30d7d[,1]
DEGs30d7d <- DEGs30d7d$gene

# Merge all Wald DEGs
DEGs_chr <- c(DEGs12h, DEGs24h, DEGs72h, DEGs7d, DEGs24h12h, DEGs72h24h, DEGs7d72h, DEGs30d7d)
all_DEGs <- unique(DEGs_chr) # select unique DEGs
all_DEGs_df <- as.data.frame(all_DEGs)
vst_limited <- vst_counts_df[match(all_DEGs_df$all_DEGs, row.names(vst_counts_df)),]

pheatmap(vst_limited, cluster_rows=TRUE, show_rownames=FALSE,
         cluster_cols=FALSE, scale = "row", main="Wald DEGs up & down, n = 1503")


###     3.2. Hypothesis testing: LIKELIHOOD RATIO TEST (LRT) ####
# Hypotheses:
# H0: no difference in models (full and reduced) 
# H1: the full model is statistically "better" than the reduced
ddsuro_lrt <- DESeq(ddsuro, test="LRT", reduced = ~1) # reduced = ~ timepoint
res_lrt <- results(ddsuro_lrt)
resultsNames(ddsuro_lrt)

res_lrt_tb <- res_lrt %>% data.frame() %>% rownames_to_column(var="gene") %>% as_tibble() %>% arrange(padj)
# write.table(res_lrt_tb, "results/res_lrt_tb.xls", row.names = FALSE)

# NOTE: If we cast 'day' as a numeric, then it's a linear regression, and the fold-change represents the expected change in log-expression for each day that passes. Thus, when filtering significant genes from the LRT we use only the FDR as our threshold.
# See https://support.bioconductor.org/p/9139197/

## Extracting significant differentially expressed genes. 
padj.cutoff_0.001 <- 0.001
padj.cutoff_0.01 <- 0.01
padj.cutoff <- 0.05

# Get significant gene lists
sig_lrt <- res_lrt %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>% 
  filter(padj < padj.cutoff)

sig_lrt_0.01 <- res_lrt %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>% 
  filter(padj < padj.cutoff_0.01)

sig_lrt_0.001 <- res_lrt %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>% 
  filter(padj < padj.cutoff_0.001)
# write.table(sig_LRT_day, "results/sig_LRT_day.xls", row.names = FALSE)

# Get sig gene lists
sig_lrt <- sig_lrt %>% 
  pull(gene)

sig_lrt_0.01 <- sig_lrt_0.01 %>% 
  pull(gene)

sig_lrt_0.001 <- sig_lrt_0.001 %>% 
  pull(gene)

length(sig_lrt)
length(sig_lrt_0.01)
length(sig_lrt_0.001)

# Compare to numbers we had from Wald test
nrow(sigDE12h)
nrow(sigDE24h)
nrow(sigDE72h)
nrow(sigDE7d)
nrow(sigDE30d)
nrow(sigDE24h12h)
nrow(sigDE72h24h)
nrow(sigDE7d72h)
nrow(sigDE30d7d)

vst_DEGs_lrt01 <- vst_counts_df[match(sig_lrt_0.01, row.names(vst_counts_df)),]
# write.table(vst_DEGs_lrt01, file = "data/vst_DEGs_lrt01.xls", sep="\t", row.names = TRUE)

# VENN diagram Wald - LRT
listDEGs_wald_lrt <- list("Wald" = all_DEGs, "LRT" = sig_lrt_0.01)
listDEGs_wald_lrtfull <- list("Wald" = all_DEGs, "LRT" = sig_lrt)

plot(venn(listDEGs_wald_lrt), quantities = TRUE, legend = TRUE, labels = TRUE, adjust_labels = TRUE, main = "DEGs Wald & LRT")

plot(euler(listDEGs_wald_lrt), quantities = TRUE, legend = TRUE, labels = TRUE, adjust_labels = TRUE, main = "DEGs Wald & LRT")

commonDEG_wald_lrt <- intersect(all_DEGs, sig_lrt_0.01)
onlyWald_DEG <- setdiff(all_DEGs, sig_lrt_0.01)

# Lets look at the expression of some genes of interest:
genes_interest <- c("Krt5", "Krt6a", "Krt14", "Krt18", "Krt20", "Upk1a", "Upk1b", "Upk3a", "Trp63", "Slpi", "Fosl1", "Cdkn1a", "Atg9b", "Dtl", "Daglb", "Ppp1r18os")

vst_genes_interest <- vst_counts_df[match(genes_interest, row.names(vst_counts_df)),]

pheatmap(vst_genes_interest, cluster_rows=TRUE, cluster_cols=FALSE, main = "Gene VST counts vs. sampletype", scale = "row")

# ("plots/relevantgenes_lrt_heatmap.", width=12, height=8)
# pheatmap(vst_DEGs_lrt, cluster_rows=FALSE, cluster_cols=FALSE)
# dev.off()

wald_lrt <- res_lrt_tb[match(all_DEGs_df$all_DEGs, res_lrt_tb$gene),]
wald_lrt <- arrange(wald_lrt, padj)

norm_wald <- normalized_counts[match(wald_lrt$gene, row.names(normalized_counts)),]
pheatmap(norm_wald[1:50,], cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=FALSE, scale = "row", main="Wald DEGs up & down normalized counts")

rld_wald <- rld_mat[match(wald_lrt$gene, row.names(rld_mat)),]
pheatmap(rld_wald[1:50,], cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=FALSE, scale = "row", main="Wald DEGs up & down rlog")


############## 4. Wald vs. LRT ##################
# Top 20 variable genes - heatmap
topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 20)
mat_Var  <- assay(vsd)[ topVarGenes, ]
mat_Var <- mat_Var - rowMeans(mat_Var)
anno <- as.data.frame(colData(vsd)[, c("timepoint","day")])

pheatmap(mat_Var, annotation_col = anno, cluster_cols = FALSE, main = "Top 20 variable genes") # clustered by row

pheatmap(mat_Var, annotation_col = anno, cluster_cols = F, cluster_rows = F, main = "Top 20 variable genes") # not clustered

# ("plots/top20var_clust_heatmap.", width=12, height=8)
# pheatmap(mat_Var, annotation_col = anno, cluster_cols = T, cluster_rows = F, main = "Top 20 variable genes")
# dev.off()

# Top 20 most variable genes LRT - heatmap
top20_lrt <- head(res_lrt_tb, 20)
mat_Var_lrt <- vst_counts_df[top20_lrt$gene,]
mat_Var_lrt <- mat_Var_lrt - rowMeans(mat_Var_lrt)

pheatmap(mat_Var_lrt, annotation_col = anno, cluster_cols = F, cluster_rows = T, main = "Top 20 variable genes LRT") # clustered by row
pheatmap(mat_Var_lrt, annotation_col = anno, cluster_cols = F, cluster_rows = F, main = "Top 20 variable genes LRT") # not clustered

# ("plots/top20var_LRT_heatmap.", width=12, height=8)
# pheatmap(mat_Var_lrt, annotation_col = anno, cluster_cols = F, cluster_rows = F, main = "Top 20 variable genes LRT")
# dev.off()

# Retrieve necessary data
res_tb <- res %>% data.frame() %>% rownames_to_column(var="gene") %>% as_tibble() %>% arrange(padj)  
padj.cutoff <- 0.05

merged_table <- rbind(res_12h_tb, res_24h_tb, res_72h_tb, res_7d_tb, res_30d_tb, res_24h_12h_tb, res_72h_24h_tb, res_7d_72h_tb, res_30d_7d_tb)

wald_df <- merged_table %>% group_by(gene) %>% filter(padj == min(padj) & padj != 1)
df_filtered <- merged_table %>% group_by(gene) %>% filter(padj == min(padj)) %>% dplyr::filter(padj < padj.cutoff) %>% arrange(padj)

lrt_df <- res_lrt_tb[match(wald_df$gene, res_lrt_tb$gene),]
lrt_filtered <- res_lrt_tb[match(row.names(data), res_lrt_tb$gene),]
lrt_filtered <- res_lrt_tb %>% filter(padj < padj.cutoff)


# Create volcano plot for Wald test
p1 <- ggplot(df_filtered, aes(x = log2FoldChange, y = -log10(padj))) + 
  geom_point(aes(color = log2FoldChange > 0 & padj < 0.05), size = 1.5, alpha = 0.8) +
  scale_color_manual(values = c("blue", "red")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  labs(title = "Volcano plot (Wald test)",
       x = "log2(fold change)",
       y = "-log10(adjusted p-value)",
       color = "Significant",
       caption = "Blue points: Downregulated genes; Red points: Upregulated genes") +
  theme_bw()

# Create volcano plot for LRT test
p2 <- ggplot(lrt_filtered, aes(x = log2FoldChange, y = -log10(padj))) + 
  geom_point(aes(color = log2FoldChange > 0 & padj < 0.05), size = 1.5, alpha = 0.8) +
  scale_color_manual(values = c("blue", "red")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  labs(title = "Volcano plot (LRT test)",
       x = "log2(fold change)",
       y = "-log10(adjusted p-value)",
       color = "Significant",
       caption = "Blue points: Downregulated genes; Red points: Upregulated genes") +
  theme_bw()

# Arrange plots side by side
grid.arrange(p1, p2, ncol = 2)

# Extract the -log10(p.adj) values for each data frame
wald_vals <- -log10(wald_df$pvalue)
lrt_vals <- -log10(lrt_df$pvalue)

# Create a data frame with the values from both tests
plot_df <- data.frame(wald_vals, lrt_vals)

# Create the scatter plot
ggplot(plot_df, aes(x = wald_vals, y = lrt_vals)) +
  geom_point() +
  xlab("-log10(padj) Wald test") +
  ylab("-log10(padj) LRT test")

# Merge the two data frames by gene names
merged_df <- merge(wald_df, lrt_df, by = "gene", suffixes = c("_wald", "_lrt"))

# Subset to keep only rows where the p.adj is not 1 and is the minimum for the gene
min_padj_df <- merged_df %>% group_by(gene) %>% filter(padj_wald == min(padj_wald), padj_wald != 1)

# Create the plot
# pdf("plots/lrt_vs_wald_plot.pdf", width=12, height=8)
ggplot(min_padj_df, aes(-log10(padj_wald), -log10(padj_lrt))) +
  geom_point() +
  xlab("-log10 Wald adj. p-value") +
  ylab("-log10 LRT adj. p-value") +
  ggtitle("Wald vs LRT test results") +
  theme_bw() +
  geom_text(aes(-log10(padj_wald), -log10(padj_lrt), label=gene), 
            hjust=0, vjust=0, nudge_x=0.3, nudge_y=0.3, size=3, data=subset(min_padj_df, padj_lrt < 0.05) %>% top_n(20, padj_lrt))
# dev.off()

############## 5. GENE CLUSTERS IDENTIFICATION ####
###     5.1. degPatterns ####
# degPatterns groups the genes based on their changes in expression across sample groups
# We could perform clustering to identify genes that change over time in a way meaningful to us. Subset results for faster cluster finding. Only DEGs.
sig_lrt_0.01 <- res_lrt %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>% 
  filter(padj < padj.cutoff_0.01)

clust_sig_genes <- sig_lrt_0.01 %>%
  arrange(padj)

#  %>% head(n=1000)

rld_LRT <- rlog(ddsuro_lrt, blind=T)
rld_mat_LRT <- assay(rld_LRT)

# Obtain rlog values for those significant genes
cluster_rlog <- rld_mat_LRT[clust_sig_genes$gene, ]

# Use the `degPatterns` function from the 'DEGreport' package to show gene clusters across sample groups
meta$index <- as.factor(meta$index)
clusters <- degPatterns(cluster_rlog, metadata = meta, time="hours", col="index", scale = TRUE)

df_clusters <- clusters$summarise
df_clusters <- arrange(df_clusters, -n_genes)
df_clusters <- df_clusters[, c(2,6)]
df_clusters <- unique(df_clusters)
df_clusters$cluster <- factor(df_clusters$cluster, levels = df_clusters$cluster)

# Genes per cluster
df_clusters %>%
  ggplot(aes(x=cluster, y=n_genes))+
  geom_col(color="black", fill="#69b3a2")+
  labs(title="degPatterns genes per cluster") +
  xlab("Cluster") + 
  ylab("Number of genes per cluster") +
  geom_hline(aes(yintercept=mean(df_clusters$n_genes)), 
             color="red",  size=0.6) +
  annotate("text", x=27, y=145, label="gene count mean")

df_clusters %>%
  ggplot(aes(x=cluster, y=n_genes)) +
  geom_line(aes(group=1)) +
  geom_point(shape=21, color="black", fill="#69b3a2", size=2) +
  labs(title="degPatterns genes per cluster") +
  xlab("Cluster") + 
  ylab("Number of genes per cluster") +
  geom_hline(aes(yintercept=mean(df_clusters$n_genes)),
             color="red", size=0.6)+
  annotate("text", x=27, y=145, label="gene count mean")

# Cutoff of genes per cluster in 200
clusters_1 <- degPatterns(cluster_rlog, metadata = meta, time="hours", col="index", minc = 200, scale = TRUE)

degPlotCluster(clusters_1$normalized, "timepoint", "timepoint", lines = TRUE, points = FALSE, boxes = TRUE, facet = TRUE, smooth = TRUE)

clusters_1_df <- clusters_1$df

cols_clust <- c("0h" = "#F8766D", "12h" = "#B79F00", "24h" = "#00BA38", "72h" = "#00BFC4", "7d" = "#619CFF", "30d" = "#F564E3")

ggplot(clusters[["normalized"]],
       aes(timepoint, value, color = timepoint, fill = timepoint)) +
  geom_boxplot() +
  geom_point(position = position_jitterdodge(dodge.width = 0.9)) +
  scale_color_manual(values = cols_clust, aesthetics = c("colour", "fill")) +
  scale_x_discrete(limits=c("0h","12h","24h","72h","7d","30d")) +
  geom_smooth(aes(group=timepoint), method = "lm")

# Explore the output:
head(clusters$df)

# After extracting a group of genes, we can perform functional analysis to explore associated functions... --> go to next script "functional_analysis.R"

###     5.2. maSigPro ####
maSigProUsersGuide()

# load data
meta_maSigPro <- read.table("meta/meta_MaSigPro", header = TRUE, row.names = 1)
meta_maSigPro_log2 <- read.table("meta/meta_MaSigPro_1group", header = TRUE, row.names = 1)

meta_maSig <- make.design.matrix(meta_maSigPro)
meta_maSig_log2 <- make.design.matrix(meta_maSigPro_log2)

# If we can use maSigPro with theta = 10:
pos_norm <- which(row.names(normalized_counts) %in% sig_lrt_0.01$gene)
sig_norm <- normalized_counts[pos_norm,]

normp_log2 <- p.vector(sig_norm, meta_maSig_log2, counts = TRUE)
normt_log2 <- T.fit(normp_log2)
get_norm_1group <- get.siggenes(normt_log2, vars = "all")
get_norm_1group$summary
# Changing this arguments in p.vector is enough. Further functions will take these arguments from a p.vector object.

# See clusters patterns
see_masig_norm <- see.genes(get_norm_1group$sig.genes, k = 7)
masig_clustertable <- as.data.frame(see_masig_norm$cut)

# Let's explore the set of genes in Group1 in more detail:
# Extract the Group 1 genes
post_clusters <- clusters$df
pattern1 <- clusters$df %>%
  dplyr::filter(cluster == 1)

###     5.3. Kmeans ###########
vst_DEGs_lrt <- read.table(file = "data/vst_DEGs_lrt01.xls", 
                           header = TRUE, row.names = 1) # retrive vst from data
vst_scaled <- t(scale(t(vst_DEGs_lrt))) # input kmeans
set.seed(123)
km.res <- kmeans(vst_scaled, centers = 3, iter.max = 10, nstart = 25) # centers = 3, 4, 5...
clustered_genes3 <- as.data.frame(km.res$cluster)
colnames(clustered_genes3)[1] ="cluster"
clustered_genes3$cluster <- paste('cluster', clustered_genes3$cluster, sep=' ')

# write.table(clustered_genes3, file = "data/clustered_genes3_01.xls", sep="\t", row.names = TRUE)

# elbow plot
fviz_nbclust(vst_scaled, kmeans, method = "wss") +
  geom_vline(xintercept = 3, linetype = 2)

# cluster plot
fviz_cluster(km.res, data = vst_scaled[, -5],
             palette = c("orange", "red", "purple", "blue", "green"),
             geom = "point",
             ellipse.type = "convex", 
             ggtheme = theme_bw()
)

# pdf("plots/kmeans0.01_elbowplot.pdf", width=12, height=8)
# fviz_nbclust(vst_scaled, kmeans, method = "wss") +
#   geom_vline(xintercept = 3, linetype = 2)
# dev.off()


# Retrieve kmeans gene clustering data
clust3_genes = read.table("data/clustered_genes3.xls", header = TRUE, row.names = 1)
clust4_genes = read.table("data/clustered_genes4.xls", header = TRUE, row.names = 1)
clust5_genes = read.table("data/clustered_genes5.xls", header = TRUE, row.names = 1)

# annotation colors
# annoCol3 <-list(cluster=c("cluster 1" ="orange", "cluster 2"="red", "cluster 3"="purple"), deg_Patt_cluster=c("patt 1"="blue", "patt 2"="orange", "patt 3"="yellow", "patt 4"="purple", "patt 5"="grey", "patt 7"="brown", "patt 9"="pink"))
annoCol3 <-list(cluster=c("cluster 1" ="orange", "cluster 2"="red", "cluster 3"="purple"))
annoCol4 <-list(cluster=c("cluster 1" ="orange", "cluster 2"="red", "cluster 3"="purple", "cluster 4"="blue"))
annoCol5 <-list(cluster=c("cluster 1" ="orange", "cluster 2"="red", "cluster 3"="purple", "cluster 4"="blue", "cluster 5" = "green"))

clusters_1_df = clusters_1_df %>% rename(degPatt_cluster = cluster)
intento1 <- merge(clust3_genes, clusters_1_df, by = "row.names")
intento1 <- intento1[-1]
row.names(intento1) <- intento1$genes
intento1 <- intento1[-2]
intento1$deg_Patt_cluster <- paste('patt', intento1$deg_Patt_cluster, sep=' ')

# clustered LRT 0.01 heatmap
pheatmap(vst_DEGs_lrt01, scale = "row", cluster_rows=TRUE, cluster_cols=FALSE, show_rownames = FALSE, annotation_row = clust3_genes, annotation_colors = annoCol3, gaps_row = 3, cutree_rows = 3)

pheatmap(vst_DEGs_lrt01, scale = "row", cluster_rows=TRUE, cluster_cols=FALSE, show_rownames = FALSE, annotation_row = clust4_genes, annotation_colors = annoCol4, cutree_rows = 4)

pheatmap(vst_DEGs_lrt01, scale = "row", cluster_rows=TRUE, cluster_cols=FALSE, show_rownames = FALSE, annotation_row = intento1, annotation_colors = annoCol5, cutree_rows = 5)

# ("plots/LRT_0.01_heatmap.", width=12, height=8)
pheatmap(vst_DEGs_lrt01, scale = "row", cluster_rows=TRUE, cluster_cols=FALSE, show_rownames = FALSE)
# dev.off()


###     5.4. Kmeans, degPatterns & MaSigpro concordance ####
# Run this code changing the cluster names to get all necessary data
kmeans_clust3 <- clust3_genes %>%
  dplyr::filter(cluster == "cluster 3")

kmeans5_clust5 <- clust5_genes %>%
  dplyr::filter(cluster == "cluster 5")

masig_clust7 <- masig_clustertable %>%
  dplyr::filter(see_masig_norm$cut == "7")

rownames(pattern1)

clust_venn3kmeans <- list("patt1" = rownames(pattern1), "patt2" = rownames(pattern2), "patt3" = rownames(pattern3), "patt4" = rownames(pattern4), "patt5" = rownames(pattern5), "patt7" = rownames(pattern7), "patt9" = rownames(pattern9), "k1" = rownames(kmeans_clust1), "k2" = rownames(kmeans_clust2), "k3" = rownames(kmeans_clust3), "mSP1" = rownames(masig_clust1), "mSP2" = rownames(masig_clust2),"mSP3" = rownames(masig_clust3),"mSP4" = rownames(masig_clust4),"mSP5" = rownames(masig_clust5),"mSP6" = rownames(masig_clust6),"mSP7" = rownames(masig_clust7))

clust_venn5kmeans <- list("patt1" = rownames(pattern1), "patt2" = rownames(pattern2), "patt3" = rownames(pattern3), "patt4" = rownames(pattern4), "patt5" = rownames(pattern5), "patt7" = rownames(pattern7), "patt9" = rownames(pattern9), "k1" = rownames(kmeans5_clust1), "k2" = rownames(kmeans5_clust2), "k3" = rownames(kmeans5_clust3), "k4" = rownames(kmeans5_clust4), "k5" = rownames(kmeans5_clust5))

plot(euler(clust_venn3kmeans), quantities = TRUE, legend = TRUE, labels = TRUE, adjust_labels = TRUE, main = "7 gene patterns & 3 clusters kmeans")

plot(euler(clust_venn5kmeans), quantities = TRUE, legend = TRUE, labels = TRUE, adjust_labels = TRUE, main = "gene patterns 5 clusters kmeans")

clust_venn_k1 <- list("patt1" = rownames(pattern1), "patt2" = rownames(pattern2), "patt3" = rownames(pattern3), "k2" = rownames(kmeans_clust2), "k3" = rownames(kmeans_clust3))

plot(venn(clust_venn_k1), quantities = TRUE, legend = TRUE, labels = TRUE, adjust_labels = TRUE, main = "gene patterns")


############## 6. Celltype deconvolution - CibersortX ####
# Get input mixture table
gene_symbols <- rownames(data)
data_ciberX <- data.frame(GeneSym = gene_symbols, data)

# write.table(data_ciberX, file='/home/sozaez/Desktop/TFM/scRNAseq/actualized_Olaya/scRNAseq_TFM_Olaya/deconvolution/bulk_rawdata_8K.tsv', quote=FALSE, sep='\t', col.names = TRUE, row.names = FALSE)
