setwd("/home/sozaez/Desktop/TFM/bulkRNAseq/Analyses/Uroth_reg_bulkRNAseq/")

# Load packages
library(grid)
library(pheatmap)
library(dslice)
library(dplyr)
library(tidyverse)
library(eulerr)

# Load data
vst_counts_df <- read.table("data/vst_counts.xls", header = TRUE, row.names = 1)
normalized_counts <- read.table("data/normalized_counts.txt", header = TRUE, row.names = 1)
rld_mat <- read.table("data/rld_mat.xls", header = TRUE, row.names = 1)
res_lrt_tb <- read.table("results/res_lrt_tb.xls", header = TRUE, row.names = 1)


############################ CHENG GEN SET ####################################
Cheng_genes_raw <- load_gmt("./data/GSEA/Datasets/Cheng_markers/Cheng_markers.gmt")
Cheng_genes <- Cheng_genes_raw[["gene_symbol"]][[1]]

LRT_Cheng <- res_lrt_tb[match(Cheng_genes, row.names(res_lrt_tb)),]

LRT_Cheng$celltype <- rep(c("basal", "intermediate","superficial","cycling"), 
             c(15, 17, 4, 4))


# Arrange genes by padj res_LRT within each celltype
LRT_Cheng <- arrange(arrange(LRT_Cheng, padj),celltype)
# Reorder celltypes
LRT_Cheng <- LRT_Cheng[c((1:cumsum(table(LRT_Cheng$celltype))[1]),
                         (cumsum(table(LRT_Cheng$celltype))[2]+1):cumsum(table(LRT_Cheng$celltype))[3],
                         (cumsum(table(LRT_Cheng$celltype))[3]+1):cumsum(table(LRT_Cheng$celltype))[4],
                         (cumsum(table(LRT_Cheng$celltype))[1]+1):cumsum(table(LRT_Cheng$celltype))[2]),]

LRT_Cheng <- slice(LRT_Cheng, c(1:7, 16:24, 33:34, 37:40))

# Get vst, normalized and rlog counts
vst_Cheng <- vst_counts_df[match(row.names(LRT_Cheng), row.names(vst_counts_df)),]
norm_Cheng <- normalized_counts[match(row.names(LRT_Cheng), row.names(normalized_counts)),]
rld_Cheng <- rld_mat[match(row.names(LRT_Cheng), row.names(rld_mat)),]

# substraction of mean 0h rlog
rlog_row_means <- rowMeans(rld_mat[,1:4])
rlog_0h_baseline <- rld_mat - rlog_row_means
rlog_0h_Cheng <- rlog_0h_baseline[match(row.names(LRT_Cheng), row.names(rlog_0h_baseline)),]

# Set pheatmap color parameters
paletteLength <- 50
myBreaks <- c(seq(min(rlog_0h_Cheng), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(rlog_0h_Cheng)/paletteLength, max(rlog_0h_Cheng), length.out=floor(paletteLength/2)))

myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)

cols <- RColorBrewer:::brewer.pal(9,"Reds")

# Annotation pheatmap
df_Cheng <- data.frame(celltype=LRT_Cheng$celltype, row.names = rownames(LRT_Cheng))

annoCheng <-list(celltype=c("basal" ="red", "intermediate"="green", "superficial"="orange", "cycling"="blue"))

# Plot heatmaps
pheatmap(norm_Cheng, cluster_rows=FALSE, cluster_cols=FALSE, scale = "row", color = myColor, annotation_row = df_Cheng, annotation_colors = annoCheng, main = "norm")

pheatmap(vst_Cheng, cluster_rows=FALSE, cluster_cols=FALSE, scale = "row", color = myColor, annotation_row = df_Cheng, annotation_colors = annoCheng, main = "vst")

pheatmap(rld_Cheng, cluster_rows=FALSE, cluster_cols=FALSE, scale = "row", color = myColor, annotation_row = df_Cheng, annotation_colors = annoCheng, main = "Cheng rlog - row scaled") # interesante...

pheatmap(rld_Cheng, cluster_rows=FALSE, cluster_cols=FALSE, color = cols, main = "Cheng rlog") # interesante... sin escalar

pheatmap(rlog_0h_Cheng, cluster_rows=FALSE, cluster_cols=FALSE, main = "Cheng genes - rlog - 0h mean subtracted", breaks = myBreaks, color = myColor, fontsize_row = 7, annotation_row = df_Cheng, annotation_colors = annoCheng, angle_col = 90, gaps_row = cumsum(table(LRT_Cheng$celltype)[c(1,3,4,2)]))


LRT_Cheng$row_colors <- ifelse(LRT_Cheng$padj < 0.05, "sig", "no_sig")

rlog_0h_Cheng_df <- data.frame(rlog_0h_Cheng)
rlog_0h_Cheng_df$status <- LRT_Cheng$row_colors
rlog_0h_Cheng_df$colors=ifelse(rlog_0h_Cheng_df$status=="sig","red","black")

Cheng_h = pheatmap(rlog_0h_Cheng, cluster_rows=FALSE, cluster_cols=FALSE, breaks = myBreaks, color = myColor, annotation_row = df_Cheng, annotation_colors = annoCheng, angle_col = 90, gaps_row = cumsum(table(LRT_Cheng$celltype)[c(1,3,4,2)]), cellwidth = 15, cellheight = 15, gaps_col = c(4,8,11,15,19), main = "Urothelial markers (Cheng et al. 2021)")

sig_cols=rlog_0h_Cheng_df[order(match(rownames(rlog_0h_Cheng_df), Cheng_h$gtable$grobs[[4]]$label)), ]$colors

Cheng_h$gtable$grobs[[4]]$gp=gpar(col=sig_cols)
Cheng_h

# pdf("plots/GSEA/Heatmaps/Cheng_rlog0h", width=12, height=8)
# Cheng_h
# dev.off()


###########################  REACTOME APOPTOSIS  ##############################
# Load .gmt file. Obtained from MSigDB Mus musculus - M2 Reactome
apoptosis_genes <-  load_gmt("./data/GSEA/Datasets/Reactome&hallmark_apoptosis/REACTOME_APOPTOSIS.v2022.1.Mm.gmt")
apoptosis_genes <- apoptosis_genes[["gene_symbol"]][[1]]

# Arrange genes by res LRT padj
LRT_apop <- res_lrt_tb[match(apoptosis_genes, row.names(res_lrt_tb)),]
LRT_apop <- arrange(LRT_apop, padj)

# Get vst, normalized and rlog counts
vst_apoptosis <- vst_counts_df[match(row.names(LRT_apop), row.names(vst_counts_df)),]
norm_apoptosis <- normalized_counts[match(row.names(LRT_apop), 
                                          row.names(normalized_counts)),]
rld_apoptosis <- rld_mat[match(row.names(LRT_apop), row.names(rld_mat)),]

# Subtract 0h counts mean
rlog_row_means <- rowMeans(rld_mat[,1:4])
rlog_0h_baseline <- rld_mat - rlog_row_means
rlog_0h_apoptosis <- rlog_0h_baseline[match(row.names(LRT_apop), row.names(rlog_0h_baseline)),]

# Set color parameters
paletteLength <- 50
br_0 <- c(seq(min(rlog_0h_apoptosis), 0, length.out=ceiling(paletteLength/2) + 1), 
          seq(max(rlog_0h_apoptosis)/paletteLength, max(rlog_0h_apoptosis), length.out=floor(paletteLength/2)))

br <- c(seq(-0.5, 0, length.out=ceiling(paletteLength/2) + 1), 
        seq(max(rld_apoptosis)/paletteLength, 0.5, length.out=floor(paletteLength/2)))

myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)

# Plot heatmaps
pheatmap(norm_apoptosis, cluster_rows=FALSE, cluster_cols=FALSE, color = myColor, show_rownames = FALSE, main = "norm apoptosis")
pheatmap(vst_apoptosis, cluster_rows=FALSE, cluster_cols=FALSE, color = myColor, show_rownames = FALSE, main = "vst apoptosis")
pheatmap(rld_apoptosis, cluster_rows=FALSE, cluster_cols=FALSE, color = myColor, show_rownames = FALSE, main = "rlog apoptosis")

pheatmap(rlog_0h_apoptosis, cluster_rows=FALSE, cluster_cols=FALSE, breaks = br_0, color = myColor, show_rownames = TRUE, main = "rlog0h apoptosis")


LRT_apop$row_colors <- ifelse(LRT_apop$padj < 0.05, "sig", "no_sig")

rlog_0h_apoptosis_df <- data.frame(rlog_0h_apoptosis)
rlog_0h_apoptosis_df$status <- LRT_apop$row_colors
rlog_0h_apoptosis_df$colors=ifelse(rlog_0h_apoptosis_df$status=="sig","red","black")

Apop_h = pheatmap(rlog_0h_apoptosis, cluster_rows=FALSE, cluster_cols=FALSE, breaks = myBreaks, color = myColor, angle_col = 90, main = "Apoptosis - rlog 0h mean subtracted")

sig_cols=rlog_0h_apoptosis_df[order(match(rownames(rlog_0h_apoptosis_df), Apop_h$gtable$grobs[[4]]$label)), ]$colors

Apop_h$gtable$grobs[[4]]$gp=gpar(col=sig_cols)
Apop_h

# pdf("plots/GSEA/Heatmaps/Cheng_rlog0h", width=12, height=8)
# Cheng_h
# dev.off()


###########################  REACTOME CELL CYCLE  #############################
# Set up list of canonical cell type markers
G1_G1_S_transition_list <-  load_gmt("./data/GSEA/Datasets/Reactome_cell_cycle/REACTOME_MITOTIC_G1_PHASE_AND_G1_S_TRANSITION.v2022.1.Mm.gmt")
S_list <-  load_gmt("./data/GSEA/Datasets/Reactome_cell_cycle/REACTOME_S_PHASE.v2022.1.Mm.gmt")
G2_and_G2_M_list <-  load_gmt("./data/GSEA/Datasets/Reactome_cell_cycle/REACTOME_MITOTIC_G2_G2_M_PHASES.v2022.1.Mm.gmt")
M_list <-  load_gmt("./data/GSEA/Datasets/Reactome_cell_cycle/REACTOME_M_PHASE.v2022.1.Mm.gmt")

G1_list <- setdiff(G1_G1_S_transition_list[["gene_symbol"]][[1]], 
                   S_list[["gene_symbol"]][[1]])
G2_list <- setdiff(G2_and_G2_M_list[["gene_symbol"]][[1]], 
                   M_list[["gene_symbol"]][[1]])

# create dataset
cell_cycle_genset <- c(G1_list, 
                       S_list[["gene_symbol"]][[1]],
                       G2_list,
                       M_list[["gene_symbol"]][[1]])

# Create list
cell_cycle_genset_list <- list("G1_list" = G1_list,
                               "S_phase" = S_list[["gene_symbol"]][[1]],
                               "G2_list" = G2_list,
                               "M_phase" = M_list[["gene_symbol"]][[1]])


# Create dataframe with phase column
phase <- rep(c("G1_phase","S_phase", "G2_phase", "M_phase"), 
             c(length(G1_list),
               length(S_list[["gene_symbol"]][[1]]), 
               length(G2_list),
               length(M_list[["gene_symbol"]][[1]])))

cecy_phase_df <- data.frame(phase)

# Order genes by res LRT padj within each phase
LRT_cell_cycle <- res_lrt_tb[match(cell_cycle_genset, row.names(res_lrt_tb)),]
LRT_cell_cycle <- cbind(cecy_phase_df, LRT_cell_cycle) 
LRT_cell_cycle$gene <- rownames(LRT_cell_cycle)
LRT_cell_cycle <- LRT_cell_cycle %>%  group_by(phase) %>% slice_min(order_by = padj, n = 15)
LRT_cell_cycle$gene <- as.character(LRT_cell_cycle$gene)
to_delete <- grepl("\\.[0-9]+$", LRT_cell_cycle$gene)
LRT_cell_cycle <- LRT_cell_cycle[!to_delete, ]

# LRT_cell_cycle$gene <- gsub("\\..*", "", LRT_cell_cycle$gene)


# Reorder phases
LRT_cell_cycle <- LRT_cell_cycle[c((1:cumsum(table(LRT_cell_cycle$phase))[1]),
                                   (cumsum(table(LRT_cell_cycle$phase))[3]+1):cumsum(table(LRT_cell_cycle$phase))[4],
                                   (cumsum(table(LRT_cell_cycle$phase))[1]+1):cumsum(table(LRT_cell_cycle$phase))[2],
                                   (cumsum(table(LRT_cell_cycle$phase))[2]+1):cumsum(table(LRT_cell_cycle$phase))[3]),]

# Select just first 8 genes per category
LRT_cell_cycle <- slice(LRT_cell_cycle, c(1:8, 16:24, 31:38, 40:48))


# In order to create the heatmap breaks
cell_cycle_genset_sub10 <- list("G1_phase" = LRT_cell_cycle[which(LRT_cell_cycle$phase == "G1_phase"),][["gene"]],
                                "S_phase" = LRT_cell_cycle[which(LRT_cell_cycle$phase == "S_phase"),][["gene"]],
                                "G2_phase" = LRT_cell_cycle[which(LRT_cell_cycle$phase == "G2_phase"),][["gene"]], 
                                "M_phase" = LRT_cell_cycle[which(LRT_cell_cycle$phase == "M_phase"),][["gene"]])


# Get vst, normalized and rlog counts
vst_cell_cycle <- vst_counts_df[match(LRT_cell_cycle$gene, row.names(vst_counts_df)),]
norm_cell_cycle <- normalized_counts[match(LRT_cell_cycle$gene, 
                                           row.names(normalized_counts)),]
rld_cell_cycle <- rld_mat[match(LRT_cell_cycle$gene, row.names(rld_mat)),]

# Subtraction of mean 0h rlog
rlog_row_means <- rowMeans(rld_mat[,1:4])
rlog_0h_baseline <- rld_mat - rlog_row_means
rlog_0h_cell_cycle <- rlog_0h_baseline[match(LRT_cell_cycle$gene, row.names(rlog_0h_baseline)),]
rlog_0h_cell_cycle$genes <- LRT_cell_cycle$gene

# row.names(rlog_0h_cell_cycle) <- ifelse(duplicated(rlog_0h_cell_cycle$genes), paste0((rlog_0h_cell_cycle$genes), "_", rep(c("G1", "S", "G2", "M"), c(10, 10, 10, 10))), rlog_0h_cell_cycle$genes)

rlog_0h_cell_cycle <- rlog_0h_cell_cycle[,-24]

# Set color parameters
paletteLength <- 50

br_0 <- c(seq(min(rlog_0h_cell_cycle), 0, length.out=ceiling(paletteLength/2) + 1),
          seq(max(rlog_0h_cell_cycle)/paletteLength, max(rlog_0h_cell_cycle), length.out=floor(paletteLength/2)))

myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)

# Heatmap annotations
df_Cell_cycle <- data.frame(phase=LRT_cell_cycle$phase, row.names = rownames(rlog_0h_cell_cycle))

annoCellCycle <-list(phase=c("G1_phase" = "red", "S_phase" = "orange","G2_phase" = "yellow","M_phase" = "green"))

# Cell cycle heatmaps
pheatmap(rlog_0h_cell_cycle, cluster_rows=FALSE, cluster_cols=FALSE, breaks = br_0, gaps_row = cumsum(summary(cell_cycle_genset_sub10)[,1]), color = myColor, annotation_row = df_Cell_cycle, annotation_colors = annoCellCycle, show_rownames = TRUE, angle_col = 90, number_color = 6 ,main = "rlog0h cell cycle")

# Mark significant genes (here are all significant)
LRT_cell_cycle$row_colors <- ifelse(LRT_cell_cycle$padj < 0.05, "sig", "no_sig")
rlog_0h_cell_cycle_df <- data.frame(rlog_0h_cell_cycle)
rlog_0h_cell_cycle_df$status <- LRT_cell_cycle$row_colors
rlog_0h_cell_cycle_df$colors=ifelse(rlog_0h_cell_cycle_df$status=="sig","red","black")

CellCycle_h = pheatmap(rlog_0h_cell_cycle, cluster_rows=FALSE, cluster_cols=FALSE, breaks = br_0, gaps_row = cumsum(summary(cell_cycle_genset_sub10)[,1]), color = myColor, annotation_row = df_Cell_cycle, annotation_colors = annoCellCycle, cellwidth = 10, cellheight = 10, gaps_col = c(4,8,11,15,19), show_rownames = TRUE, angle_col = 90, number_color = 6 ,main = "Reactome DB cell cycle markers")

sig_cols=rlog_0h_cell_cycle_df[order(match(rownames(rlog_0h_cell_cycle_df), CellCycle_h$gtable$grobs[[4]]$label)), ]$colors

CellCycle_h$gtable$grobs[[4]]$gp=gpar(col=sig_cols)
CellCycle_h

# pdf("plots/GSEA/Heatmaps/Cell_cycle_rlog0h", width=12, height=8)
# CellCycle_h
# dev.off()

###########################  DREAM COMPLEX  ##################################

# Set up list of canonical cell type markers
DREAM_rep_list <- load_gmt("./data/GSEA/Datasets/DREAM/DREAM_repressors.gmt")
DREAM_act_list <- load_gmt("./data/GSEA/Datasets/DREAM/DREAM_activators.gmt")
DREAM_targets_list <- load_gmt("./data/GSEA/Datasets/DREAM/maria_ramal_genesets.gmt")


# create dataset
DREAM_geneset <- c(DREAM_rep_list[["gene_symbol"]][[1]],
                   DREAM_act_list[["gene_symbol"]][[1]],
                   DREAM_targets_list[["gene_symbol"]][[1]])

# Create dataframe with geneset type column
DREAM_complex <- rep(c("DREAM_activators", "DREAM_repressors", "DREAM_targets"), 
                     c(length(DREAM_act_list[["gene_symbol"]][[1]]),
                       length(DREAM_rep_list[["gene_symbol"]][[1]]),
                       length(DREAM_targets_list[["gene_symbol"]][[1]])))

DREAM_complex_df <- data.frame(DREAM_complex)

# Order genes by res LRT padj within each phase
LRT_DREAM <- res_lrt_tb[match(DREAM_geneset, row.names(res_lrt_tb)),]
LRT_DREAM <- cbind(DREAM_complex_df, LRT_DREAM) 
LRT_DREAM$gene <- rownames(LRT_DREAM)
LRT_DREAM <- LRT_DREAM %>%  group_by(DREAM_complex) %>% slice_min(order_by = padj, n = 30)


# Get vst, normalized and rlog counts
vst_DREAM <- vst_counts_df[match(LRT_DREAM$gene, row.names(vst_counts_df)),]
norm_DREAM <- normalized_counts[match(LRT_DREAM$gene, 
                                      row.names(normalized_counts)),]
rld_DREAM <- rld_mat[match(LRT_DREAM$gene, row.names(rld_mat)),]

# Subtract 0h counts mean
rlog_row_means <- rowMeans(rld_mat[,1:4])
rlog_0h_baseline <- rld_mat - rlog_row_means
rlog_0h_DREAM <- rlog_0h_baseline[match(LRT_DREAM$gene, row.names(rlog_0h_baseline)),]

# In order to create the heatmap breaks
DREAM_genset_sub10 <- list("DREAM_activators" = LRT_DREAM[which(LRT_DREAM$DREAM_complex == "DREAM_activators"),][["gene"]],
                           "DREAM_repressors" = LRT_DREAM[which(LRT_DREAM$DREAM_complex == "DREAM_repressors"),][["gene"]],
                           "DREAM_targets" = LRT_DREAM[which(LRT_DREAM$DREAM_complex == "DREAM_targets"),][["gene"]])


# Set color parameters
paletteLength <- 50
br_0 <- c(seq(min(rlog_0h_DREAM), 0, length.out=ceiling(paletteLength/2) + 1), 
          seq(max(rlog_0h_DREAM)/paletteLength, max(rlog_0h_DREAM), length.out=floor(paletteLength/2)))

br <- c(seq(-4, 0, length.out=ceiling(paletteLength/2) + 1), 
        seq(max(rld_collagen)/paletteLength, 4, length.out=floor(paletteLength/2)))

myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)

# Heatmap annotations
df_DREAM <- data.frame(DREAM_complex=LRT_DREAM$DREAM_complex, row.names = rownames(rlog_0h_DREAM))

annoDREAM <-list(DREAM_complex=c("DREAM_activators" = "green", "DREAM_repressors" = "red", "DREAM_targets" = "blue"))


# Plot heatmaps
pheatmap(norm_DREAM, cluster_rows=FALSE, cluster_cols=FALSE, color = myColor, show_rownames = TRUE, main = "norm DREAM")
pheatmap(vst_DREAM, cluster_rows=FALSE, cluster_cols=FALSE, color = myColor, show_rownames = TRUE, main = "vst DREAM")
pheatmap(rld_DREAM, cluster_rows=FALSE, cluster_cols=FALSE, color = myColor, show_rownames = TRUE, main = "rlog DREAM")

pheatmap(rlog_0h_DREAM, cluster_rows=FALSE, cluster_cols=FALSE, breaks = br_0, gaps_row = cumsum(summary(DREAM_genset_sub10)[,1]), color = myColor, show_rownames = TRUE, annotation_row = df_DREAM, annotation_colors = annoDREAM, main = "rlog0h DREAM")

pheatmap(rlog_0h_DREAM, cluster_rows=FALSE, cluster_cols=FALSE, breaks = br, gaps_row = cumsum(summary(DREAM_genset_sub10)[,1]), color = myColor, show_rownames = TRUE, annotation_row = df_DREAM, annotation_colors = annoDREAM, main = "rlog0h DREAM")


LRT_DREAM$row_colors <- ifelse(LRT_DREAM$padj < 0.05, "sig", "no_sig")

rlog_0h_DREAM_df <- data.frame(rlog_0h_DREAM)
rlog_0h_DREAM_df$status <- LRT_DREAM$row_colors
rlog_0h_DREAM_df$colors=ifelse(rlog_0h_DREAM_df$status=="sig","red","black")

DREAM_h = pheatmap(rlog_0h_DREAM, cluster_rows=FALSE, cluster_cols=FALSE, breaks = br_0, gaps_row = cumsum(summary(DREAM_genset_sub10)[,1]), color = myColor, show_rownames = TRUE, annotation_row = df_DREAM, annotation_colors = annoDREAM, main = "DREAM")

sig_cols=rlog_0h_DREAM_df[order(match(rownames(rlog_0h_DREAM_df), DREAM_h$gtable$grobs[[4]]$label)), ]$colors

DREAM_h$gtable$grobs[[4]]$gp=gpar(col=sig_cols)
DREAM_h

pdf("plots/GSEA/Heatmaps/DREAM_rlog0h", width=12, height=8)
DREAM_h
dev.off()


###########################  REACTOME COLLAGENS  #############################
# Set up list of canonical cell type markers
collagen_assembly <-  load_gmt("./data/GSEA/Datasets/Collagens/REACTOME_ASSEMBLY_OF_COLLAGEN_FIBRILS_AND_OTHER_MULTIMERIC_STRUCTURES.v2022.1.Mm.gmt")
collagen_biosynthesis <-  load_gmt("./data/GSEA/Datasets/Collagens/REACTOME_COLLAGEN_BIOSYNTHESIS_AND_MODIFYING_ENZYMES.v2022.1.Mm.gmt")
collagen_trimerization <-  load_gmt("./data/GSEA/Datasets/Collagens/REACTOME_COLLAGEN_CHAIN_TRIMERIZATION.v2022.1.Mm.gmt")
collagen_formation <-  load_gmt("./data/GSEA/Datasets/Collagens/REACTOME_COLLAGEN_FORMATION.v2022.1.Mm.gmt")

# create dataset
collagen_geneset <- c(collagen_assembly[["gene_symbol"]][[1]],
                       collagen_biosynthesis[["gene_symbol"]][[1]],
                       collagen_trimerization[["gene_symbol"]][[1]],
                       collagen_formation[["gene_symbol"]][[1]])

collagen_geneset <- unique(collagen_geneset)

# Arrange genes by res LRT padj
LRT_collagen <- res_lrt_tb[match(collagen_geneset, row.names(res_lrt_tb)),]
LRT_collagen <- arrange(LRT_collagen, padj)
LRT_collagen <- LRT_collagen %>% slice_min(order_by = padj, n =16)

# Get vst, normalized and rlog counts
vst_collagen <- vst_counts_df[match(row.names(LRT_collagen), row.names(vst_counts_df)),]
norm_collagen <- normalized_counts[match(row.names(LRT_collagen), 
                                          row.names(normalized_counts)),]
rld_collagen <- rld_mat[match(row.names(LRT_collagen), row.names(rld_mat)),]

# Subtract 0h counts mean
rlog_row_means <- rowMeans(rld_mat[,1:4])
rlog_0h_baseline <- rld_mat - rlog_row_means
rlog_0h_collagen <- rlog_0h_baseline[match(row.names(LRT_collagen), row.names(rlog_0h_baseline)),]

# Set color parameters
paletteLength <- 50
br_0 <- c(seq(min(rlog_0h_collagen), 0, length.out=ceiling(paletteLength/2) + 1), 
          seq(max(rlog_0h_collagen)/paletteLength, max(rlog_0h_collagen), length.out=floor(paletteLength/2)))

br <- c(seq(-4, 0, length.out=ceiling(paletteLength/2) + 1), 
        seq(max(rld_collagen)/paletteLength, 4, length.out=floor(paletteLength/2)))

myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)

# Plot heatmaps
pheatmap(norm_collagen, cluster_rows=FALSE, cluster_cols=FALSE, color = myColor, show_rownames = FALSE, main = "norm collagen")
pheatmap(vst_collagen, cluster_rows=FALSE, cluster_cols=FALSE, color = myColor, show_rownames = FALSE, main = "vst collagen")
pheatmap(rld_collagen, cluster_rows=FALSE, cluster_cols=FALSE, color = myColor, show_rownames = FALSE, main = "rlog collagen")

pheatmap(rlog_0h_collagen, cluster_rows=FALSE, cluster_cols=FALSE, breaks = br_0, color = myColor, show_rownames = TRUE, main = "rlog0h collagen")


LRT_collagen$row_colors <- ifelse(LRT_collagen$padj < 0.05, "sig", "no_sig")

rlog_0h_collagen_df <- data.frame(rlog_0h_collagen)
rlog_0h_collagen_df$status <- LRT_collagen$row_colors
rlog_0h_collagen_df$colors=ifelse(rlog_0h_collagen_df$status=="sig","red","black")

Collagen_h = pheatmap(rlog_0h_collagen, cluster_rows=FALSE, cluster_cols=FALSE, breaks = br_0, cellwidth = 20, cellheight = 20, gaps_col = c(4,8,11,15,19), angle_col = 45, color = myColor, show_rownames = TRUE, main = "Reactome BD collagens")

sig_cols=rlog_0h_collagen_df[order(match(rownames(rlog_0h_apoptosis_df), Collagen_h$gtable$grobs[[4]]$label)), ]$colors

Collagen_h$gtable$grobs[[4]]$gp=gpar(col=sig_cols)
Collagen_h

# pdf("plots/GSEA/Heatmaps/collagens_rlog0h", width=12, height=8)
# Collagen_h
# dev.off()


###########################  HEME #############################################
# Set up list of blood genesets: reactome, hallmarks and wikipathways
H_coagulation_list <-  load_gmt("./data/GSEA/Datasets/Reactome&hallmark_blood/HALLMARK_COAGULATION.v2022.1.Mm.gmt")
H_Heme_metabolism_list <-  load_gmt("./data/GSEA/Datasets/Reactome&hallmark_blood/HALLMARK_HEME_METABOLISM.v2022.1.Mm.gmt")
R_blood_biosynthesis_list <-  load_gmt("./data/GSEA/Datasets/Reactome&hallmark_blood/REACTOME_BLOOD_GROUP_SYSTEMS_BIOSYNTHESIS.v2022.1.Mm.gmt")
R_heme_degradation_list <-  load_gmt("./data/GSEA/Datasets/Reactome&hallmark_blood/REACTOME_HEME_DEGRADATION.v2022.1.Mm.gmt")
# WP_heme_biosynthesis_list <-  load_gmt("./data/GSEA/Datasets/Reactome&hallmark_blood/WP_HEME_BIOSYNTHESIS.v2022.1.Mm.gmt")


# create dataset
heme_genset <- c(H_coagulation_list[["gene_symbol"]][[1]], 
                 H_Heme_metabolism_list[["gene_symbol"]][[1]],
                 R_blood_biosynthesis_list[["gene_symbol"]][[1]])

# Create list
heme_geneset_list <- list("H_coagulation" = H_coagulation_list[["gene_symbol"]][[1]],
                               "H_heme_metabolism" = H_Heme_metabolism_list[["gene_symbol"]][[1]],
                               "R_blood_biosynthesis" = R_blood_biosynthesis_list[["gene_symbol"]][[1]])


# Create dataframe with db_source column
db_source <- rep(c("H_coagulation","H_heme_metabolism", "R_blood_biosynthesis"), 
             c(length(H_coagulation_list[["gene_symbol"]][[1]]),
               length(H_Heme_metabolism_list[["gene_symbol"]][[1]]),
               length(R_blood_biosynthesis_list[["gene_symbol"]][[1]])))

heme_db_source_df <- data.frame(db_source)

# Order genes by res LRT padj within each db_source
LRT_heme <- res_lrt_tb[match(heme_genset, row.names(res_lrt_tb)),]
LRT_heme <- cbind(heme_db_source_df, LRT_heme) 
LRT_heme$gene <- rownames(LRT_heme)
LRT_heme <- LRT_heme %>%  group_by(db_source) %>% slice_min(order_by = padj, n = 10)


# In order to create the heatmap breaks
heme_sub10 <- list("H_coagulation" = LRT_heme[which(LRT_heme$db_source == "H_coagulation"),][["gene"]],
                                "H_heme_metabolism" = LRT_heme[which(LRT_heme$db_source == "H_heme_metabolism"),][["gene"]],
                                "R_blood_biosynthesis" = LRT_heme[which(LRT_heme$db_source == "R_blood_biosynthesis"),][["gene"]])


# Get vst, normalized and rlog counts
vst_heme <- vst_counts_df[match(LRT_heme$gene, row.names(vst_counts_df)),]
norm_heme <- normalized_counts[match(LRT_heme$gene, 
                                           row.names(normalized_counts)),]
rld_heme <- rld_mat[match(LRT_heme$gene, row.names(rld_mat)),]

# Subtraction of mean 0h rlog
rlog_row_means <- rowMeans(rld_mat[,1:4])
rlog_0h_baseline <- rld_mat - rlog_row_means
rlog_0h_heme <- rlog_0h_baseline[match(LRT_heme$gene, row.names(rlog_0h_baseline)),]
rlog_0h_heme$genes <- LRT_heme$gene
rlog_0h_heme <- rlog_0h_heme[-24]

# Set color parameters
paletteLength <- 50
br_0 <- c(seq(min(rlog_0h_heme), 0, length.out=ceiling(paletteLength/2) + 1), 
          seq(max(rlog_0h_heme)/paletteLength, max(rlog_0h_heme), length.out=floor(paletteLength/2)))

myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)

# Heatmap annotations
df_heme <- data.frame(db_source=LRT_heme$db_source, row.names = rownames(rlog_0h_heme))

annoheme <-list(db_source=c("H_coagulation" = "red", "H_heme_metabolism" = "orange","R_blood_biosynthesis" = "yellow"))

# Plot heatmaps
pheatmap(rlog_0h_heme, cluster_rows=FALSE, cluster_cols=FALSE, breaks = br_0, gaps_row = cumsum(summary(heme_sub10)[,1]), color = myColor, show_rownames = TRUE, angle_col = 90, annotation_row = df_heme, annotation_colors = annoheme, number_color = 6, main = "rlog0h heme")


# Mark significant genes (here are all significant)
LRT_heme$row_colors <- ifelse(LRT_heme$padj < 0.05, "sig", "no_sig")
rlog_0h_heme_df <- data.frame(rlog_0h_heme)
rlog_0h_heme_df$status <- LRT_heme$row_colors
rlog_0h_heme_df$colors=ifelse(rlog_0h_heme_df$status=="sig","red","black")

Heme_h = pheatmap(rlog_0h_heme, cluster_rows=FALSE, cluster_cols=FALSE, breaks = br_0, gaps_row = cumsum(summary(heme_sub10)[,1]), color = myColor, show_rownames = TRUE, angle_col = 90, annotation_row = df_heme, annotation_colors = annoheme, number_color = 6, main = "Blood pathways")

sig_cols=rlog_0h_heme_df[order(match(rownames(rlog_0h_heme_df), Heme_h$gtable$grobs[[4]]$label)), ]$colors

Heme_h$gtable$grobs[[4]]$gp=gpar(col=sig_cols)
Heme_h

pdf("plots/GSEA/Heatmaps/Heme_rlog0h", width=12, height=8)
Heme_h
dev.off()

###########################  REACTOME KERATINITATION  ##############################
# Load keratinitation geneset
keratinitation_genes_list <-  load_gmt("./data/GSEA/Datasets/Keratinitation/REACTOME_KERATINIZATION.v2022.1.Mm.gmt")

kera_genset <- keratinitation_genes_list[["gene_symbol"]][[1]]

# Arrange genes by res LRT padj
LRT_kera <- res_lrt_tb[match(kera_genset, row.names(res_lrt_tb)),]
LRT_kera <- arrange(LRT_kera, padj)
LRT_kera <- LRT_kera %>% slice_min(order_by = padj, n = 30)


# Get vst, normalized and rlog counts
vst_kera <- vst_counts_df[match(row.names(LRT_kera), row.names(vst_counts_df)),]
norm_kera <- normalized_counts[match(row.names(LRT_kera), 
                                          row.names(normalized_counts)),]
rld_kera <- rld_mat[match(row.names(LRT_kera), row.names(rld_mat)),]

# Subtract 0h counts mean
rlog_row_means <- rowMeans(rld_mat[,1:4])
rlog_0h_baseline <- rld_mat - rlog_row_means
rlog_0h_kera <- rlog_0h_baseline[match(row.names(LRT_kera), row.names(rlog_0h_baseline)),]

# Set color parameters
paletteLength <- 50
br_0 <- c(seq(min(rlog_0h_kera), 0, length.out=ceiling(paletteLength/2) + 1), 
          seq(max(rlog_0h_kera)/paletteLength, max(rlog_0h_kera), length.out=floor(paletteLength/2)))

br <- c(seq(-0.5, 0, length.out=ceiling(paletteLength/2) + 1), 
        seq(max(rld_apoptosis)/paletteLength, 0.5, length.out=floor(paletteLength/2)))

myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)

# Plot heatmaps
pheatmap(norm_kera, cluster_rows=FALSE, cluster_cols=FALSE, color = myColor, show_rownames = FALSE, main = "norm kera")
pheatmap(vst_kera, cluster_rows=FALSE, cluster_cols=FALSE, color = myColor, show_rownames = FALSE, main = "vst kera")
pheatmap(rld_kera, cluster_rows=FALSE, cluster_cols=FALSE, color = myColor, show_rownames = FALSE, main = "rlog kera")

pheatmap(rlog_0h_kera, cluster_rows=FALSE, cluster_cols=FALSE, breaks = br_0, color = myColor, show_rownames = TRUE, main = "rlog0h kera")


LRT_kera$row_colors <- ifelse(LRT_apop$padj < 0.05, "sig", "no_sig")

rlog_0h_kera_df <- data.frame(rlog_0h_apoptosis)
rlog_0h_kera_df$status <- LRT_apop$row_colors
rlog_0h_kera_df$colors=ifelse(rlog_0h_kera_df$status=="sig","red","black")

Kera_h = pheatmap(rlog_0h_kera, cluster_rows=FALSE, cluster_cols=FALSE, breaks = br_0, color = myColor, show_rownames = TRUE, main = "Keratinitation")

sig_cols=rlog_0h_kera_df[order(match(rownames(rlog_0h_kera_df), Kera_h$gtable$grobs[[4]]$label)), ]$colors

Kera_h$gtable$grobs[[4]]$gp=gpar(col=sig_cols)
Kera_h

# pdf("plots/GSEA/Heatmaps/Kera_rlog0h", width=12, height=8)
# Kera_h
# dev.off()


###########################  INFLAMMATORY RESPONSE  ##############################
# Load inflammatory response geneset
R_cytokine_signaling_list <-  load_gmt("./data/GSEA/Datasets/Inflammation/REACTOME_CYTOKINE_SIGNALING_IN_IMMUNE_SYSTEM.v2022.1.Mm.gmt")

H_inflammatory_response_list <-  load_gmt("./data/GSEA/Datasets/Inflammation/HALLMARK_INFLAMMATORY_RESPONSE.v2022.1.Mm.gmt")

WP_inflammatory_response_list <-  load_gmt("./data/GSEA/Datasets/Inflammation/WP_INFLAMMATORY_RESPONSE_PATHWAY.v2022.1.Mm.gmt")

WP_cytokines_list <-  load_gmt("./data/GSEA/Datasets/Inflammation/WP_CYTOKINES_AND_INFLAMMATORY_RESPONSE.v2022.1.Mm.gmt")

# create dataset
inflam_genset <- c(R_cytokine_signaling_list[["gene_symbol"]][[1]], 
                   H_inflammatory_response_list[["gene_symbol"]][[1]],
                   WP_inflammatory_response_list[["gene_symbol"]][[1]],
                   WP_cytokines_list[["gene_symbol"]][[1]])

# Create list
inflam_geneset_list <- list("R_cytokine_signaling" = R_cytokine_signaling_list[["gene_symbol"]][[1]],
                          "H_inflammatory_response" = H_inflammatory_response_list[["gene_symbol"]][[1]],
                          "WP_inflammatory_response" = WP_inflammatory_response_list[["gene_symbol"]][[1]],
                          "WP_cytokines" =               WP_cytokines_list[["gene_symbol"]][[1]])

plot(euler(inflam_geneset_list), quantities = TRUE, legend = TRUE, labels = TRUE, adjust_labels = TRUE, main = "Reactome vs Seurat")

# Create dataframe with db_source column
db_source <- rep(c("R_cytokine_signaling","H_inflammatory_response", "WP_inflammatory_response", "WP_cytokines"), 
                 c(length(R_cytokine_signaling_list[["gene_symbol"]][[1]]),
                   length(H_inflammatory_response_list[["gene_symbol"]][[1]]),
                   length(WP_inflammatory_response_list[["gene_symbol"]][[1]]),
                   length(WP_cytokines_list[["gene_symbol"]][[1]])))

inflam_db_source_df <- data.frame(db_source)

# Order genes by res LRT padj within each db_source
LRT_inflam <- res_lrt_tb[match(inflam_genset, row.names(res_lrt_tb)),]
LRT_inflam <- cbind(inflam_db_source_df, LRT_inflam) 
LRT_inflam$gene <- rownames(LRT_inflam)
LRT_inflam <- LRT_inflam %>%  group_by(db_source) %>% slice_min(order_by = padj, n = 10)
LRT_inflam$gene <- as.character(LRT_inflam$gene)
to_delete <- grepl("\\.[0-9]+$", LRT_inflam$gene)
LRT_inflam <- LRT_inflam[!to_delete, ]


# In order to create the heatmap breaks
inflam_sub10 <- list("H_inflammatory_response" = LRT_inflam[which(LRT_inflam$db_source == "H_inflammatory_response"),][["gene"]],
                     "R_cytokine_signaling" = LRT_inflam[which(LRT_inflam$db_source == "R_cytokine_signaling"),][["gene"]],
                     "WP_cytokines" = LRT_inflam[which(LRT_inflam$db_source == "WP_cytokines"),][["gene"]],
                     "WP_inflammatory_response" = LRT_inflam[which(LRT_inflam$db_source == "WP_inflammatory_response"),][["gene"]])

# Fix
db_source <- rep(c("H_inflammatory_response", "R_cytokine_signaling", "WP_cytokines", "WP_inflammatory_response"), 
                 c(length(LRT_inflam[which(LRT_inflam$db_source == "H_inflammatory_response"),][["gene"]]),
                   length(LRT_inflam[which(LRT_inflam$db_source == "R_cytokine_signaling"),][["gene"]]),
                   length(LRT_inflam[which(LRT_inflam$db_source == "WP_cytokines"),][["gene"]]),
                   length(LRT_inflam[which(LRT_inflam$db_source == "WP_inflammatory_response"),][["gene"]])))

# Get vst, normalized and rlog counts
vst_inflam <- vst_counts_df[match(LRT_inflam$gene, row.names(vst_counts_df)),]
norm_inflam <- normalized_counts[match(LRT_inflam$gene, 
                                     row.names(normalized_counts)),]
rld_inflam <- rld_mat[match(LRT_inflam$gene, row.names(rld_mat)),]

# Subtraction of mean 0h rlog
rlog_row_means <- rowMeans(rld_mat[,1:4])
rlog_0h_baseline <- rld_mat - rlog_row_means
rlog_0h_inflam <- rlog_0h_baseline[match(LRT_inflam$gene, row.names(rlog_0h_baseline)),]
rlog_0h_inflam$genes <- LRT_inflam$gene
rlog_0h_inflam <- rlog_0h_inflam[-24]

# Set color parameters
paletteLength <- 50
br_0 <- c(seq(min(rlog_0h_inflam), 0, length.out=ceiling(paletteLength/2) + 1), 
          seq(max(rlog_0h_inflam)/paletteLength, max(rlog_0h_inflam), length.out=floor(paletteLength/2)))

myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)

# Heatmap annotations
df_inflam <- data.frame(db_source=db_source, row.names = rownames(rlog_0h_inflam))

annoinflam <-list(db_source=c("H_inflammatory_response" = "red", "R_cytokine_signaling" = "orange", "WP_cytokines" = "yellow", "WP_inflammatory_response" = "green"))

# Plot heatmaps
pheatmap(rlog_0h_inflam, cluster_rows=FALSE, cluster_cols=FALSE, breaks = br_0, gaps_row = cumsum(summary(inflam_sub10)[,1]), color = myColor, show_rownames = TRUE, angle_col = 90, annotation_row = df_inflam, annotation_colors = annoinflam, number_color = 6, main = "rlog0h inflammation")

# Mark significant genes (here are all significant)
LRT_inflam$row_colors <- ifelse(LRT_inflam$padj < 0.05, "sig", "no_sig")
rlog_0h_inflam_df <- data.frame(rlog_0h_inflam)
rlog_0h_inflam_df$status <- LRT_inflam$row_colors
rlog_0h_inflam_df$colors=ifelse(rlog_0h_inflam_df$status=="sig","red","black")

Inflam_h = pheatmap(rlog_0h_inflam, cluster_rows=FALSE, cluster_cols=FALSE, breaks = br_0, gaps_row = cumsum(summary(inflam_sub10)[,1]), color = myColor, show_rownames = TRUE, angle_col = 90, annotation_row = df_inflam, annotation_colors = annoinflam, number_color = 6, main = "Inflammation pathways")

sig_cols=rlog_0h_inflam_df[order(match(rownames(rlog_0h_inflam_df), Inflam_h$gtable$grobs[[4]]$label)), ]$colors

Inflam_h$gtable$grobs[[4]]$gp=gpar(col=sig_cols)
Inflam_h

pdf("plots/GSEA/Heatmaps/Inflam_rlog0h", width=12, height=8)
Inflam_h
dev.off()


###########################  REACTOME COLLAGENS  #############################
# Set up list of canonical cell type markers
collagen_assembly <-  load_gmt("./data/GSEA/Datasets/Immunity/REACTOME_DOWNSTREAM_SIGNALING_EVENTS_OF_B_CELL_RECEPTOR_BCR.v2022.1.Mm.gmt")
collagen_biosynthesis <-  load_gmt("./data/GSEA/Datasets/Immunity/REACTOME_SIGNALING_BY_THE_B_CELL_RECEPTOR_BCR.v2022.1.Mm.gmt")
collagen_trimerization <-  load_gmt("./data/GSEA/Datasets/Immunity/REACTOME_IMMUNOREGULATORY_INTERACTIONS_BETWEEN_A_LYMPHOID_AND_A_NON_LYMPHOID_CELL.v2022.1.Mm.gmt")
collagen_formation <-  load_gmt("./data/GSEA/Datasets/Immunity/REACTOME_CHEMOKINE_RECEPTORS_BIND_CHEMOKINES.v2022.1.Mm.gmt")
collagen_formation <-  load_gmt("./data/GSEA/Datasets/Immunity/REACTOME_DEFENSINS.v2022.1.Mm.gmt")
collagen_formation <-  load_gmt("./data/GSEA/Datasets/Immunity/REACTOME_TCR_SIGNALING.v2022.1.Mm.gmt")
collagen_formation <-  load_gmt("./data/GSEA/Datasets/Immunity/REACTOME_BINDING_AND_UPTAKE_OF_LIGANDS_BY_SCAVENGER_RECEPTORS.v2022.1.Mm.gmt")

# create dataset
collagen_geneset <- c(collagen_assembly[["gene_symbol"]][[1]],
                      collagen_biosynthesis[["gene_symbol"]][[1]],
                      collagen_trimerization[["gene_symbol"]][[1]],
                      collagen_formation[["gene_symbol"]][[1]])

collagen_geneset <- unique(collagen_geneset)

# Arrange genes by res LRT padj
LRT_collagen <- res_lrt_tb[match(collagen_geneset, row.names(res_lrt_tb)),]
LRT_collagen <- arrange(LRT_collagen, padj)
LRT_collagen <- LRT_collagen %>% slice_min(order_by = padj, n =30)

# Get vst, normalized and rlog counts
vst_collagen <- vst_counts_df[match(row.names(LRT_collagen), row.names(vst_counts_df)),]
norm_collagen <- normalized_counts[match(row.names(LRT_collagen), 
                                         row.names(normalized_counts)),]
rld_collagen <- rld_mat[match(row.names(LRT_collagen), row.names(rld_mat)),]

# Subtract 0h counts mean
rlog_row_means <- rowMeans(rld_mat[,1:4])
rlog_0h_baseline <- rld_mat - rlog_row_means
rlog_0h_collagen <- rlog_0h_baseline[match(row.names(LRT_collagen), row.names(rlog_0h_baseline)),]

# Set color parameters
paletteLength <- 50
br_0 <- c(seq(min(rlog_0h_collagen), 0, length.out=ceiling(paletteLength/2) + 1), 
          seq(max(rlog_0h_collagen)/paletteLength, max(rlog_0h_collagen), length.out=floor(paletteLength/2)))

br <- c(seq(-4, 0, length.out=ceiling(paletteLength/2) + 1), 
        seq(max(rld_collagen)/paletteLength, 4, length.out=floor(paletteLength/2)))

myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)

# Plot heatmaps
pheatmap(norm_collagen, cluster_rows=FALSE, cluster_cols=FALSE, color = myColor, show_rownames = FALSE, main = "norm collagen")
pheatmap(vst_collagen, cluster_rows=FALSE, cluster_cols=FALSE, color = myColor, show_rownames = FALSE, main = "vst collagen")
pheatmap(rld_collagen, cluster_rows=FALSE, cluster_cols=FALSE, color = myColor, show_rownames = FALSE, main = "rlog collagen")

pheatmap(rlog_0h_collagen, cluster_rows=FALSE, cluster_cols=FALSE, breaks = br_0, color = myColor, show_rownames = TRUE, main = "rlog0h collagen")


LRT_collagen$row_colors <- ifelse(LRT_collagen$padj < 0.05, "sig", "no_sig")

rlog_0h_collagen_df <- data.frame(rlog_0h_collagen)
rlog_0h_collagen_df$status <- LRT_collagen$row_colors
rlog_0h_collagen_df$colors=ifelse(rlog_0h_collagen_df$status=="sig","red","black")

Collagen_h = pheatmap(rlog_0h_collagen, cluster_rows=FALSE, cluster_cols=FALSE, breaks = br_0, color = myColor, show_rownames = TRUE, main = "Collagens")

sig_cols=rlog_0h_collagen_df[order(match(rownames(rlog_0h_apoptosis_df), Collagen_h$gtable$grobs[[4]]$label)), ]$colors

Collagen_h$gtable$grobs[[4]]$gp=gpar(col=sig_cols)
Collagen_h

pdf("plots/GSEA/Heatmaps/collagens_rlog0h", width=12, height=8)
Collagen_h
dev.off()





# Euler diagram - Cell cycle genesets: Reactome vs Seurat
Seurat_cell_cycle <-  load_gmt("./data/GSEA/Datasets/Seurat_cell_cycle/g2m_s.gmt")

Seurat_list <- list("g2m" = Seurat_cell_cycle[["gene_symbol"]][[1]], "s" = Seurat_cell_cycle[["gene_symbol"]][[3]])

list_Reactome_Seurat <- list("R_G0_G1" = cell_cycle_genset_list$G0_early_G1, 
                             "R_G1" = cell_cycle_genset_list$G1_phase, 
                             "R_G1_S" = cell_cycle_genset_list$G1_S_transition, 
                             "R_S" = cell_cycle_genset_list$S_phase, 
                             "R_G2" = cell_cycle_genset_list$G2_phase, 
                             "R_M" = cell_cycle_genset_list$M_phase, 
                             "S_G2_M" = Seurat_list$g2m)

plot(venn(list_Reactome_Seurat), quantities = TRUE, legend = TRUE, labels = TRUE, adjust_labels = TRUE, main = "Reactome vs Seurat")

plot(euler(list_Reactome_Seurat), quantities = TRUE, legend = TRUE, labels = TRUE, adjust_labels = TRUE, main = "Reactome vs Seurat")



###########################  REACTOME NEUTROPHIL DEGRANULATION  ##############################
# Load neutrophil geneset
neutrophil_degranulation_genes_list <-  load_gmt("./data/GSEA/Datasets/Immunity/REACTOME_NEUTROPHIL_DEGRANULATION.v2022.1.Mm.gmt")

neutr_genset <- neutrophil_degranulation_genes_list[["gene_symbol"]][[1]]

# Arrange genes by res LRT padj
LRT_neutr <- res_lrt_tb[match(neutr_genset, row.names(res_lrt_tb)),]
LRT_neutr <- arrange(LRT_neutr, padj)
LRT_neutr <- LRT_neutr %>% slice_min(order_by = padj, n = 21)

# Get vst, normalized and rlog counts
vst_neutr <- vst_counts_df[match(row.names(LRT_neutr), row.names(vst_counts_df)),]
norm_neutr <- normalized_counts[match(row.names(LRT_neutr), 
                                     row.names(normalized_counts)),]
rld_neutr <- rld_mat[match(row.names(LRT_neutr), row.names(rld_mat)),]

# Subtract 0h counts mean
rlog_row_means <- rowMeans(rld_mat[,1:4])
rlog_0h_baseline <- rld_mat - rlog_row_means
rlog_0h_neutr <- rlog_0h_baseline[match(row.names(LRT_neutr), row.names(rlog_0h_baseline)),]

# Set color parameters
paletteLength <- 50
br_0 <- c(seq(min(rlog_0h_neutr), 0, length.out=ceiling(paletteLength/2) + 1), 
          seq(max(rlog_0h_neutr)/paletteLength, max(rlog_0h_neutr), length.out=floor(paletteLength/2)))

br <- c(seq(-2, 0, length.out=ceiling(paletteLength/2) + 1), 
        seq(max(rlog_0h_neutr)/paletteLength, 2, length.out=floor(paletteLength/2)))

myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)

# Plot heatmaps
pheatmap(norm_neutr, cluster_rows=FALSE, cluster_cols=FALSE, color = myColor, show_rownames = FALSE, main = "norm neutr")
pheatmap(vst_neutr, cluster_rows=FALSE, cluster_cols=FALSE, color = myColor, show_rownames = FALSE, main = "vst neutr")
pheatmap(rld_neutr, cluster_rows=FALSE, cluster_cols=FALSE, color = myColor, show_rownames = FALSE, main = "rlog neutr")

pheatmap(rlog_0h_neutr, cluster_rows=FALSE, cluster_cols=FALSE, breaks = br_0, color = myColor, show_rownames = TRUE, main = "rlog0h neutr")

LRT_neutr$row_colors <- ifelse(LRT_neutr$padj < 0.05, "sig", "no_sig")

rlog_0h_neutr_df <- data.frame(rlog_0h_neutr)
rlog_0h_neutr_df$status <- LRT_neutr$row_colors
rlog_0h_neutr_df$colors=ifelse(rlog_0h_neutr_df$status=="sig","red","black")

Neutr_h = pheatmap(rlog_0h_neutr, cluster_rows=FALSE, cluster_cols=FALSE, breaks = br, color = myColor, cellwidth = 15, cellheight = 15, gaps_col = c(4,8,11,15,19), angle_col = 90, show_rownames = TRUE, main = "Reactome DB neutrophil degranulation")

sig_cols=rlog_0h_neutr_df[order(match(rownames(rlog_0h_neutr_df), Neutr_h$gtable$grobs[[4]]$label)), ]$colors

Neutr_h$gtable$grobs[[4]]$gp=gpar(col=sig_cols)
Neutr_h

# pdf("plots/GSEA/Heatmaps/neutr_rlog0h", width=12, height=8)
# Neutr_h
# dev.off()