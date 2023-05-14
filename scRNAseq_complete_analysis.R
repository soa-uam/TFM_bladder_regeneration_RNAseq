rm(list = ls())
setwd("/home/sozaez/Desktop/TFM/scRNAseq/actualized_Olaya/scRNAseq_TFM_Olaya")

# harmony2 Seurat object - integration of 5 datasets:
#   - Yu. Homeostasis (normal conditions)
#   - TM (Tabula Muris). Homeostasis
#   - Baker. Homeostasis
#   - Cheng. 3 conditions: Cyclophosphamide(chronic & accute) + Homeostasis
#   - Sfakianos Nature Comm. 2021. Cancer

library(dplyr)
library(Seurat)
library(patchwork)
library(pheatmap)
library(ggplot2)
library(eulerr)
library(tidyverse)
library(clustree)
library(biomaRt)
library(viridis)

#### Data loading ####
# harmony2 <- readRDS("/home/sozaez/Desktop/TFM/scRNAseq/actualized_Olaya/scRNAseq_TFM_Olaya/harmony2.rds")
# harmony2

uroth <- readRDS("/home/sozaez/Desktop/TFM/scRNAseq/actualized_Olaya/scRNAseq_TFM_Olaya/uroth_Sandra.rds")
uroth # harmony subset w/o Sfakianos and clusters 9 and 11

uroth_Cheng <- readRDS("/home/sozaez/Desktop/TFM/scRNAseq/actualized_Olaya/scRNAseq_TFM_Olaya/uroth_Cheng_res0.3.rds")
uroth_res0.3 # uroth subset only Cheng urothelial cells w/o CP chronic

sub_condition_markers <- read.csv("results/1_DGE_Acute_vs_Homeostasis_Finmarkers/sub_condition_markers.csv")
colnames(sub_condition_markers)[1] <- "gene"

celltype_markers2 <- read.csv("results/2_DGE_Celltype_findallmarkers/celltype_markers2.csv")

Cheng_Uroth_Condition <- read.csv("results/3_DGE_Celltype_subCondition_findmarkers/markers_celltype_conditionCheng-Urothelial cells.csv")
colnames(Cheng_Uroth_Condition)[1] <- "gene"

Uroth_Condition <- read.csv("results/3_DGE_Celltype_subCondition_findmarkers/markers_celltype_conditionUrothelial cells.csv")
colnames(Uroth_Condition)[1] <- "gene"

Endothelial_Condition <- read.csv("results/3_DGE_Celltype_subCondition_findmarkers/markers_celltype_conditionEndothelial cells.csv")
colnames(Endothelial_Condition)[1] <- "gene"

Immune_Condition <- read.csv("results/3_DGE_Celltype_subCondition_findmarkers/markers_celltype_conditionImmune cells.csv")
colnames(Immune_Condition)[1] <- "gene"

SMCs_Condition <- read.csv("results/3_DGE_Celltype_subCondition_findmarkers/markers_celltype_conditionSMCs.csv")
colnames(SMCs_Condition)[1] <- "gene"

Stromal_Condition <- read.csv("results/3_DGE_Celltype_subCondition_findmarkers/markers_celltype_conditionStromal cells.csv")
colnames(Stromal_Condition)[1] <- "gene"

################## UMAP ###############################################
# DimPlot(harmony2, reduction = "umap", label = TRUE)
# DimPlot(harmony2, reduction = "umap", label = TRUE, group.by = "Condition")
# DimPlot(harmony2, reduction = "umap", label = TRUE, split.by = "Condition")
# DimPlot(harmony2, reduction = "umap", label = TRUE, group.by = "subID")
# DimPlot(harmony2, reduction = "umap", label = FALSE, group.by = "ID")

# pdf(file="plots/umap_uroth.pdf")
# DimPlot(uroth, reduction = "umap", label = TRUE)
# dev.off()

# Subset without Sfakianos, Baker-stromal and Baker-urothelial cells
# uroth <- subset(harmony2, subset = ( seurat_clusters == '9' | seurat_clusters == '11' |
#                                        ID == 'Sfakianos'), invert = TRUE)
# 
# saveRDS(uroth, file = "/home/sozaez/Desktop/TFM/scRNAseq/actualized_Olaya/scRNAseq_TFM_Olaya/uroth_Sandra.rds")

DimPlot(uroth, reduction = "umap", label = TRUE)
DimPlot(uroth, reduction = "umap", label = FALSE, group.by = "ID")
DimPlot(uroth, reduction = "umap", label = TRUE, group.by = "Condition")
DimPlot(uroth, reduction = "umap", label = TRUE, group.by = "CellType")
DimPlot(uroth, reduction = "umap", label = FALSE, group.by = "sub_Condition")
DimPlot(uroth, reduction = "umap", label = TRUE, group.by = "CellType", split.by = "sub_Condition")

DimPlot(uroth, reduction = "umap", label = TRUE, split.by = "Condition")
DimPlot(uroth, reduction = "umap", label = TRUE, split.by = "subID")
levels(uroth)

# pdf("plots/umap_uroth_bueno.pdf", width=12, height=8)
# DimPlot(uroth, reduction = "umap", label = TRUE)
# dev.off()

############################ DGEs ############################
# Access the data from the condition and subID
condition_data <- uroth@meta.data$Condition
condition_type <- uroth@meta.data$subID

# Create a new slot that combines the two sets of data based the previous conditions
uroth$sub_Condition <- ifelse(condition_data == "Homeostasis", condition_data,
                   ifelse(condition_type == "Cheng Cyclophosphamide-Acute", "Cyclophosphamide_Acute", 
                          ifelse(condition_type == "Cheng Cyclophosphamide-Chronic", "Cyclophosphamide_Chronic", NA)))

# saveRDS(uroth, file = "/home/sozaez/Desktop/TFM/scRNAseq/actualized_Olaya/scRNAseq_TFM_Olaya/uroth_Sandra.rds")

# DGE by sub_Condition: Cyclophosphamide_Acute vs. Homeostasis
Idents(uroth) <- "sub_Condition"

# Select all cells in clusters 
sub_condition_markers <- FindMarkers(uroth, ident.1 = "Cyclophosphamide_Acute", ident.2= "Homeostasis", min.pct=0.25, logfc.threshold=0.25)

sub_condition_markers <- sub_condition_markers %>%
  dplyr::arrange(p_val_adj)

write.csv(sub_condition_markers,
          file = "results/1_DGE_Acute_vs_Homeostasis_Finmarkers/sub_condition_markers.csv",
          quote = FALSE,
          row.names = TRUE)

write.csv(sub_condition_markers,
          file = "results/1_DGE_Acute_vs_Homeostasis_Finmarkers/sub_condition_markers.csv",
          quote = FALSE,
         row.names = TRUE)


# DGE by Condition: Cheng Cyclophosphamide-Acute (subID) vs. Homeostasis (subID) (by cluster)
Idents(uroth) <- seurat_clusters
uroth$cluster_subID <- paste(Idents(uroth), uroth$subID, sep = "_")
uroth$celltype <- Idents(uroth)
Idents(uroth) <- "cluster_subID"
output <- "DEG_CP_acute.csv"

clust_of_interest <- c("_Cheng Homeostais", "_Baker DMS Enriched",
                       "_Baker DMS Only", "_Baker Nuclei",
                       "_Baker Whole Bladder", "_TM", "_Yu")

for (i in 0:19){ # or however many clusters you have
  try({
    ident1 <- paste0(i,"_Cheng Cyclophosphamide-Acute")
    ident2 <- paste0(i, clust_of_interest)
    condition.diffgenes <- FindMarkers(uroth, ident.1 = ident1, ident.2= ident2, min.pct=0.25, logfc.threshold=0.25)
    write.csv(condition.diffgenes, file=paste0(output,i,".csv"))
  })
}

# view results
head(condition.diffgenes)

AcuteCP_vs_Homeo_markers <- condition.diffgenes %>%
  dplyr::arrange(p_val_adj)

# write.csv(AcuteCP_vs_Homeo_markers, 
#           file = "results/AcuteCP_vs_Homeo_markers.csv", 
#           quote = FALSE, 
#           row.names = TRUE)


# DEGs by Condition: Cyclophosphamide (Acute & Chronic) vs. Homeostasis
Idents(uroth) = "Condition"
overall_Condition_markers <- FindMarkers(uroth, ident.1 = "Cyclophosphamide", ident.2= "Homeostasis", min.pct=0.25, logfc.threshold=0.25)

overall_Condition_markers <- overall_Condition_markers %>%
  dplyr::arrange(p_val_adj)

write.csv(overall_Condition_markers,
          file = "results/overall_Condition_markers.csv",
          quote = FALSE,
          row.names = TRUE)


# DGE by CellType
Idents(uroth) <- "CellType"
celltype_markers <- FindAllMarkers(uroth, log2FC.threshold = 0.2,
                                   min.pct = 0.1, min.diff.pct = 0.2, only.pos = TRUE)

celltype_markers2 <- FindAllMarkers(uroth, log2FC.threshold = 0.25,
                                   min.pct = 0.1, only.pos = TRUE)

# write.csv(celltype_markers2, 
#           file = "results/2_DGE_Celltype_findallmarkers/celltype_markers2.csv", 
#           quote = FALSE, 
#           row.names = FALSE)


# DGE by celltype and Condition
# Access the data from the condition and damage_type slots
Celltype_data <- uroth@meta.data$CellType

# Create a new vector that combines the two sets of data based on your desired conditions and assign the new vector to a new slot
uroth$celltype_condition <- paste(Celltype_data, uroth$sub_Condition, sep = "_")
Idents(uroth) <- "celltype_condition"
output <- "markers_celltype_condition"

celltype_names <- unique(uroth$CellType)

for (i in celltype_names){
  try({
    ident1 <- paste0(i,"_Cyclophosphamide_Acute")
    ident2 <- paste0(i, "_Homeostasis")
    celltype_condition_markers <- FindMarkers(uroth, ident.1 = ident1, ident.2= ident2, min.pct=0.25, logfc.threshold=0.25)
    write.csv(celltype_condition_markers, file=paste0(output,i,".csv"))
  })
}

# view results
head(celltype_condition_markers)

celltype_condition_markers <- celltype_condition_markers %>%
  dplyr::arrange(p_val_adj)

write.csv(celltype_condition_markers, 
          file = "results/3_DGE_Celltype_subCondition_findmarkers/celltype_condition_markers.csv", 
          quote = FALSE, 
          row.names = TRUE)


uroth$celltype_uroth <- recode(uroth$CellType,
                               "Cheng-Urothelial cells" = "Urothelial cells",
                               "Urothelial cells" = "Urothelial cells")


# DGE by sub_Condition: Cyclophosphamide_Acute vs. Homeostasis
Idents(uroth) <- "celltype_uroth"



###################### Euler plot: DEGs CellTypes ######################
Urothelial_DEGs <- celltype_markers2[celltype_markers2$cluster == "Urothelial cells",]
Endothelial_DEGs <- celltype_markers2[celltype_markers2$cluster == "Endothelial cells",]
Baker_Uroth_DEGs <- celltype_markers2[celltype_markers2$cluster == "Baker-Urothelial cells",]
Baker_Stromal_DEGs <- celltype_markers2[celltype_markers2$cluster == "SMCs",]
Immune_cells_DEGs <- celltype_markers2[celltype_markers2$cluster == "Immune cells",]
Stromal_cells_DEGs <- celltype_markers2[celltype_markers2$cluster == "Stromal cells",]
Cheng_uroth_DEGs <- celltype_markers2[celltype_markers2$cluster == "Cheng-Urothelial cells",]
SMCs_DEGs <- celltype_markers2[celltype_markers2$cluster == "SMCs",]


# Euler plot: DEGs CellTypes by Acute vs Homeostasis condition
padj.cutoff <- 0.05

Urothelial_acute_up <- Uroth_Condition %>% dplyr::filter(p_val_adj < padj.cutoff & avg_log2FC >= 0) %>% pull(gene)
Urothelial_acute_down <- Uroth_Condition %>% dplyr::filter(p_val_adj < padj.cutoff & avg_log2FC < 0) %>% pull(gene)

Endothelial_acute_up <- Endothelial_Condition %>% dplyr::filter(p_val_adj < padj.cutoff & avg_log2FC >= 0) %>% pull(gene)
Endothelial_acute_down <- Endothelial_Condition %>% dplyr::filter(p_val_adj < padj.cutoff & avg_log2FC < 0) %>% pull(gene)

Immune_cells_acute_up <- Immune_Condition %>% dplyr::filter(p_val_adj < padj.cutoff & avg_log2FC >= 0) %>% pull(gene)
Immune_cells_acute_down <- Immune_Condition %>% dplyr::filter(p_val_adj < padj.cutoff & avg_log2FC < 0) %>% pull(gene)

Stromal_cells_acute_up <- Stromal_Condition %>% dplyr::filter(p_val_adj < padj.cutoff & avg_log2FC >= 0) %>% pull(gene)
Stromal_cells_acute_down <- Stromal_Condition %>% dplyr::filter(p_val_adj < padj.cutoff & avg_log2FC < 0) %>% pull(gene)

Cheng_uroth_acute_up <- Cheng_Uroth_Condition %>% dplyr::filter(p_val_adj < padj.cutoff & avg_log2FC >= 0) %>% pull(gene)
Cheng_uroth_acute_down <- Cheng_Uroth_Condition %>% dplyr::filter(p_val_adj < padj.cutoff & avg_log2FC < 0) %>% pull(gene)

SMCs_acute_up <- SMCs_Condition %>% dplyr::filter(p_val_adj < padj.cutoff & avg_log2FC >= 0) %>% pull(gene)
SMCs_acute_down <- SMCs_Condition %>% dplyr::filter(p_val_adj < padj.cutoff & avg_log2FC < 0) %>% pull(gene)


Uroth_damage_up <- intersect(Urothelial_DEGs$gene, Urothelial_acute_up)
Uroth_damage_down <- intersect(Urothelial_DEGs$gene, Urothelial_acute_down)
Cheng_uroth_damage_up <- intersect(Cheng_uroth_DEGs$gene, Cheng_uroth_acute_up)
Cheng_uroth_damage_down <- intersect(Cheng_uroth_DEGs$gene, Cheng_uroth_acute_down)
Endothelial_damage_up <- intersect(Endothelial_DEGs$gene, Endothelial_acute_up)
Endothelial_damage_down <- intersect(Endothelial_DEGs$gene, Endothelial_acute_down)
Immune_damage_up <- intersect(Immune_cells_DEGs$gene, Immune_cells_acute_up)
Immune_damage_down <- intersect(Immune_cells_DEGs$gene, Immune_cells_acute_down)
Stromal_damage_up <- intersect(Stromal_cells_DEGs$gene, Stromal_cells_acute_up)
Stromal_damage_down <- intersect(Stromal_cells_DEGs$gene, Stromal_cells_acute_down)
SMCs_damage_up <- intersect(SMCs_DEGs$gene, SMCs_acute_up)
SMCs_damage_down <- intersect(SMCs_DEGs$gene, SMCs_acute_down)

# Reassign name "Urothelial" to union "Urothelial" and "Cheng-uroth"
All_urothelial_up <- union(Uroth_damage_up, Cheng_uroth_damage_up)
All_urothelial_down <- union(Uroth_damage_down, Cheng_uroth_damage_down)

# Make list
list_DEGs_celltype_vs_condition_up <- list("Acute markers" = sub_condition_markers$gene, "Urothelial" = All_urothelial_up, "Endothelial" = Endothelial_damage_up, "Stromal" = Stromal_damage_up, "Immune" = Immune_damage_up, "SMCs" = SMCs_damage_up)

list_DEGs_celltype_vs_condition_down <- list("Acute markers" = sub_condition_markers$gene, "Urothelial" = All_urothelial_down, "Endothelial" = Endothelial_damage_down, "Stromal" = Stromal_damage_down, "Immune" = Immune_damage_down, "SMCs" = SMCs_damage_down)

plot(euler(list_DEGs_celltype_vs_condition_up), quantities = TRUE, legend = TRUE, adjust_labels = FALSE, main = "Acute vs Homeostasis Celltype DEGs - up")

plot(euler(list_DEGs_celltype_vs_condition_down), quantities = TRUE, legend = TRUE, adjust_labels = FALSE, main = "Acute vs Homeostasis Celltype DEGs - down")

# pdf("plots/acute_vs_homeostasis_celltype_DEGs_down.pdf", width=12, height=8)
# plot(euler(list_DEGs_celltype_vs_condition_down), quantities = TRUE, legend = TRUE, adjust_labels = FALSE, main = "Acute vs Homeostasis Celltype DEGs - down")
# dev.off()


###################### bulkRNAseq data vs scRNAseq data ######################
# Kmeans gene clustering data from LRT 0.01 (3839 genes)
kmeans_cl1 <- read.table("data/kmeans_cluster1.xls", header = TRUE)
kmeans_cl2 <- read.table("data/kmeans_cluster2.xls", header = TRUE)
kmeans_cl3 <- read.table("data/kmeans_cluster3.xls", header = TRUE)

kmeans_cl1 <- kmeans_cl1[,1]
kmeans_cl2 <- kmeans_cl2[,1]
kmeans_cl3 <- kmeans_cl3[,1]

OVERALL_up <- unique(c(sub_condition_markers$gene, All_urothelial_up, Endothelial_damage_up, Stromal_damage_up, Immune_damage_up, SMCs_damage_up))

OVERALL_down <- unique(c(All_urothelial_down, sub_condition_markers$gene, Endothelial_damage_down, Stromal_damage_down, Immune_damage_down, SMCs_damage_down))

all_kmeans3 <- read.table("/home/sozaez/Desktop/TFM/bulkRNAseq/Analyses/Uroth_reg_bulkRNAseq/data/clustered_genes3.xls", header = TRUE)
all_kmeans3 <- row.names(all_kmeans3)

# Comparison list
LRT_clusters <- list("Cluster1" = kmeans_cl1, 
                     "Cluster2" = kmeans_cl2, 
                     "Cluster3" = kmeans_cl3)


Uroth_celltypes_up <- list("Celltype1" = sub_condition_markers$gene, "Celltype2" = All_urothelial_up, "Celltype3" = Endothelial_damage_up, "Celltype4" = Stromal_damage_up, "Celltype5" = Immune_damage_up, "Celltype6" = SMCs_damage_up, "Celltype8" = OVERALL_up)

Uroth_celltypes_down <- list("Celltype1" = sub_condition_markers$gene, "Celltype2" = All_urothelial_down, "Celltype3" = Endothelial_damage_down, "Celltype4" = Stromal_damage_down, "Celltype5" = Immune_damage_down, "Celltype6" = SMCs_damage_down, "Celltype7" = OVERALL_down)


# First we create an empty matrix to store the results: 3x6 (clustersxcelltypes)
LRT_sc_matrix_up <- matrix(0, nrow = length(LRT_clusters), ncol = length(Uroth_celltypes_up))
LRT_sc_matrix_down <- matrix(0, nrow = length(LRT_clusters), ncol = length(Uroth_celltypes_down))

# Then we define the row and column names
rownames(LRT_sc_matrix_up) <- c("Cluster1", "Cluster2", "Cluster3")
colnames(LRT_sc_matrix_up) <- c("Celltype1", "Celltype2", "Celltype3", "Celltype4", "Celltype5", "Celltype 6", "Celltype 7")

rownames(LRT_sc_matrix_down) <- c("Cluster1", "Cluster2", "Cluster3")
colnames(LRT_sc_matrix_down) <- c("Celltype1", "Celltype2", "Celltype3", "Celltype4", "Celltype5", "Celltype6", "Celltype 7")


# Loop through the lists and calculate the intersection and store it in the result matrix
for(i in seq_along(LRT_clusters)){
  for(j in seq_along(Uroth_celltypes_up)){
    inter <- intersect(LRT_clusters[[i]], Uroth_celltypes_up[[j]])
    LRT_sc_matrix_up[i, j] <- length(inter) / length(OVERALL_up) * 100
  }
}

for(i in seq_along(LRT_clusters)){
  for(j in seq_along(Uroth_celltypes_down)){
    inter <- intersect(LRT_clusters[[i]], Uroth_celltypes_down[[j]])
    LRT_sc_matrix_down[i, j] <- length(inter) / length(OVERALL_down) * 100
  }
}

# Rename rows and columns
rownames(LRT_sc_matrix_up) <- c("Late responsive", "Transiently repressed", "Early responsive")
colnames(LRT_sc_matrix_up) <- c("Global acute", "Urothelial", "Endothelial", "Stromal", "Immune", "SMCs", "OVERALL")
LRT_sc_matrix_up <- LRT_sc_matrix_up[c(3,1,2),]


rownames(LRT_sc_matrix_down) <- c("Late responsive", "Transiently repressed", "Early responsive")
colnames(LRT_sc_matrix_down) <- c("Global acute", "Urothelial", "Endothelial", "Stromal", "Immune", "SMCs", "OVERALL")
LRT_sc_matrix_down <- LRT_sc_matrix_down[c(3,1,2),]


# write.csv(LRT_sc_matrix_up,
#           file = "results/LRT_sc_matrix_up.xls",
#           quote = FALSE,
#           row.names = TRUE)


###################### Wald ######################
Wald_up <- readRDS("data/list_all_Wald_DEup.RData")
Wald_down <- readRDS("data/list_all_Wald_DEdown.RData")
all_wald_DEGs <- read.table("data/All_Wald_DEGs.xls", header = TRUE)

Uroth_celltypes_up <- list("Celltype1" = sub_condition_markers$gene, "Celltype2" = All_urothelial_up, "Celltype3" = Endothelial_damage_up, "Celltype4" = Stromal_damage_up, "Celltype5" = Immune_damage_up, "Celltype6" = SMCs_damage_up, "Celltype7" = OVERALL_up)

Uroth_celltypes_down <- list("Celltype1" = sub_condition_markers$gene, "Celltype2" = All_urothelial_down, "Celltype3" = Endothelial_damage_down, "Celltype4" = Stromal_damage_down, "Celltype5" = Immune_damage_down, "Celltype6" = SMCs_damage_down, "Celltype7" = OVERALL_down)


# First we create an empty matrix to store the results: 3x6 (clustersxcelltypes)
Wald_scRNA_up <- matrix(0, nrow = length(Wald_up), ncol = length(Uroth_celltypes_up))
Wald_scRNA_down <- matrix(0, nrow = length(Wald_down), ncol = length(Uroth_celltypes_down))

# Loop through the lists and calculate the intersection and store it in the result matrix
for(i in seq_along(Wald_up)){
  for(j in seq_along(Uroth_celltypes_up)){
    inter <- intersect(Wald_up[[i]], Uroth_celltypes_up[[j]])
    Wald_scRNA_up[i, j] <- length(inter) / length(OVERALL_up) * 100
  }
}

for(i in seq_along(Wald_down)){
  for(j in seq_along(Uroth_celltypes_down)){
    inter <- intersect(Wald_down[[i]], Uroth_celltypes_down[[j]])
    Wald_scRNA_down[i, j] <- length(inter) / length(OVERALL_down) * 100
  }
}

# Then we define the row and column names
rownames(Wald_scRNA_up) <- c("12h", "24h", "72h", "7d")
colnames(Wald_scRNA_up) <- c("Global acute", "Urothelial", "Endothelial", "Stromal", "Immune", "SMCs", "OVERALL")
rownames(Wald_scRNA_down) <- c("12h", "24h", "72h", "7d")
colnames(Wald_scRNA_down) <- c("Global acute", "Urothelial", "Endothelial", "Stromal", "Immune", "SMCs", "OVERALL")

Wald_scRNA_up
Wald_scRNA_down


# write.csv(Wald_scRNA_up,
#           file = "results/Wald_scRNA_up.xls",
#           quote = FALSE,
#           row.names = TRUE)


################################# HEATMAPS #################################
# Define the maximum and minimum values for the color scale
max_val <- round(max(LRT_sc_matrix_up[,c(-1,-7)]))
min_val <- round(min(LRT_sc_matrix_up[,c(-1,-7)]))

# Create the pheatmap
pheatmap(LRT_sc_matrix_up,
         scale = "none",
         color = colorRampPalette(c("white", "blue"))(100),
         breaks = seq(min_val, max_val, length.out = 101),
         legend_breaks = seq(min_val, max_val, length.out = 5),
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         show_colnames = TRUE,
         show_rownames = TRUE,
         display_numbers = round(LRT_sc_matrix_up, 2),
         number_color = "grey",
         main = "bulkRNAseq LRT vs scRNAseq - Acute",
         border_color = "gray50"
)


# Define the maximum and minimum values for the color scale
max_val <- round(max(LRT_sc_matrix_down[,c(-1,-7)]))
min_val <- round(min(LRT_sc_matrix_down[,c(-1,-7)]))

# Create the pheatmap
pheatmap(LRT_sc_matrix_down,
         scale = "none",
         color = colorRampPalette(c("white", "red"))(100),
         breaks = seq(min_val, max_val, length.out = 101),
         legend_breaks = seq(min_val, max_val, length.out = 5),
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         show_colnames = TRUE,
         show_rownames = TRUE,
         display_numbers = round(LRT_sc_matrix_down, 2),
         number_color = "black",
         main = "bulkRNAseq LRT vs scRNAseq - Acute",
         border_color = "gray50"
)

# Define the maximum and minimum values for the color scale
max_val <- round(max(Wald_scRNA_up[,c(-1,-7)]))
min_val <- round(min(Wald_scRNA_up[,c(-1,-7)]))

# Create the pheatmap
pheatmap(Wald_scRNA_up,
         scale = "none",
         color = colorRampPalette(c("white", "blue"))(100),
         breaks = seq(min_val, max_val, length.out = 101),
         legend_breaks = seq(min_val, max_val, length.out = 5),
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         show_colnames = TRUE,
         show_rownames = TRUE,
         display_numbers = round(Wald_scRNA_up, 2),
         number_color = "black",
         main = "bulkRNAseq Wald vs scRNAseq - Acute",
         border_color = "gray50"
)



# Define the maximum and minimum values for the color scale
max_val <- round(max(Wald_scRNA_down[,c(-1,-7)]))
min_val <- round(min(Wald_scRNA_down[,c(-1,-7)]))

# Create the pheatmap
pheatmap(Wald_scRNA_down,
         scale = "none",
         color = colorRampPalette(c("white", "red"))(100),
         breaks = seq(min_val, max_val, length.out = 101),
         legend_breaks = seq(min_val, max_val, length.out = 5),
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         show_colnames = TRUE,
         show_rownames = TRUE,
         display_numbers = round(Wald_scRNA_down, 2),
         number_color = "black",
         main = "bulkRNAseq Wald vs scRNAseq - Acute",
         border_color = "gray50"
)


################################# DOTPLOTS #################################
LRT_sc_matrix_up.df <- as.data.frame(LRT_sc_matrix_up)

LRT_sc_matrix_up.df %>%
  rownames_to_column(var = "pattern") %>%
  gather(celltype, overlap, -pattern) %>%
  mutate(celltype = factor(celltype, levels = c("Immune","SMCs","Stromal","Endothelial", "Urothelial", "Global acute", "OVERALL"))) %>%
  ggplot() +
  geom_point(aes(x = pattern, y = celltype, size = overlap, color=overlap)) +
  scale_size_continuous(breaks=0:5)+
  scale_color_gradient(low="lightblue", high="darkblue")+
  labs(title = "LRT bulk vs. scRNAseq - up")+
  theme(axis.text = element_text(size = 10))+
  coord_fixed(ratio = 0.5) # Set the aspect ratio to 1:1 # Adjust plot margins


LRT_sc_matrix_down.df <- as.data.frame(LRT_sc_matrix_down)
LRT_sc_matrix_down.df %>%
  rownames_to_column(var = "pattern") %>%
  gather(celltype, overlap, -pattern) %>%
  mutate(celltype = factor(celltype, levels = c("Immune","SMCs","Stromal","Endothelial", "Urothelial", "Global acute", "OVERALL"))) %>%
  ggplot() +
  geom_point(aes(x = pattern, y = celltype, size = overlap, color=overlap)) +
  scale_size_continuous(breaks=0:5) +
  scale_color_gradient(low="lightblue", high="darkblue")+
  labs(title = "LRT bulk vs. scRNAseq - down")+
  theme(axis.text = element_text(size = 10))+
  coord_fixed(ratio = 0.5)


Wald_scRNA_up.df <- as.data.frame(Wald_scRNA_up)
Wald_scRNA_up.df %>%
  rownames_to_column(var = "pattern") %>%
  gather(celltype, overlap, -pattern) %>%
  mutate(celltype = factor(celltype, levels = c("Immune","SMCs","Stromal","Endothelial", "Urothelial", "Global acute", "OVERALL"))) %>%
  ggplot() +
  geom_point(aes(x = pattern, y = celltype, size = overlap, color=overlap)) +
  scale_size_continuous(breaks=0:5) +
  scale_color_gradient(low="lightblue", high="darkblue")+
  labs(title = "Wald bulk vs. scRNAseq - up")+
  theme(axis.text = element_text(size = 10))+
  coord_fixed(ratio = 0.5)


Wald_scRNA_down.df <- as.data.frame(Wald_scRNA_down)
Wald_scRNA_down.df %>%
  rownames_to_column(var = "pattern") %>%
  gather(celltype, overlap, -pattern) %>%
  mutate(celltype = factor(celltype, levels = c("Immune","SMCs","Stromal","Endothelial", "Urothelial", "Global acute", "OVERALL"))) %>%
  ggplot() +
  geom_point(aes(x = pattern, y = celltype, size = overlap, color=overlap)) +
  scale_size_continuous(breaks=0:5) +
  scale_color_gradient(low="lightblue", high="darkblue")+
  labs(title = "Wald bulk vs. scRNAseq - down")+
  theme(axis.text = element_text(size = 10))+
  coord_fixed(ratio = 0.5)


######################## FEATURE PLOTS ################################
library("viridis")

DEGs_damage_up <- lapply(list_DEGs_celltype_vs_condition_up, function(x) x[1:5])
DEGs_damage_down <- lapply(list_DEGs_celltype_vs_condition_down, function(x) x[1:5])

pal <- viridis(n=10, option='D')
umap_theme <- theme(
  axis.line=element_blank(),
  axis.text.x=element_blank(),
  axis.text.y=element_blank(),
  axis.ticks=element_blank(),
  axis.title.x=element_blank(),
  axis.title.y=element_blank(),
  panel.background=element_blank(),
  panel.border=element_blank(),
  panel.grid.major=element_blank(),
  panel.grid.minor=element_blank(),
  plot.background=element_blank()
)

# Create feature plots
plot_list <- FeaturePlot(
  uroth,
  features=unlist(DEGs_damage_up),
  combine=FALSE, 
  cols=pal,
  max.cutoff='q98'
)
for(i in 1:length(plot_list)){
  plot_list[[i]] <- plot_list[[i]] + umap_theme 
}
png("/home/sozaez/Desktop/TFM/Figuras/damage_up_featurePlot.png", width=35, height=20, units='in', res=200)
CombinePlots(plot_list, ncol=6)
dev.off()

# Create feature plots
plot_list <- FeaturePlot(
  uroth,
  features=unlist(DEGs_damage_down),
  combine=FALSE, 
  cols=pal,
  max.cutoff='q98'
)
for(i in 1:length(plot_list)){
  plot_list[[i]] <- plot_list[[i]] + umap_theme 
}
png("/home/sozaez/Desktop/TFM/Figuras/damage_down_featurePlot.png", width=35, height=20, units='in', res=200)
CombinePlots(plot_list, ncol=6)
dev.off()

FeaturePlot(uroth, features = c("Mki67", "Krt14", "Trp63", "Krt5", "Krt8", "Krt15", "Krt6a", "Psca", "Upk3a", "Krt20"))

# Damage & celltype specific markers: Uroth Cheng, Uroth, Endothelial, SMCs, Stromal, Immune

# Cheng_uroth damage up
FeaturePlot(uroth, features = c("Cenpm", "Tubb4b", "Lsm2", "Ccng1", "1810037I17Rik"))
VlnPlot(uroth, features = c("Cenpm", "Tubb4b", "Lsm2", "Ccng1", "1810037I17Rik"))
# Cheng_uroth damage down
FeaturePlot(uroth, features = c("Stmn1", "Hmgb2", "Lockd", "Top2a", "Cenpf"))
# Uroth damage up
FeaturePlot(uroth, features = c("Sptssb", "2200002D01Rik", "Sprr1a",  "Krt18", "Upk3a"))
# Uroth damage down
FeaturePlot(uroth, features = c("Ly6d", "Gsta4", "Krt15", "Gstm1", "Sfn"  ))
# Endothelial damage up
FeaturePlot(uroth, features = c("Tm4sf1",  "Ctla2a",  "Tmem252", "Cd200", "Apold1"))
# Endothelial damage down
FeaturePlot(uroth, features = c("Cldn5",  "Egfl7",  "Esam", "Pecam1", "Plvap"))
# SMCs damage up
FeaturePlot(uroth, features = c("Rasl11a", "Rrad", "Gm13889", "Mustn1",  "Sncg"))
# SMCs damage down
FeaturePlot(uroth, features = c("Myh11", "Tagln", "Acta2", "Actg2", "Mylk" ))
# Stromal up
FeaturePlot(uroth, features = c("Lgals1", "Ifi27l2a", "Rbp1", "Serpina3n", "Bglap2"))
# Stromal down
FeaturePlot(uroth, features = c("Dcn", "Col1a2", "Col1a1", "Mgp",    "Sparc"))
# Immune up
FeaturePlot(uroth, features = c("Ccl4", "H2-Aa",  "H2-Ab1", "C1qb",   "C1qa"))
# Immune down
FeaturePlot(uroth, features = c("Cd74", "H2.Ab1", "H2.Aa",  "H2.Eb1", "Apoe"))

           
# DAMAGE FEATURE PLOTS UP
FeaturePlot(uroth, reduction = "umap", features = c("Cenpm", "Sprr1a","Krt18", "Upk3a","Tm4sf1", "Rrad", "Rbp1", "Ccl4","H2-Aa"))
FeaturePlot(uroth, reduction = "umap", features = c("Cenpm", "Sprr1a","Krt18", "Upk3a"), split.by = "sub_Condition")
VlnPlot(uroth, features = c("Cenpm", "Sptssb", "2200002D01Rik"), split.by = "sub_Condition")

# DAMAGE FEATURE PLOTS DOWN
FeaturePlot(uroth, reduction = "umap", features = c("Hmgb2", "Ly6d", "Krt15", "Gsta4","Gstm1", "Cldn5", "Acta2", "Mylk", "Col1a1", "H2.Ab1", "Apoe"))
VlnPlot(uroth, features = c("Dcn", "Col1a2", "Col1a1", "Mgp",    "Sparc"), split.by = "sub_Condition")

# With the markers genes: Mki67, Krt14, Palp, Upk3a, Krt5, Krt20
VlnPlot(uroth, features = c("Mki67", "Krt14", "Upk3a", "Krt5", "Krt20"))



########################## Subset urothelial cells ######################
uroth_onlyUroth <- subset(uro_Sandra, subset = (celltype_uroth == "urothelial cells"))

saveRDS(uroth_onlyUroth, file = "./uroth_onlyUroth.rds")

DimPlot(uroth_onlyUroth, reduction = "umap", label = TRUE, split.by = "sub_Condition")

only_uroth <- readRDS("/home/sozaez/Desktop/TFM/scRNAseq/actualized_Olaya/scRNAseq_TFM_Olaya/deconvolution/uroth_onlyUroth_edit.rds")

DimPlot(only_uroth, reduction = "umap", label = TRUE, split.by = "sub_Condition")

################## CELL CLUSTERING
only_uroth <- FindNeighbors(only_uroth, dims = 1:20)
only_uroth <- FindClusters(only_uroth, resolution = c(0,0.001, 0.002, 0.003, 0.004, 0.005, 0.01))

clustree(only_uroth, prefix = "SCT_snn_res.")

################## UMAP
only_uroth <- only_uroth %>%
  RunUMAP(reduction = "pca", dims = 1:20) %>% 
  FindNeighbors(reduction = "pca", dims = 1:20) %>% 
  FindClusters(resolution = 0.01)  %>%
  identity() # 5 clusters

DimPlot(only_uroth, reduction = "umap", label = TRUE)
DimPlot(only_uroth, reduction = "umap", label = TRUE, split.by = "sub_Condition")


########################## Subset Cheng cells ######################
uroth$celltype_uroth <- recode(uroth$CellType,
                               "Cheng-Urothelial cells" = "Urothelial cells",
                               "Urothelial cells" = "Urothelial cells")


# DGE by sub_Condition: Cyclophosphamide_Acute vs. Homeostasis
Idents(uroth) <- "celltype_uroth"

uroth_Cheng <- subset(uroth, subset = (ID == "Cheng" & celltype_uroth == "Urothelial cells"))
uroth_Cheng <- subset(uroth_Cheng, subset = (sub_Condition == "Cyclophosphamide_Chronic"), invert = TRUE)

saveRDS(uroth_Cheng, file = "./uroth_Cheng.rds")

DimPlot(uroth_Cheng, reduction = "umap", label = TRUE)

Idents(uroth_Cheng) <- "seurat_clusters"


################## CELL CLUSTERING
uroth_Cheng <- FindNeighbors(uroth_Cheng, dims = 1:20)
uroth_Cheng <- FindClusters(uroth_Cheng, resolution = c(0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1))

clustree(uroth_Cheng, prefix = "RNA_snn_res.")

# pdf("/home/sozaez/Desktop/clustree.pdf",width=8,height=10)
# clustree(uroth_res0.3, prefix = "RNA_snn_res.")
# dev.off()

################## UMAP
uroth_Cheng <- uroth_Cheng %>%
  RunUMAP(reduction = "pca", dims = 1:20) %>% 
  FindNeighbors(reduction = "pca", dims = 1:20) %>% 
  FindClusters(resolution = 0.03)  %>%
  identity() # 5 clusters

uroth_Cheng <- subset(uroth_Cheng, subset = ( seurat_clusters == '3' | seurat_clusters == '4' | seurat_clusters == '5' ), invert = TRUE)

saveRDS(uroth_Cheng, file = "./3clust_uroth_Cheng.rds")
saveRDS(uroth_Cheng, file = "./3clust_uroth_Cheng.rds")

DimPlot(uroth_Cheng, reduction = "umap", label = TRUE)
DimPlot(uroth_Cheng, reduction = "umap", label = TRUE, split.by = "sub_Condition")

FeaturePlot(uroth_Cheng, features = c("Mki67", "Krt14", "Trp63", "Krt5", "Krt8", "Krt6a", "Psca", "Upk3a", "Krt20"))

uroth_markers <- FindAllMarkers(uroth_Cheng, logfc.threshold = 0.5, min.pct = 0.1, only.pos = TRUE)

# write.csv(uroth_markers,'All_markers_Cheng3clusters.csv')

# Define a list of marker genes for each cell type
uroth_markers <- list("Basal" = c("Mki67", "Krt14", "Trp63", "Cenpa", "Itga6", "Itgb4", "Col17a1", "Col18a1"), 
  "Intermediate" = c("Krt4", "Krt5", "Krt6a", "Krt6b", "Cdh3", "Cdh1", "S100a14", "S100a16"), "Luminal" = c("Psca", "Gata3", "Foxa1", "Krt20", "Upk1a", "Upk1b", "Upk2", "Upk3a"))

# Calculate module scores for each cell type
uroth_Cheng <- AddModuleScore(object = uroth_Cheng, features = uroth_markers, prefix = "uroth_Sub_celltype", assay = "RNA")

# Find markers that are differentially expressed between cell types
uroth_Cheng <- FindAllMarkers(object = uroth_Cheng, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA")

FeaturePlot(uroth_Cheng, features = c("Mki67", "Krt14", "Trp63", "Cenpa", "Itgb4", "Col17a1", "Col18a1"))
FeaturePlot(uroth_Cheng, features = c("Krt4", "Krt5", "Krt6a", "Krt6b", "Cdh3", "Cdh1", "Trp63", "S100a14", "S100a16"))
FeaturePlot(uroth_Cheng, features = c("Psca", "Gata3", "Foxa1", "Krt20", "Upk1a", "Upk1b", "Upk2", "Upk3a"))

uroth_Cheng <- readRDS("/home/sozaez/Desktop/TFM/scRNAseq/actualized_Olaya/scRNAseq_TFM_Olaya/3clust_uroth_Cheng.rds")


#### Annotate cell clusters ####
cluster_annotations <- list(
  '0' = 'C0-Intermediate',
  '1' = 'C1-Intermediate',
  '2' = 'C2-Resting basal',
  '3' = 'C3-Luminal',
  '4' = 'C4-Luminal acute',
  '5' = 'C5-Basal proliferative',
  '6' = 'C6-Intermediate high',
  '7' = 'C7-Intermediate low'
)

celltype_annotations <- list(
  '0' = 'Intermediate',
  '1' = 'Intermediate',
  '2' = 'Resting basal',
  '3' = 'Luminal',
  '4' = 'Luminal acute',
  '5' = 'Basal proliferative',
  '6' = 'Intermediate high',
  '7' = 'Intermediate low'
)


# Add CellType to metadata
uroth_Cheng$CellType_cluster <- unlist(cluster_annotations[uroth_Cheng$seurat_clusters])
uroth_Cheng$CellType <- unlist(celltype_annotations[uroth_Cheng$seurat_clusters])


uroth_Cheng_res0.3 <- subset(uroth_Cheng, subset = ( seurat_clusters == '6' | seurat_clusters == '7' ), invert = TRUE)

saveRDS(uroth_Cheng_res0.3, file = "./uroth_Cheng_res0.3.rds")



# Visualize canonical marker genes ####

list_DEGs_celltype_vs_condition_up <- list_DEGs_celltype_vs_condition_up[c(2,3,4,5,6)]
list_DEGs_celltype_vs_condition_down <- list_DEGs_celltype_vs_condition_down[c(2,3,4,5,6)]

list_damage_up <- lapply(list_DEGs_celltype_vs_condition_up, head, n = 5)
list_damage_down <- lapply(list_DEGs_celltype_vs_condition_down, head, n = 5)

pal <- viridis(n=10, option="plasma")

umap_theme <- theme(
  axis.line=element_blank(),
  axis.text.x=element_blank(),
  axis.text.y=element_blank(),
  axis.ticks=element_blank(),
  axis.title.x=element_blank(),
  axis.title.y=element_blank(),
  panel.background=element_blank(),
  panel.border=element_blank(),
  panel.grid.major=element_blank(),
  panel.grid.minor=element_blank(),
  plot.background=element_blank()
)


# Create feature plots up
plot_list_up <- FeaturePlot(
  uroth,
  features=unlist(list_damage_up),
  combine=FALSE, 
  cols=pal,
  max.cutoff='q98'
)

for(i in 1:length(plot_list_up)){
  plot_list_up[[i]] <- plot_list_up[[i]] + umap_theme 
}

png('damage_up.png', width=35, height=20, units='in', res=200)
CombinePlots(plot_list_up, ncol=6)
dev.off()


# Create feature plots
plot_list_down <- FeaturePlot(
  uroth,
  features=unlist(list_damage_down),
  combine=FALSE, 
  cols=pal,
  max.cutoff='q98'
)

for(i in 1:length(plot_list_down)){
  plot_list_down[[i]] <- plot_list_down[[i]] + umap_theme 
}
png('Featplot_damage_down.png', width=35, height=20, units='in', res=200)
CombinePlots(plot_list_down, ncol=6)
dev.off()


all_damage <- c("Sptssb", "2200002D01Rik", "Ctla2a", "Rbp1", "Serpina3n", "Ccl4", "Gm13889", "krt15", "Sfn", "Egfl7", "Col1a1", "Cd74", "Myh11", "Mylk")


# Create feature plots
plot_list <- FeaturePlot(
  uroth,
  features=all_damage,
  combine=FALSE, 
  cols=pal,
  max.cutoff='q98'
)

for(i in 1:length(plot_list)){
  plot_list[[i]] <- plot_list[[i]] + umap_theme 
}

png('all_damage_up.png', width=35, height=20, units='in', res=200)
CombinePlots(plot_list, ncol=6)
dev.off()

all_damage_down <- c("Ly6d", "Sfn", "Cldn5", "Dcn", "Col1a1", "Cd74", "H2.Ab1", "Myh1", "Mylk", "Egfl7", "Col1a1", "Cd74", "Myh11", "Mylk")

# Urothelial - Ly6d, Sfn
# Endothelial cells - Cldn5
# Stromal cells - Dcn, Col1a1
# Immune - Cd74, H2.Ab1
# SMCs - Myh1, Mylk

# Create feature plots
plot_list <- FeaturePlot(
  uroth,
  features=all_damage,
  combine=FALSE, 
  cols=pal,
  max.cutoff='q98'
)

for(i in 1:length(plot_list)){
  plot_list[[i]] <- plot_list[[i]] + umap_theme 
}

png('all_damage_up.png', width=35, height=20, units='in', res=200)
CombinePlots(plot_list, ncol=6)
dev.off()

combined_damage_feature <- CombinePlots(plot_list, ncol=6)

png('combined_damage.png', width=35, height=20, units='in', res=200)
CombinePlots(plots = list(combined_damage_feature, p3))
dev.off()

# Violin plots
# plot distributions for marker genes:
p <- VlnPlot(
  uroth,
  group.by = 'celltype_uroth',
  split.by='sub_Condition',
  features=unlist(list_damage_up),
  pt.size = 0, ncol=5
)

png(('vlnPlot_damage_up.png'), width=30, height=1*length(unlist(list_damage_up)), units='in', res=200)
print(p)
dev.off()


# plot distributions for marker genes:
p2 <- VlnPlot(
  uroth,
  group.by = 'celltype_uroth',
  split.by='sub_Condition',
  features=unlist(list_damage_down),
  pt.size = 0, ncol=5
)

png(('vlnPlot_damage_down.png'), width=30, height=1*length(unlist(list_damage_up)), units='in', res=200)
print(p2)
dev.off()


# Violin plots
# plot distributions for marker genes:
p3 <- VlnPlot(
  uroth,
  group.by = 'celltype_uroth',
  split.by='sub_Condition',
  features=all_damage,
  pt.size = 0, ncol=5
)

png(('Vlnplot_all_damage.png'), width=30, height=1*length(all_damage), units='in', res=200)
print(p3)
dev.off()



# # Cell Sequence and Cell Label
# write.table(uroth_Cheng@active.ident, file='deconvolution/sub_Convert_Barcode_Label_8K.tsv', quote=FALSE, sep='\t', col.names = TRUE)

# # Gene counts per cell
# write.table(uroth_Cheng@assays[["RNA"]]@counts, file='deconvolution/sub_Gene_Count_per_Cell_8K.tsv', quote=FALSE, sep='\t', col.names = TRUE)



############################## CybersortX input #############################
uroth_ensmbl <- readRDS("/home/sozaez/Desktop/TFM/scRNAseq/actualized_Olaya/scRNAseq_TFM_Olaya/uroth_Ensembl.rds")

uroth$celltype_uroth <- recode(uroth$CellType,
                               "Cheng-Urothelial cells" = "Urothelial cells",
                               "Urothelial cells" = "Urothelial cells")


# DGE by sub_Condition: Cyclophosphamide_Acute vs. Homeostasis
Idents(uroth) <- "celltype_uroth"

uroth_subset <- subset(uroth[, sample(colnames(uroth), size = 8000, replace = F)])

saveRDS(uroth_subset, file = "/home/sozaez/Desktop/TFM/scRNAseq/actualized_Olaya/scRNAseq_TFM_Olaya/deconvolution/uroth_subset8K.rds")


############################ ENSEMBL ID ###############################
# We use the following code to eliminate "." from EnsemblIDs

# Get the gene symbols from the Seurat object
gene_symbols <- rownames(uroth_subset@assays$RNA@data)
uroth_genes <- as.data.frame(gene_symbols)

# extract the row names of the counts matrix
gene_names <- rownames(uroth_subset@assays[["RNA"]]@counts)

# extract genes that start with ENSMUSG, have numbers in the middle, and end with ".numbers"
cut_genes <- gsub("(ENSMUSG\\d+)\\.\\d+", "\\1", uroth_genes$gene_symbols)
gene_names <- gsub("(ENSMUSG\\d+)\\.\\d+", "\\1", gene_names)

# create new column in original data frame with cut gene names
uroth_genes$cut_genes <- cut_genes

# update the `rownames` slot of the `rowData` object with the modified gene names as well
rowData(uroth_subset)$ensembl_gene_id <- gene_names

# Compare with bulkRNAseq genes
data <- read.table("/home/sozaez/Desktop/TFM/bulkRNAseq/Analyses/Uroth_reg_bulkRNAseq/data/data_uro.txt", header = TRUE, row.names = 1)

length(which(uroth_genes$gene_symbols %in% row.names(data)))
length(which(uroth_genes$cut_genes %in% row.names(data)))

length(setdiff(uroth_genes$gene_symbols, row.names(data)))
length(setdiff(uroth_genes$cut_genes, row.names(data)))

length(setdiff(uroth_genes$gene_symbols, row.names(data)))
length(setdiff(uroth_genes$cut_genes, row.names(data)))

# Replace the old gene names with the new gene names in the Seurat object
rownames(uroth_subset@assays$RNA@data) <- uroth_genes$cut_genes
# set the modified gene names as the new row names of the counts matrix
rownames(uroth_subset@assays[["RNA"]]@counts) <- gene_names


saveRDS(uroth, file = "/home/sozaez/Desktop/TFM/scRNAseq/actualized_Olaya/scRNAseq_TFM_Olaya/uroth_Ensembl.rds")

# saveRDS(uroth_subset, file = "/home/sozaez/Desktop/TFM/scRNAseq/actualized_Olaya/scRNAseq_TFM_Olaya/deconvolution/uroth_subset8K.rds")

uroth_subset <- readRDS("/home/sozaez/Desktop/TFM/scRNAseq/actualized_Olaya/scRNAseq_TFM_Olaya/deconvolution/uroth_subset8K.rds")

# Cell Sequence and Cell Label
write.table(uroth_subset@active.ident, file='deconvolution/sub_Convert_Barcode_Label_8K.tsv', quote=FALSE, sep='\t', col.names = TRUE)

# Gene counts per cell
write.table(uroth_subset@assays[["RNA"]]@counts, file='deconvolution/sub_Gene_Count_per_Cell_8K.tsv', quote=FALSE, sep='\t', col.names = TRUE)

# Access the data in the RNA assay with the 2000 variable features
rna_vars <- uroth_subset@assays$RNA@counts[match(uroth_subset@assays$RNA@var.features, row.names(uroth_subset@assays$RNA@counts)),]

# Write the data to a TSV file
write.table(rna_vars, file = 'deconvolution/sub_Gene_Count_per_Cell_8K_2Kvars.tsv', quote = FALSE, sep = '\t', col.names = TRUE)



#################### CIBERSORT X RESTULS #######################
library(tidyr)
library(ggplot2)

# Read in the CSV file
cibersortx_output <- read.csv("/home/sozaez/Downloads/CIBERSORTx_Job28_Results.csv", header = TRUE)

cibersortx_output_clean <- cibersortx_output[c(-7,-8,-9)]
colnames(cibersortx_output_clean) <- c("Mixture", "Intermediate", "Resting basal", "Luminal acute", "Basal proliferative", "Luminal")

cibersortx_output_clean <- select(cibersortx_output_clean, c("Mixture", "Basal proliferative", "Resting basal", "Intermediate", "Luminal", "Luminal acute"))


# Reshape the data into long format
cibersortx_output_long <- reshape2::melt(cibersortx_output_clean, variable.name = "cell_type", value.name = "proportion", na.rm = TRUE)

angle <- 45
# Create a stacked bar chart
ggplot(cibersortx_output_long, aes(x = Mixture, y = proportion, fill = cell_type)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(labels = c("Basal proliferative", "Resting basal", "Intermediate", "Luminal", "Luminal acute"), values = c("#1F78B4", "#33A02C", "#E31A1C", "#FF7F00", "#6A3D9A")) +
  theme(axis.text.x = element_text(angle = angle, hjust = 1))


# Add the RMSE column to the long format data frame
cibersortx_output_long$rmse <- cibersortx_output$RMSE

# geom_errorbar(aes(ymin = proportion - rmse, ymax = proportion + rmse), width = 0.2)+
# scale_y_continuous(labels = scales::percent_format()) +
# scale_size(range = c(1, 3)) # Resize the error bars based on RMSE




################# HEATMAP #####################33
# Cargar la tabla de resultados cibersort.csv
cibersortx <- read.csv("/home/sozaez/Downloads/CIBERSORTx_Job24_Results.csv", header = TRUE, row.names = 1)

cibersortx_clean <- cibersortx[c(-6, -7,-8)]
colnames(cibersortx_clean) <- c("Urothelial cells", "Stromal cells", "SMCs", "Endothelial cells", "Immune cells")

cibersortx_clean <- select(cibersortx_clean, c("Urothelial cells", "Endothelial cells", "Stromal cells", "SMCs", "Immune cells"))


# Set pheatmap color parameters
paletteLength <- 20
# myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
myColor <- colorRampPalette(c("white", "lightblue", "blue"))(paletteLength)

myColor <- colorRampPalette(brewer.pal(9,"Blues"))(23)
# length(breaks) == length(paletteLength) + 1
# use floor and ceiling to deal with even/odd length pallettelengths
myBreaks <- c(seq(min(cibersortx_clean), 0.05, length.out=ceiling(paletteLength/2) + 1), 0.79,
              seq(max(cibersortx_clean)/paletteLength, max(cibersortx_clean), length.out=floor(paletteLength/2)))


# Create the heatmap with the updated breaks and colors
pheatmap::pheatmap(cibersortx_clean, 
                   cluster_rows = FALSE, 
                   cluster_cols = FALSE, 
                   display_numbers = round(cibersortx_clean, 3), 
                   fontsize_row = 8,
                   fontsize_col = 8, 
                   angle_col = 45,
                   show_rownames = TRUE, 
                   show_colnames = TRUE, 
                   breaks = myBreaks,
                   color=myColor,
                   number_color = ifelse(cibersortx_clean >= 0.6, "white", "black"))


# Add the numbers to the heatmap with different colors
# geom_text(aes(label = round(value, 2), color = ifelse(value >= mean(value), "black", "white"))) +



########################### Tabla 2
# create a new data frame to hold the results
upreg_damage <- data.frame(matrix(nrow = 50, ncol = 6))

# loop over each vector in the list and extract the first 50 occurrences
for (i in 1:6) {
  upreg_damage[, i] <- list_DEGs_celltype_vs_condition_up[[i]][1:50]
}


# set the column names of the new table to match the vector names in the list
colnames(upreg_damage) <- names(list_DEGs_celltype_vs_condition_up)

# view the resulting table
View(upreg_damage)



# create a new data frame to hold the results
downreg_damage <- data.frame(matrix(nrow = 50, ncol = 6))

# loop over each vector in the list and extract the first 50 occurrences
for (i in 1:6) {
  downreg_damage[, i] <- list_DEGs_celltype_vs_condition_down[[i]][1:50]
}

# set the column names of the new table to match the vector names in the list
colnames(downreg_damage) <- names(list_DEGs_celltype_vs_condition_down)

# view the resulting table
View(downreg_damage)


# create a logical vector indicating which genes belong to each string group
uroth_cl1_up <- upreg_damage$Urothelial %in% unlist(strsplit(kmeans_cl1, " "))
uroth_cl2_up <- upreg_damage$Urothelial %in% unlist(strsplit(kmeans_cl2, " "))
uroth_cl3_up <- upreg_damage$Urothelial %in% unlist(strsplit(kmeans_cl3, " "))

endo_cl1_up <- upreg_damage$Endothelial %in% unlist(strsplit(kmeans_cl1, " "))
endo_cl2_up <- upreg_damage$Endothelial %in% unlist(strsplit(kmeans_cl2, " "))
endo_cl3_up <- upreg_damage$Endothelial %in% unlist(strsplit(kmeans_cl3, " "))

strom_cl1_up <- upreg_damage$Stromal %in% unlist(strsplit(kmeans_cl1, " "))
strom_cl2_up <- upreg_damage$Stromal %in% unlist(strsplit(kmeans_cl2, " "))
strom_cl3_up <- upreg_damage$Stromal %in% unlist(strsplit(kmeans_cl3, " "))

immune_cl1_up <- upreg_damage$Immune %in% unlist(strsplit(kmeans_cl1, " "))
immune_cl2_up <- upreg_damage$Immune %in% unlist(strsplit(kmeans_cl2, " "))
immune_cl3_up <- upreg_damage$Immune %in% unlist(strsplit(kmeans_cl3, " "))

SMCs_cl1_up <- upreg_damage$SMCs %in% unlist(strsplit(kmeans_cl1, " "))
SMCs_cl2_up <- upreg_damage$SMCs %in% unlist(strsplit(kmeans_cl2, " "))
SMCs_cl3_up <- upreg_damage$SMCs %in% unlist(strsplit(kmeans_cl3, " "))

# Print the genes that belong to each string group
cat("Uroth in string 1:", upreg_damage$Urothelial[uroth_cl1_up], "\n")
cat("Uroth in string 2:", upreg_damage$Urothelial[uroth_cl2_up], "\n")
cat("Uroth in string 3:", upreg_damage$Urothelial[uroth_cl3_up], "\n")

cat("Endoth in string 1:", upreg_damage$Endothelial[endo_cl1_up], "\n")
cat("Endoth in string 2:", upreg_damage$Endothelial[endo_cl2_up], "\n")
cat("Endoth in string 3:", upreg_damage$Endothelial[endo_cl3_up], "\n")

cat("Stromal in string 1:", upreg_damage$Stromal[strom_cl1_up], "\n")
cat("Stromal in string 2:", upreg_damage$Stromal[strom_cl2_up], "\n")
cat("Stromal in string 3:", upreg_damage$Stromal[strom_cl3_up], "\n")

cat("Immune in string 1:", upreg_damage$Immune[immune_cl1_up], "\n")
cat("Immune in string 2:", upreg_damage$Immune[immune_cl2_up], "\n")
cat("Immune in string 3:", upreg_damage$Immune[immune_cl3_up], "\n")

cat("SMCs in string 1:", upreg_damage$SMCs[SMCs_cl1_up], "\n")
cat("SMCs in string 2:", upreg_damage$SMCs[SMCs_cl2_up], "\n")
cat("SMCs in string 3:", upreg_damage$SMCs[SMCs_cl3_up], "\n")


cat("Urothelial in string 1:",unlist(strsplit(head(list_DEGs_celltype_vs_condition_up$Urothelial), " "), "\n"))


# Logical vector indicating which genes belong to each string group
uroth_cl1_down <- downreg_damage$Urothelial %in% unlist(strsplit(kmeans_cl1, " "))
uroth_cl2_down <- downreg_damage$Urothelial %in% unlist(strsplit(kmeans_cl2, " "))
uroth_cl3_down <- downreg_damage$Urothelial %in% unlist(strsplit(kmeans_cl3, " "))

endo_cl1_down <- downreg_damage$Endothelial %in% unlist(strsplit(kmeans_cl1, " "))
endo_cl2_down <- downreg_damage$Endothelial %in% unlist(strsplit(kmeans_cl2, " "))
endo_cl3_down <- downreg_damage$Endothelial %in% unlist(strsplit(kmeans_cl3, " "))

strom_cl1_down <- downreg_damage$Stromal %in% unlist(strsplit(kmeans_cl1, " "))
strom_cl2_down <- downreg_damage$Stromal %in% unlist(strsplit(kmeans_cl2, " "))
strom_cl3_down <- downreg_damage$Stromal %in% unlist(strsplit(kmeans_cl3, " "))

immune_cl1_down <- downreg_damage$Immune %in% unlist(strsplit(kmeans_cl1, " "))
immune_cl2_down <- downreg_damage$Immune %in% unlist(strsplit(kmeans_cl2, " "))
immune_cl3_down <- downreg_damage$Immune %in% unlist(strsplit(kmeans_cl3, " "))

SMCs_cl1_down <- downreg_damage$SMCs %in% unlist(strsplit(kmeans_cl1, " "))
SMCs_cl2_down <- downreg_damage$SMCs %in% unlist(strsplit(kmeans_cl2, " "))
SMCs_cl3_down <- downreg_damage$SMCs %in% unlist(strsplit(kmeans_cl3, " "))

# Print the genes that belong to each string group
cat("Uroth in string 1:", downreg_damage$Urothelial[uroth_cl1_down], "\n")
cat("Uroth in string 2:", downreg_damage$Urothelial[uroth_cl2_down], "\n")
cat("Uroth in string 3:", downreg_damage$Urothelial[uroth_cl3_down], "\n")

cat("Endoth in string 1:", downreg_damage$Endothelial[endo_cl1_down], "\n")
cat("Endoth in string 2:", downreg_damage$Endothelial[endo_cl2_down], "\n")
cat("Endoth in string 3:", downreg_damage$Endothelial[endo_cl3_down], "\n")

cat("Stromal in string 1:", downreg_damage$Stromal[strom_cl1_down], "\n")
cat("Stromal in string 2:", downreg_damage$Stromal[strom_cl2_down], "\n")
cat("Stromal in string 3:", downreg_damage$Stromal[strom_cl3_down], "\n")

cat("Immune in string 1:", downreg_damage$Immune[immune_cl1_down], "\n")
cat("Immune in string 2:", downreg_damage$Immune[immune_cl2_down], "\n")
cat("Immune in string 3:", downreg_damage$Immune[immune_cl3_down], "\n")

cat("SMCs in string 1:", downreg_damage$SMCs[SMCs_cl1_down], "\n")
cat("SMCs in string 2:", downreg_damage$SMCs[SMCs_cl2_down], "\n")
cat("SMCs in string 3:", downreg_damage$SMCs[SMCs_cl3_down], "\n")