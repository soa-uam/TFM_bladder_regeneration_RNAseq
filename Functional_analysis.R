setwd("/home/sozaez/Desktop/TFM/bulkRNAseq/Analyses/Uroth_reg_bulkRNAseq/")

# remotes::update_packages("cli")
library(dplyr)
library(devtools)
library(factoextra)
library(ggplot2)
library(ConsensusClusterPlus)


######################################################################
#################### FUNCTIONAL ANALYSIS #############################
######################################################################

library(org.Mm.eg.db)
library(DOSE)
library(tidyverse)
library(AnnotationHub)
library(ensembldb)
library(RcppEigen)
library(BSgenome.Mmusculus.UCSC.mm39)
library(annotables) # genome version mm38... We need mm39
library(clusterProfiler)
library(pathview)
library(enrichplot)
library(GOSemSim)
library(biomaRt)
library(EnsDb.Mmusculus.v79)


#### 1. Running clusterProfiler
## Explore the grcm38 table loaded by the annotables library
grcm38

# This database doesn't correspond to the last version that we are using.
# Being so, we retrieve our updated annotation data as follows:
ah <- AnnotationHub()
q <- query(ah, "OrgDb")
id <- q$ah_id[length(q)]
id <- q$ah_id[[11]] # 11 org.Mm.eg.db.sqlite
Mm <- ah[[id]]

# Mm dataset info
keytypes(Mm)
columns(Mm)

orgdb <- query(ah, c("mm39", "maintainer@bioconductor.org"))[[1]]
keytypes(orgdb)
columns(orgdb)
head(keys(orgdb, keytype="GENEID"))
ls(Mm)
help("ENTREZID")

egid <- head(keys(orgdb, "GENEID"))
egid <- head(keys(orgdb, "ENTREZID"))
egid <- head(keys(orgdb, "ENTREZID"))
egid <- head(keys(orgdb, "ENTREZID"))

select(orgdb, egid, c("SYMBOL", "GENENAME"), "ENTREZID")

# load data
biomart_ids <- read.table("data/biomart_ids.txt", header = TRUE, row.names = 1)
res_lrt_tb <- read.table(file = "results/res_lrt_tb.xls", header = TRUE) # res_LRT_0.01 3839 significant genes, but here are all loaded
clustered_genes3 <- read.table("data/clustered_genes3.xls", header = TRUE)
# 3 kmeans clusters from LRT_0.01
res_ids1 <- read.table("results/k3_clust1_enrichGO.xls", header = TRUE)
res_ids2 <- read.table("results/k3_clust2_enrichGO.xls", header = TRUE)
res_ids3 <- read.table("results/k3_clust3_enrichGO.xls", header = TRUE)

res_12h <- read.table("results/sigDE12h.xls", header = TRUE)
res_24h <- read.table("results/sigDE24h.xls", header = TRUE)
res_72h <- read.table("results/sigDE72h.xls", header = TRUE)
res_7d <- read.table("results/sigDE7d.xls.xls", header = TRUE)
res_30d <- read.table("results/sigDE30d.xls", header = TRUE)

res_24h12h <- read.table("results/sigDE24h12h.xls", header = TRUE)
res_72h24h <- read.table("results/sigDE72h24h.xls", header = TRUE)
res_7d72h <- read.table("results/sigDE7d72h.xls", header = TRUE)
res_30d7d <- read.table("results/sigDE30d7d.xls", header = TRUE)
wald_res <- read.table("results/wald_res.xls", header = TRUE)
# For a faster executing, load data and skip to line 123 :)


################### LRT - FUNCTIONAL ANALYSIS ###################

# Extract Cluster groups lists:
num_clust1 <- which(clustered_genes3$cluster == "cluster 1")
num_clust2 <- which(clustered_genes3$cluster == "cluster 2")
num_clust3 <- which(clustered_genes3$cluster == "cluster 3")

genes1 <- rownames(clustered_genes3)[num_clust1]
genes2 <- rownames(clustered_genes3)[num_clust2]
genes3 <- rownames(clustered_genes3)[num_clust3]

mgi <- which(biomart_ids$mgi_symbol %in% genes1) # repeat for each cluster
enb <- which(biomart_ids$ensembl_gene_id %in% genes1)
bio_ens1 <- biomart_ids$ensembl_gene_id[mgi]
bio_ens2 <- biomart_ids$ensembl_gene_id[enb]

u1 <- union(bio_ens1, bio_ens2)

## Return the IDs for the gene symbols in the DE results
idx_cl1 <- grcm38$ensgene %in% u1
idx_cl2 <- grcm38$ensgene %in% u2
idx_cl3 <- grcm38$ensgene %in% u3

ids_cl1 <- grcm38[idx_cl1, ]
ids_cl2 <- grcm38[idx_cl2, ]
ids_cl3 <- grcm38[idx_cl3, ]

# retrieve missing IDs (due to version conflict in annotables package)
cl1_leftgenes <- setdiff(u1, ids_cl1$ensgene)
cl2_leftgenes <- setdiff(u2, ids_cl2$ensgene)
cl3_leftgenes <- setdiff(u3, ids_cl3$ensgene)

# remove duplicate IDs prior to assessing enriched GO terms. Gene names can map to more than one Ensembl ID
non_duplicates1 <- which(duplicated(ids_cl1$symbol) == FALSE)
non_duplicates2 <- which(duplicated(ids_cl2$symbol) == FALSE)
non_duplicates3 <- which(duplicated(ids_cl3$symbol) == FALSE)

ids_cl1 <- ids_cl1[non_duplicates1, ]
ids_cl2 <- ids_cl2[non_duplicates2, ]
ids_cl3 <- ids_cl3[non_duplicates3, ]

## Merge the IDs with the results
res_ids1 <- inner_join(res_lrt_tb, ids_cl1, by=c("gene"="symbol"))
res_ids2 <- inner_join(res_lrt_tb, ids_cl2, by=c("gene"="symbol"))
res_ids3 <- inner_join(res_lrt_tb, ids_cl3, by=c("gene"="symbol"))

# write.table(res_ids3, file = "results/k3_clust3_enrichGO.xls", sep="\t", row.names = TRUE)

## Create background dataset for hypergeometric testing using all genes tested for significance in the results
e <- which(biomart_ids$mgi_symbol %in% res_lrt_tb$gene)
a <- which(biomart_ids$ensembl_gene_id %in% res_lrt_tb$gene)
i <- biomart_ids$ensembl_gene_id[e]
o <- biomart_ids$ensembl_gene_id[a]
u <- union(i, o)

allLRT_genes <- as.character(u)
idx_all <- grcm38$ensgene %in% u
ids_all <- grcm38[idx_all, ]
allLRT_leftgenes <- setdiff(u, ids_all$ensgene)
non_duplicates_all <- which(duplicated(ids_all$symbol) == FALSE)
ids_all <- ids_all[non_duplicates_all, ]
res_all <- inner_join(res_lrt_tb, ids_all, by=c("gene"="symbol"))
allLRT_genes <- as.character(res_all$ensgene)
allLRT_genes <- na.omit(allLRT_genes)
allLRT_genes <- union(allLRT_genes, allLRT_leftgenes) # aparecen 51743 genes de los 56980 genes iniciales... usar Biomart mejor para hacer la anotación

## Extract significant results
sigLRT_genes1 <- as.character(res_ids1$ensgene)
sigLRT_genes1 <- union(sigLRT_genes1, cl1_leftgenes)
sigLRT_genes2 <- as.character(res_ids2$ensgene)
sigLRT_genes2 <- union(sigLRT_genes2, cl2_leftgenes)
sigLRT_genes3 <- as.character(res_ids3$ensgene)
sigLRT_genes3 <- union(sigLRT_genes3, cl3_leftgenes)

## Run GO enrichment analysis 

# INCLUIR EJECUCIÓN DATASET Mm - primero guardar database en formato correcto...
ego1 <- enrichGO(gene = sigLRT_genes1, 
                     universe = allLRT_genes,
                     keyType = "ENSEMBL",
                     OrgDb = Mm, # org.Mm.eg.db used in tutorial, but not updated
                     ont = "BP", 
                     pAdjustMethod = "BH", 
                     qvalueCutoff = 0.05, 
                     readable = TRUE)

## Output results from GO analysis to a table
cluster_summary1 <- data.frame(ego1)

# write.csv(cluster_summary, "results/clusterProfiler_cluster3.xls")


#### 2. VISUALIZING clusterProfiler results
## Dotplot 
dotplot(ego1, showCategory=50, color = "p.adjust", title = "Cluster 1")
dotplot(ego2, showCategory=50, color = "p.adjust", title = "Cluster 2")
dotplot(ego3, showCategory=50, color = "p.adjust", title = "Cluster 3")

# pdf("plots/dotplot_cluster1.pdf", width=12, height=8)
# dotplot(ego1_test, showCategory=20, color = "p.adjust", title = "Cluster 1")
# dev.off()

# GeneRatio = genes of interest in the gene set / total genes of interest

# top 10 sig genes
pos_gene <- which(row.names(res_lrt) %in% genes3)
k1 <- res_lrt_tb[pos_gene,] %>% arrange(padj)
top10_k3 <- head(k1, 10)

## Enrichmap clusters the 50 most significant (by padj) GO terms to visualize relationships between terms
d <- godata('org.Mm.eg.db', ont="BP")
d_test <- godata(OrgDb = Mm, ont="BP")
  
ego_p2 <- pairwise_termsim(ego2_test, method="Wang", semData = d)
ego_p2_JC <- pairwise_termsim(ego2_test, method="Resnik", semData = d)
ego_p2_Lin <- pairwise_termsim(ego2_test, method="Lin", semData = d)
ego_p2_Rel <- pairwise_termsim(ego2_test, method="Rel", semData = d)
ego_p2_Jiang <- pairwise_termsim(ego2_test, method="Jiang", semData = d)

ego_p2_JC <- pairwise_termsim(ego2_test, method="JC", semData = d)
ego_p2_test <- pairwise_termsim(ego2_test, method="Wang", semData = d_test)

emapplot(ego_p2, showCategory = 50)
emapplot(ego_p2_test, showCategory = 50)
emapplot(ego_p2_JC, showCategory = 50)
emapplot(ego_p2_Lin, showCategory = 50)
emapplot(ego_p2_Rel, showCategory = 50)
emapplot(ego_p2_Jiang, showCategory = 50)

treeplot(ego_p2) # varía un poco dependiendo del método de medición usado... stick to Wang


# ridgeplot(ego1)
# upsetplot(ego1)

# terms <- ego1$Description[1:5]
# p <- pmcplot(terms, 2010:2020)
# p2 <- pmcplot(terms, 2010:2020, proportion=FALSE)
# plot_grid(p, p2, ncol=2)


# CATEGORY NETPLOT
## To color genes by log2 fold changes, we need to extract the log2 fold changes from our results table creating a named vector

# Take just sig Wald genes, not the whole list. 
# Get allDEGs indexes - identities of wald_res table

DEGswald <- read.table("results/DEGs_wald.xls", header = TRUE)

post_wald_genes <- which(DEGswald$all_DEGs %in% row.names(wald_res))
sig_wald_log <- wald_res[post_wald_genes,]

post_cluster1_genes <- which(genes1 %in% row.names(wald_res))
sig_cluster1_log <- wald_res[post_cluster1_genes,]
cluster1_foldchanges <- sig_cluster1_log$log2FoldChange
names(cluster1_foldchanges) <- sig_cluster1_log$gene

post_cluster2_genes <- which(genes2 %in% row.names(wald_res))
sig_cluster2_log <- wald_res[post_cluster2_genes,]
cluster2_foldchanges <- sig_cluster2_log$log2FoldChange
names(cluster2_foldchanges) <- sig_cluster2_log$gene

post_cluster3_genes <- which(genes3 %in% row.names(wald_res))
sig_cluster3_log <- wald_res[post_cluster3_genes,]
cluster3_foldchanges <- sig_cluster3_log$log2FoldChange
names(cluster3_foldchanges) <- sig_cluster3_log$gene


e <- which(biomart_ids$mgi_symbol %in% row.names(sig_wald_log))
a <- which(biomart_ids$ensembl_gene_id %in% row.names(sig_wald_log)) #No hay
i <- biomart_ids$ensembl_gene_id[e]

# para los clusters de kmeans
e <- which(biomart_ids$mgi_symbol %in% row.names(sig_cluster1_log))
a <- which(biomart_ids$ensembl_gene_id %in% row.names(sig_cluster1_log)) #No hay
i <- biomart_ids$ensembl_gene_id[e]
idx_wald <- grcm38$ensgene %in% i

ids_wald <- grcm38[idx_wald, ]
wald_leftgenes <- setdiff(i, ids_wald$ensgene)

non_duplicates_wald <- which(duplicated(ids_wald$symbol) == FALSE)
ids_wald <- ids_wald[non_duplicates_wald, ]
wald_res$gene <- row.names(wald_res)
res_ids_wald <- inner_join(wald_res, ids_wald, by=c("gene"="symbol"))

sigWALD_genes <- as.character(res_ids_wald$ensgene)
sigWALD_genes <- union(sigWALD_genes, wald_leftgenes)

sig_wald_foldchanges <- sig_wald_log$log2FoldChange
sig_wald_log["gene"] <- row.names(sig_wald_log)
names(sig_wald_foldchanges) <- sig_wald_log$gene

ego_wald <- enrichGO(gene = sigWALD_genes, 
                      universe = allLRT_genes,
                      keyType = "ENSEMBL",
                      OrgDb = Mm, 
                      ont = "BP", 
                      pAdjustMethod = "BH", 
                      qvalueCutoff = 0.05, 
                      readable = TRUE)

heatplot(ego_wald, foldChange = sig_wald_foldchanges, showCategory = 5)

dotplot(ego_wald, showCategory=20, color = "p.adjust") # for this, better make up and down regulated separatelly

## Cnetplot details the genes associated with one or more terms - by default gives the top 5 significant terms (by padj)
cnetplot(ego1_test, 
         categorySize = "p.adjust", 
         showCategory = 5, 
         foldChange=cluster1_foldchanges, 
         vertex.label.font=6,
         layout = "fr")

## If some of the high fold changes are getting drowned out due to a large range, you could set a maximum fold change value
sig_wald_foldchanges <- ifelse(sig_wald_foldchanges > 2, 2, sig_wald_foldchanges)
sig_wald_foldchanges <- ifelse(sig_wald_foldchanges < -2, -2, sig_wald_foldchanges)

# If you are interested in significant processes that are not among the top five, you can subset your ego dataset to only display these processes:

## Subsetting the ego results without overwriting original `ego` variable
egowald_subset <- ego_wald

egowald_subset@result <- ego_wald@result[c(1,3,4,8,9),]

## Plotting terms of interest
cnetplot(egowald_subset, 
         categorySize="pvalue", 
         foldChange=sig_wald_foldchanges, 
         showCategory = 5,
         vertex.label.font=6)


# Also try:
# categorys <- c("pre-malignant neoplasm", "intestinal disease",
#                "breast ductal carcinoma", "non-small cell lung carcinoma")
# cnetplot(x2, showCategory = categorys)





################### WALD TEST - FUNCTIONAL ANALYSIS ################### 

res_24h12h <- read.table("results/sigDE24h12h.xls", header = TRUE)
res_72h24h <- read.table("results/sigDE72h24h.xls", header = TRUE)
res_7d72h <- read.table("results/sigDE7d72h.xls", header = TRUE)
res_30d7d <- read.table("results/sigDE30d7d.xls", header = TRUE)
res_12h <- read.table("results/sigDE12h.xls", header = TRUE)
res_24h <- read.table("results/sigDE24h.xls", header = TRUE)
res_72h <- read.table("results/sigDE72h.xls", header = TRUE)
wald_res <- read.table("results/wald_res.xls", header = TRUE)


lfc.cutoff <- 0

res_12h_up <- res_12h %>%
  dplyr::filter(log2FoldChange > lfc.cutoff) %>% 
  arrange(padj)
res_12h_down <- res_12h %>%
  dplyr::filter(log2FoldChange < lfc.cutoff) %>% 
  arrange(padj)

res_24h_up <- res_24h %>%
  dplyr::filter(log2FoldChange > lfc.cutoff) %>% 
  arrange(padj)
res_24h_down <- res_24h %>%
  dplyr::filter(log2FoldChange < lfc.cutoff) %>% 
  arrange(padj)

res_72h_up <- res_72h %>%
  dplyr::filter(log2FoldChange > lfc.cutoff) %>% 
  arrange(padj)
res_72h_down <- res_72h %>%
  dplyr::filter(log2FoldChange < lfc.cutoff) %>% 
  arrange(padj)

res_24h12h_up <- res_24h12h %>%
  dplyr::filter(log2FoldChange > lfc.cutoff) %>% 
  arrange(padj)
res_24h12h_down <- res_24h12h %>%
  dplyr::filter(log2FoldChange < lfc.cutoff) %>% 
  arrange(padj)

res_72h24h_up <- res_72h24h %>%
  dplyr::filter(log2FoldChange > lfc.cutoff) %>% 
  arrange(padj)
res_72h24h_down <- res_72h24h %>%
  dplyr::filter(log2FoldChange < lfc.cutoff) %>% 
  arrange(padj)

res_7d72h_up <- res_7d72h %>%
  dplyr::filter(log2FoldChange > lfc.cutoff) %>% 
  arrange(padj)
res_7d72h_down <- res_7d72h %>%
  dplyr::filter(log2FoldChange < lfc.cutoff) %>% 
  arrange(padj)

res_30d7d_up <- res_30d7d %>%
  dplyr::filter(log2FoldChange > lfc.cutoff) %>% 
  arrange(padj)
res_30d7d_down <- res_30d7d %>%
  dplyr::filter(log2FoldChange < lfc.cutoff) %>% 
  arrange(padj)


biomart_ids <- read.table("data/biomart_ids.txt", header = TRUE, row.names = 1)

wald_res$gene <- rownames(wald_res)

e <- which(biomart_ids$mgi_symbol %in% res_72h24h_down$gene)
a <- which(biomart_ids$ensembl_gene_id %in% res_72h24h_down$gene)
i <- biomart_ids$ensembl_gene_id[e]
o <- biomart_ids$ensembl_gene_id[a]

u_12 <- union(i, o)

ids_cl1 <- grcm38$ensgene %in% u_12
ids_cl1 <- grcm38[ids_cl1, ]
leftout_genes <- setdiff(u_12, ids_cl1$ensgene)
non_duplicates1 <- which(duplicated(ids_cl1$symbol) == FALSE)
ids_cl1 <- ids_cl1[non_duplicates1, ]
res_ids1 <- inner_join(wald_res, ids_cl1, by=c("gene"="symbol"))

# e <- which(biomart_ids$mgi_symbol %in% wald_res$gene)
# a <- which(biomart_ids$ensembl_gene_id %in% wald_res$gene)
# i <- biomart_ids$ensembl_gene_id[e]
# o <- biomart_ids$ensembl_gene_id[a]
# u <- union(i, o)
# 
# all_wald_genes <- as.character(u)

# write.table(all_wald_genes, file = "results/all_wald_enrichGO.xls", sep="\t", row.names = TRUE)

## Extract significant results
sig_genes <- as.character(res_ids1$ensgene)
length(sig_genes)
sig_genes <- union(sig_genes, leftout_genes)
# length(sig_genes)
# length(res_12h_up$gene)

# write.table(sig_genes, file = "results/72h_enrichGO.xls", sep="\t", row.names = TRUE)

# res_ids24h12h <- read.table("results/24h12h_enrichGO.xls", header = TRUE)
# res_ids72h24h <- read.table("results/72h24h_enrichGO.xls", header = TRUE)
# res_ids7d72h <- read.table("results/7d72h_enrichGO.xls", header = TRUE)
# res_ids30d7d <- read.table("results/30d7d_enrichGO.xls", header = TRUE)
# res_ids12h <- read.table("results/12h_enrichGO.xls", header = TRUE)
# res_ids24h <- read.table("results/24h_enrichGO.xls", header = TRUE)
# res_ids72h <- read.table("results/72h_enrichGO.xls", header = TRUE)
# 
# res_ids24h12h <- as.character(res_ids24h12h$x)
# res_ids72h24h <- as.character(res_ids72h24h$x)
# res_ids12h <- as.character(res_ids12h$x)
# res_ids24h <- as.character(res_ids24h$x)
# res_ids72h <- as.character(res_ids72h$x)

# all_wald_genes <- read.table("results/all_wald_enrichGO.xls", header = TRUE)
# all_wald_genes <- as.character(all_wald_genes$x)

ego72h24h_down <- enrichGO(gene = sig_genes, 
                           universe = all_wald_genes,
                           keyType = "ENSEMBL",
                           OrgDb = Mm, 
                           ont = "BP", # Biological Process
                           pAdjustMethod = "BH", # Benjamini-Hochberg
                           qvalueCutoff = 0.05, 
                           readable = TRUE)

cluster_summary <- data.frame(ego12h_up)
# write.table(cluster_summary, file = "results/cluster_summary30d7d.xls", sep="\t", row.names = TRUE)

dotplot(ego12h_up, showCategory=20, color = "p.adjust", title = "12h vs 0h up")
dotplot(ego12h_down, showCategory=20, color = "p.adjust", title = "12h vs 0h down")

dotplot(ego24h_up, showCategory=20, color = "p.adjust", title = "24h vs 0h up")
dotplot(ego24h_down, showCategory=20, color = "p.adjust", title = "24h vs 0h down")

dotplot(ego72h_up, showCategory=20, color = "p.adjust", title = "72h vs 0h up")
dotplot(ego72h_down, showCategory=20, color = "p.adjust", title = "72h vs 0h down")

dotplot(ego24h12h_up, showCategory=20, color = "p.adjust", title = "24h vs 12h up")
dotplot(ego24h12h_down, showCategory=20, color = "p.adjust", title = "24h vs 12h down")

dotplot(ego72h24h_up, showCategory=20, color = "p.adjust", title = "72h vs 24h up")
dotplot(ego72h24h_down, showCategory=20, color = "p.adjust", title = "72h vs 24h down")

# pdf("plots/GSEA/24h_up.pdf", width=8, height=9)
# dotplot(ego24h_up, showCategory=20, color = "p.adjust", title = "24h vs 0h up")
# dev.off()

