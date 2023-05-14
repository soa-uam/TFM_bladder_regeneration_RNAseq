# Transcriptomic reprogramming during bladder regenerative response to chemically-induced damage

This repository contains the code for the analysis of RNA sequencing data (bulkRNAseq and scRNAseq) from a study on bladder regeneration in mice. The RNA sequencing data was obtained from bulkRNAseq of bladder tissue samples collected at 0h, 12h, 24h, 72h, 7 days and 30 days after cyclophosphamide injection. The analysis includes quality control, alignment, quantification, differential expression, functional annotation, and cell type deconvolution. In addition, we also analyzed single-cell RNA sequencing data from a previous thesis carried out in our group.

## Table of contents
* [bulkRNAseq DGE analysis](#bulk_DESeq2_DGE.R)
* [bulkRNAseq GSEA](#Functional_analysis.R)
* [GSEA heatmaps representation](#GSEA_heatmaps.R)
* [scRNAseq analysis](#scRNAseq_complete_analysis.R)
* [python CIBERSORTx input conversion script](#CIBERSORTx_conversion_script.py)


## bulkRNA-seq DGE analysis
- `bulk_DESeq2_DGE.R`: R script that uses the DESeq2 package to analyze the bulkRNA-seq data from the study on bladder regeneration in mice. It performs differential expression analysis and uses clustering algorithms to identify genetic patterns.


## bulkRNA-seq GSEA
- `Functional_analysis.R`: R script that uses the clusterProfiler package to perform functional analysis of the differentially expressed genes from the bulkRNAseq data. It performs enrichment analysis of gene ontology terms and generates dot plots of the results.


## GSEA heatmaps representation
- `GSEA_heatmaps.R`: R script that processes the clusterProfiler and GSEA software pathway results from the differential expression analysis of the bulkRNAseq data. It creates heatmaps of their gene expression across the different time points for selectet pathways.


## scRNAseq analysis
- `scRNAseq_complete_analysis.R`: R script that uses the Seurat package to analyze the scRNA-seq integrated bladder data from a previous work in our group. This script performs cell annotation, damage markers identification, validation of bulkRNA-seq results and preprocessing of data for CIBERSORTx input.


## CIBERSORTx input conversion script
- `CIBERSORTx_conversion_script.py`: Python script that processes the output data from Seurat "barcode and celltype label" and "counts per celltype" to create a signature matrix for CIBERSORTx. It converts the single-cell barcodes to cell labels and writes a new file with gene counts per cell type. The output file has cell types as columns and genes as rows.





