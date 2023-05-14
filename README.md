# Transcriptomic reprogramming during bladder regenerative response to chemically-induced damage

This repository contains the code for the analysis of RNA sequencing data (bulkRNAseq and scRNAseq) from a study on bladder regeneration in mice. The RNA sequencing data was obtained from bulkRNAseq of bladder tissue samples collected at 0h, 12h, 24h, 72h, 7 days and 30 days after cyclophosphamide injection. The analysis includes quality control, alignment, quantification, differential expression, functional annotation, and cell type deconvolution. In addition, we also analyzed single-cell RNA sequencing data from a previous thesis carried out in our group.

## Table of contents
* [bulkRNAseq DGE analysis](#bulk_DESeq2_DGE.R)
* [bulkRNAseq enrichGO and GSEA functional analysis](#Functional_analysis.R)
* [code for GSEA heatmaps representation](#GSEA_heatmaps.R)
* [scRNAseq complete analysis](#scRNAseq_complete_analysis.R)
* [python CIBERSORTx input conversion script](#CIBERSORTx_conversion_script.py)
