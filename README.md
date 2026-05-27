

# MYC-Glucocorticoid Transcriptional Reprogramming Pipeline

## Overview
RNA-seq analysis pipeline investigating how MYC oncogene status 
modulates the transcriptional response to glucocorticoid drugs 
(Budesonide and Prednisolone) in Double Hit Lymphoma.

**Species:** Mus musculus  
**Drugs:** Budesonide, Prednisolone  
**Timepoints:** 2h, 4h, 8h, 16h  
**Conditions:** MYC ON / MYC OFF (OHT-inducible system)  

## Pipeline Steps

| Step | Script | Description |
|------|--------|-------------|
| 1 | `1_DEAnalysis/1_DEAnalysis.R` | DESeq2 differential expression (RAW + shrunk LFC) |
| 2 | `2_interaction/2_interaction.R` | Drug × MYC interaction model |
| 3 | `3_enrichment/3_enrichment.R` | GSEA + ORA (HALLMARK, KEGG, WikiPathways, C2:CGP) |
| 4 | `4_categories/4_categories.R` | Gene categorization (enhanced, suppressed, independent, additive) |
| 5 | `5_ORA_categories/5_ORA_Categories.R` | ORA per category + heatmaps |
| 6 | `6_scatterplots/6_scatterplots.R` | LFC scatter plots per category |
| 7 | `7_DEG_heatmap/8_DEG_heatmap.R` | DEG heatmap (Pearson clustering on drug/DMSO) |
| 8 | `8_TCGA_survival/survival_analysis.R` | TCGA-BRCA survival analysis |

## Requirements

```r
BiocManager::install(c(
  "DESeq2", "apeglm", "fgsea", "clusterProfiler",
  "msigdbr", "org.Mm.eg.db", "ComplexUpset",
  "TCGAbiolinks", "survminer", "survival"
))

install.packages(c(
  "tidyverse", "pheatmap", "ggplot2", "patchwork"
))
```

## Repository Structure

```txt
.
├── data/
│   ├── metadata.txt
│   └── annotation.txt
│
└── scripts/
    ├── 1_DEAnalysis/
    │   └── 1_DEAnalysis.R
    ├── 2_interaction/
    │   └── 2_interaction.R
    ├── 3_enrichment/
    │   └── 3_enrichment.R
    ├── 4_categories/
    │   └── 4_categories.R
    ├── 5_ORA_categories/
    │   └── 5_ORA_Categories.R
    ├── 6_scatterplots/
        └── 6_scatterplots.R
 ```  

> **Note:** The `gene_counts.tsv` raw counts matrix is not tracked in this repository due to file size. 
> It should be placed in `data/` before running the pipeline.


## Input Files
- `gene_counts.tsv` — raw gene counts matrix
- `metadata.txt` — sample metadata (SampleID, group, drug, OHT, Batch)
- `annotation.txt` — gene annotation (gene_id, symbol)

## Author
Sina Ravanbakhsh — Research Fellow, European Institute of Oncology (IEO), Milan  
Amati Lab | MYC & Cancer Biology
