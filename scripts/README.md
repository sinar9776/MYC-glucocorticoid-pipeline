
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
