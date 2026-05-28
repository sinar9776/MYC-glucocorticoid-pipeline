
## Pipeline Steps

| Step | Script | Description |
|------|--------|-------------|
| 1 | `1_DEAnalysis/1_DEAnalysis.R` | DESeq2 differential expression (RAW + shrunk LFC) |
| 2 | `2_interaction/2_interaction.R` | Drug × MYC interaction model |
| 3 | `3_enrichment/3_enrichment.R` | GSEA + ORA (HALLMARK, KEGG, WikiPathways, C2:CGP) |
| 4 | `4_categories/4_categories.R` | Gene categorization (enhanced, suppressed, independent, additive) |
| 5 | `5_ORA_categories/5_ORA_Categories.R` | ORA per category + heatmaps |
| 6 | `6_scatterplots/6_scatterplots.R` | LFC scatter plots per category |


## Requirements

```r
BiocManager::install(c(
  "DESeq2", "apeglm", "fgsea", "clusterProfiler",
  "msigdbr", "org.Mm.eg.db", "ComplexUpset"
))

install.packages(c(
  "tidyverse", "pheatmap", "ggplot2", "patchwork"
))
```
