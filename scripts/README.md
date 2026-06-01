
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

# Pipeline Explaination

## 1. Differential Expression Analysis

This script performs differential gene expression analysis using DESeq2 across 
all experimental conditions. Two glucocorticoid drugs (Budesonide and 
Prednisolone) were tested at four timepoints (2h, 4h, 8h, 16h) in two MYC 
states (ON and OFF, controlled by an OHT-inducible system), against a DMSO 
vehicle control.

**Contrasts generated:**
- Drug vs DMSO (baseline drug effect)
- Drug+OHT vs OHT (drug effect with MYC off)
- Drug+OHT vs DMSO (combined effect)
- DMSO+OHT vs DMSO (pure MYC-off effect)

**Two result types are produced for each contrast:**
- **RAW**: standard DESeq2 results with Wald test statistics
- **SHRUNK**: log2 fold changes shrunk with apeglm to reduce noise from 
  lowly expressed genes, recommended for visualization and ranking

**Outputs:**
- Full result tables (`tT_*.tsv`), upregulated (`up_*.tsv`) and 
  downregulated (`down_*.tsv`) gene lists for each contrast
- Dispersion plots, PCA plots per comparison, MA plots, and volcano plots
- Expressed gene universe file used as background for all enrichment analyses

**Key parameters:**
- Adjusted p-value cutoff: 0.01
- Batch effect correction included in the design formula
- LFC shrinkage method: apeglm

## 2. Drug × MYC Interaction Analysis

This script identifies genes whose response to glucocorticoid treatment 
is significantly different depending on MYC status — in other words, genes 
where MYC actively modifies the drug response rather than acting independently.

**Statistical model:**

For each drug-timepoint combination, a DESeq2 model is fitted on four groups:
DMSO, DMSO+OHT, Drug, Drug+OHT using the design:

~ Batch + drug + OHT + drug:OHT

The interaction term `drug:OHT` captures the deviation from additivity — 
a significant interaction means the drug effect in the MYC-off condition 
is not simply the sum of the drug effect and the MYC-off effect separately.

**Interpretation:**
- Positive interaction: drug response is amplified when MYC is off
- Negative interaction: drug response is dampened when MYC is off

**Outputs:**
- Interaction result tables (`tT_int_*.tsv`) per drug-timepoint
- Upregulated and downregulated interaction gene lists
- MA plots and volcano plots for each comparison

This analysis forms the statistical backbone for the gene categorization 
in Script 4, where interaction results are combined with the main 
effect results to classify gene behavior.

## 3. Functional Enrichment Analysis

This script performs pathway and gene set enrichment analysis on the 
differential expression results from Scripts 1 and 2, identifying 
which biological processes are transcriptionally regulated by 
glucocorticoids and how MYC status modulates these programs.

**Two complementary approaches are used:**

**GSEA (Gene Set Enrichment Analysis):**  
Ranks all expressed genes by their DESeq2 test statistic and tests 
whether genes belonging to a pathway are collectively enriched at 
the top (activated) or bottom (repressed) of the ranked list. 
Does not require an arbitrary significance cutoff and captures 
subtle coordinated changes across a pathway.

**ORA (Over-Representation Analysis):**  
Tests whether significantly differentially expressed genes 
(from the up/down gene lists) overlap with pathway gene sets 
more than expected by chance, using a hypergeometric test.

**Gene set collections tested:**
- HALLMARK (50 well-defined biological signatures)
- GO Biological Process (GOBP)
- KEGG pathways
- MSigDB C2:CP (curated canonical pathways)

**Outputs:**
- GSEA and ORA result tables per contrast and gene set collection
- Bar plots showing top enriched pathways per contrast
- Structured NES heatmaps summarizing pathway activity across all 
  drug-timepoints simultaneously (3 heatmap layouts per drug per collection)

## 4. Gene Categorization

This script classifies differentially expressed genes into functional 
categories based on how their response to glucocorticoid treatment 
changes with MYC status. This is the central analytical framework 
of the project, translating statistical results into biologically 
interpretable gene classes.

**Categories defined:**

<img width="866" height="466" alt="Picture 1" src="https://github.com/user-attachments/assets/06c8a4ff-8866-4bf7-ba50-e1cb4d6dfb4b" />
<img width="549" height="454" alt="Screenshot 2026-06-01 at 14 23 43" src="https://github.com/user-attachments/assets/454fd2e3-cab9-40f6-b9c8-05fc34936c29" />



| Category | Definition |
|----------|-----------|
| **Enhanced** | Gene is differentially expressed in both MYC-on and MYC-off conditions, with a stronger effect when MYC is off. The interaction term confirms amplification. |
| **Suppressed** | Gene responds to the drug but the effect is reduced or lost when MYC is off. MYC is required for the full drug response. |
| **Switched** | Gene changes direction between MYC-on and MYC-off conditions (e.g. upregulated with MYC, downregulated without). |
| **Independent** | Gene responds consistently to the drug regardless of MYC status, with no significant interaction. MYC does not modulate this response. |
| **Additive** | Subset of independent genes that are also regulated by OHT alone, suggesting both MYC loss and drug treatment independently contribute to the same direction of change. |

**Logic:**
Categories are assigned by combining the significance and direction of 
three contrasts per gene: drug/DMSO, drug+OHT/OHT, and the interaction 
term from Script 2. Higher timepoints take priority in the assignment.

**Outputs:**
- UpSet plots showing overlap between gene lists across conditions
- Per-category gene tables with ENSEMBL IDs and gene symbols
- Additive gene lists (up and down) per drug-timepoint

## 5. Over-Representation Analysis per Gene Category

This script runs functional enrichment analysis specifically on the 
gene categories defined in Script 4, asking what biological processes 
characterize each class of MYC-modulated glucocorticoid response.

**Gene set collections:**
- HALLMARK
- KEGG
- WikiPathways (WP)
- C2:CGP (Chemical and Genetic Perturbations — useful for identifying 
  known drug signatures and oncogene-related gene sets)

**Category handling:**
- **Enhanced and suppressed** categories are analyzed separately 
  for up and down genes, preserving directionality information
- **Independent and additive** categories are also kept separate 
  since their up/down distinction reflects the direction of drug 
  response and has direct biological meaning

**Heatmap outputs:**
For each gene set collection, a summary heatmap is produced per category 
showing -log10(adjusted p-value) across all drug-timepoints simultaneously 
(Budesonide and Prednisolone at 2h, 4h, 8h, 16h). This allows direct 
visual comparison of how pathway enrichment evolves over time and 
between the two drugs.

**Table outputs:**
Each result table includes the genes involved in each significant pathway 
(both ENSEMBL IDs and gene symbols), making it straightforward to 
identify which specific genes are driving enrichment.

## 6. LFC Scatter Plots

This script produces scatter plots that visualize the relationship 
between gene expression fold changes across pairs of contrasts, 
with genes colored by their assigned category from Script 4. 
These plots provide an intuitive visual validation of the 
categorization logic and reveal the overall structure of 
MYC-glucocorticoid transcriptional interactions.

**Three plot types are generated per drug-timepoint:**

**Plot type 1 — Drug/DMSO vs Drug+OHT/OHT:**  
X-axis: drug effect with MYC on  
Y-axis: drug effect with MYC off  
Points on the diagonal indicate genes that respond equally in both 
MYC states (independent). Points above the diagonal are enhanced 
(stronger response with MYC off); points below are suppressed.

**Plot type 2 — OHT/DMSO vs Drug+OHT/DMSO:**  
X-axis: pure MYC-off effect  
Y-axis: combined drug + MYC-off effect  
Reveals the relationship between baseline MYC loss and the 
combined treatment effect.

**Plot type 3 — Drug/DMSO vs Drug+OHT/DMSO:**  
X-axis: drug effect alone  
Y-axis: combined drug + MYC-off effect  
Shows how MYC loss shifts the drug response at the whole-transcriptome level.

**Each plot type is produced in four versions:**
- All categories (main view)
- Independent genes only
- Shifted baseline genes only
- Additive genes only

Individual per-drug plots and combined panels across all 
drug-timepoints are saved as PDF and high-resolution TIFF.
