## ============================================================
## 2_interaction.R — Interaction analysis (drug-OHT, per timepoint)
## Reference: drug = DMSO | OHT = OFF
## Outputs:
##   tables/: tT_int_<label>.tsv, up_int_<label>.tsv, down_int_<label>.tsv
##   figures/: MA_<label>.pdf/.tiff, Volcano_<label>.pdf/.tiff
## ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(DESeq2)
  library(ggplot2)
})

set.seed(1)
message(">>> 2_interaction: baseline-OHT interaction mode (per timepoint)")

## -------------------- Directories --------------------
setwd("/hpcnfs/data/BA/glucocorticoid_keep/sinaR/m.rezaei")
base_dir    <- getwd()
data_dir    <- file.path(base_dir, "data")
results_dir <- file.path(base_dir, "results", "2_interaction")
tbl_dir     <- file.path(results_dir, "tables")
fig_dir     <- file.path(results_dir, "figures")
dir.create(tbl_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

## -------------------- Parameters --------------------
padj_cutoff <- 0.01
lfc_cutoff  <- 0

## ============================================================
## PART 1 — Interaction DESeq2 Calculation
## ============================================================

analysis_level <- "gene"   # or "isoform"

## ---- Load metadata ----
metadata <- read.delim(file.path(base_dir, "metadata.txt"), check.names = FALSE)
stopifnot(all(c("SampleID", "group", "drug", "OHT", "Batch") %in% colnames(metadata)))

metadata$OHT   <- factor(metadata$OHT, levels = c("OFF", "ON"))
metadata$Batch <- factor(metadata$Batch)
metadata$group <- factor(metadata$group)

## ---- Load counts ----
if (analysis_level == "gene") {
  count_file <- file.path(base_dir, "gene_counts.tsv")
  id_col <- "gene_id"
} else {
  count_file <- file.path(base_dir, "isoform_counts.tsv")
  id_col <- "isoform_id"
}

counts    <- readr::read_tsv(count_file, show_col_types = FALSE) %>% as.data.frame()
rownames(counts) <- counts[[1]]
count_mat <- counts[, intersect(colnames(counts), metadata$SampleID), drop = FALSE]
metadata  <- metadata[match(colnames(count_mat), metadata$SampleID), ]
stopifnot(all(colnames(count_mat) == metadata$SampleID))

## ---- Load annotation ----
annotation <- read.delim(file.path(base_dir, "annotation.txt"), check.names = FALSE)
sym_col <- if ("symbol" %in% names(annotation)) "symbol" else if ("gene_name" %in% names(annotation)) "gene_name" else NA
ann_keep_cols <- unique(c("gene_id", id_col, sym_col))
ann_keep_cols <- ann_keep_cols[ann_keep_cols %in% names(annotation)]
annotation2 <- annotation[, ann_keep_cols, drop = FALSE]
if (!is.na(sym_col)) names(annotation2)[names(annotation2) == sym_col] <- "SYMBOL"
if (analysis_level == "gene") {
  annotation2 <- annotation2 %>%
    group_by(gene_id) %>%
    summarise(SYMBOL = dplyr::first(na.omit(SYMBOL)), .groups = "drop")
}
if (!"gene_id" %in% names(annotation2)) stop("Annotation must include gene_id column.")

## ---- Derive drug-timepoint labels (non-DMSO, OHT==OFF base groups) ----
# e.g. "Budesonide_2", "Budesonide_4", "Prednisolone_8" ...
drug_time_labels <- unique(metadata$group[
  metadata$drug != "DMSO" & metadata$OHT == "OFF"
])
drug_time_labels <- as.character(drug_time_labels)
message("Drug-timepoint groups to process: ", paste(drug_time_labels, collapse = ", "))

## ---- Loop: one DESeq2 per drug-timepoint ----
for (label in drug_time_labels) {
  
  dr  <- unique(metadata$drug[metadata$group == label])   # e.g. "Budesonide"
  oht_label <- paste0(label, "_OHT")                      # e.g. "Budesonide_2_OHT"
  
  if (!oht_label %in% levels(metadata$group)) {
    message("WARNING: No OHT group found for ", label, " — skipping.")
    next
  }
  
  message("\n>>> Processing: ", label)
  
  ## Subset: DMSO + DMSO_OHT + this drug (OFF) + this drug (ON)
  keep      <- metadata$group %in% c("DMSO", "DMSO_OHT", label, oht_label)
  meta_sub  <- metadata[keep, ]
  count_sub <- count_mat[, keep, drop = FALSE]
  
  ## Build clean drug_fac: DMSO vs the drug name
  meta_sub$drug_fac <- factor(
    ifelse(meta_sub$drug == "DMSO", "DMSO", dr),
    levels = c("DMSO", dr)
  )
  
  ## Drop Batch if only one level in this subset
  design_formula <- if (length(unique(meta_sub$Batch)) > 1) {
    ~ Batch + drug_fac + OHT + drug_fac:OHT
  } else {
    ~ drug_fac + OHT + drug_fac:OHT
  }
  message("Design: ", deparse(design_formula))
  
  dds <- DESeqDataSetFromMatrix(
    countData = round(as.matrix(count_sub)),
    colData   = meta_sub,
    design    = design_formula
  )
  dds <- DESeq(dds, parallel = FALSE)
  
  ## Get interaction term
  res_names        <- resultsNames(dds)
  interaction_term <- grep("drug_fac.*\\.OHTON$", res_names, value = TRUE)
  if (length(interaction_term) == 0) {
    message("WARNING: No interaction term found for ", label, " — skipping.")
    next
  }
  message("Interaction term: ", interaction_term)
  
  ## Extract results
  res    <- results(dds, name = interaction_term)
  res_df <- as.data.frame(res)
  res_df[[id_col]] <- sub("\\..*$", "", rownames(res_df))
  res_df <- merge(res_df, annotation2, by = id_col, all.x = TRUE, sort = FALSE)
  
  res_df$sig <- ifelse(
    !is.na(res_df$padj) & res_df$padj < padj_cutoff & res_df$log2FoldChange >  lfc_cutoff, "up",
    ifelse(
      !is.na(res_df$padj) & res_df$padj < padj_cutoff & res_df$log2FoldChange < -lfc_cutoff, "down",
      "ns"
    )
  )
  
  write.table(res_df,
              file.path(tbl_dir, paste0("tT_int_", label, ".tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(subset(res_df, sig == "up"),
              file.path(tbl_dir, paste0("up_int_", label, ".tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(subset(res_df, sig == "down"),
              file.path(tbl_dir, paste0("down_int_", label, ".tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  
  message("Saved tables for: ", label,
          " | up=", sum(res_df$sig == "up"),
          " | down=", sum(res_df$sig == "down"))
}

message("\n>>> Part 1 complete — Interaction DESeq2 results saved.")
message("Tables → ", tbl_dir)

## ============================================================
## PART 2 — Plotting (MA & Volcano)
## ============================================================

col_scale <- c("up" = "#D62728", "down" = "#1F77B4", "ns" = "gray80")
files <- list.files(tbl_dir, pattern = "^tT_int_.*\\.tsv$", full.names = TRUE)
message("Plotting ", length(files), " result files...")

for (fp in files) {
  nm     <- sub("^tT_int_", "", tools::file_path_sans_ext(basename(fp)))
  res_df <- read.delim(fp, check.names = FALSE)
  res_df$cat <- factor(res_df$sig, levels = c("ns", "down", "up"))
  plot_df <- res_df %>% arrange(cat)
  
  ## --- MA Plot ---
  p_ma <- ggplot(plot_df, aes(x = log10(baseMean + 1), y = log2FoldChange, color = cat)) +
    geom_point(size = 0.6, alpha = 0.85) +
    scale_color_manual(values = col_scale,
                       breaks = c("up", "down", "ns"),
                       labels = c("Up", "Down", "NS"),
                       name = "DE category") +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(title = paste("MA (interaction):", nm),
         x = "log10(baseMean + 1)", y = "log2(Fold Change)") +
    theme_bw(10)
  
  ## --- Volcano Plot ---
  p_vol <- ggplot(plot_df, aes(x = log2FoldChange,
                               y = -log10(pmax(padj, .Machine$double.xmin)),
                               color = cat)) +
    geom_point(size = 0.6, alpha = 0.8) +
    scale_color_manual(values = col_scale,
                       breaks = c("up", "down", "ns"),
                       labels = c("Up", "Down", "NS"),
                       name = "DE category") +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed") +
    labs(title = paste("Volcano (interaction):", nm),
         x = "log2(Fold Change)", y = "-log10(adjusted p-value)") +
    theme_bw(10)
  
  ## --- Save plots ---
  ggsave(file.path(fig_dir, paste0("MA_",      nm, ".pdf")),  p_ma,  width = 6,   height = 5)
  ggsave(file.path(fig_dir, paste0("MA_",      nm, ".tiff")), p_ma,  width = 6,   height = 5,
         dpi = 600, device = "tiff", compression = "lzw")
  ggsave(file.path(fig_dir, paste0("Volcano_", nm, ".pdf")),  p_vol, width = 6.5, height = 5)
  ggsave(file.path(fig_dir, paste0("Volcano_", nm, ".tiff")), p_vol, width = 6.5, height = 5,
         dpi = 600, device = "tiff", compression = "lzw")
  
  message("Saved plots for: ", nm)
}

message("\n>>> Part 2 complete — MA & Volcano plots saved.")
message("Figures saved in: ", fig_dir)