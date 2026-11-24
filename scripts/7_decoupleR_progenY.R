## ============================================================
## 4_decoupler_TF.R
## TF activity (decoupleR)
## Input: gene_counts.tsv + metadata.txt
## Design: drug × timepoint × OHT (Mus musculus)
## ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(DESeq2)
  library(decoupleR)
  library(dorothea)
  library(pheatmap)
  library(ggrepel)
})

set.seed(42)
message(">>> 4_decoupler_TF.R starting...")

## ============================================================
## 0. Directories
## ============================================================
setwd("/hpcnfs/data/BA/glucocorticoid_keep/sinaR/m.rezaei")
base_dir <- getwd()
data_dir <- file.path(base_dir, "data")
res_dir <- file.path(base_dir, "results_0.05", "4_decoupler_TF")
tbl_dir <- file.path(res_dir, "tables")
fig_dir <- file.path(res_dir, "figures")
dir.create(tbl_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

padj_cutoff <- 0.05
min_targets <- 5
drugs <- c("Budesonide", "Prednisolone")
timepoints <- c("2", "4", "8", "16")

## ============================================================
## 1. Load data
## ============================================================
message("Loading counts and metadata...")

metadata <- read.delim(file.path(base_dir, "metadata.txt"), check.names = FALSE) %>%
  mutate(
    Batch = factor(Batch),
    Timepoint = factor(Timepoint),
    OHT = factor(OHT, levels = c("OFF", "ON")),
    drug = factor(drug),
    group = factor(group)
  )
rownames(metadata) <- metadata$SampleID

counts_raw <- readr::read_tsv(file.path(base_dir, "gene_counts.tsv"), show_col_types = FALSE) %>% as.data.frame()

## Strip version from gene IDs
strip_version <- function(x) sub("\\..*$", "", x)
counts_raw[[1]] <- strip_version(counts_raw[[1]])
rownames(counts_raw) <- counts_raw[[1]]
counts_raw <- counts_raw[, -1, drop = FALSE]

## Align samples
shared_samples <- intersect(colnames(counts_raw), metadata$SampleID)
counts_raw <- counts_raw[, shared_samples, drop = FALSE]
metadata <- metadata[shared_samples, , drop = FALSE]
stopifnot(all(colnames(counts_raw) == rownames(metadata)))
message("  Samples: ", ncol(counts_raw), "   Genes: ", nrow(counts_raw))

## ============================================================
## 2. VST normalisation
## ============================================================
message("Running DESeq2 VST normalisation...")

dds <- DESeqDataSetFromMatrix(
  countData = round(counts_raw),
  colData = metadata,
  design = ~ Batch + group
)
dds <- dds[rowSums(counts(dds) >= 10) >= 4, ]
vsd <- vst(dds, blind = FALSE)

## Matrix: genes × samples
annotation <- read.delim(file.path(base_dir, "annotation.txt"), check.names = FALSE) %>%
  distinct(gene_id, .keep_all = TRUE)
id2sym <- setNames(annotation$symbol, strip_version(annotation$gene_id))

mat_vst <- assay(vsd)
rownames(mat_vst) <- ifelse(
  !is.na(id2sym[rownames(mat_vst)]) & id2sym[rownames(mat_vst)] != "",
  id2sym[rownames(mat_vst)],
  rownames(mat_vst)
)

## Remove duplicate symbols; keep highest variance
mat_vst <- mat_vst[order(-apply(mat_vst, 1, var)), , drop = FALSE]
mat_vst <- mat_vst[!duplicated(rownames(mat_vst)), , drop = FALSE]

## ============================================================
## 3. Build contrast stat matrices from shrunk DESeq2 tables
## ============================================================
message("Loading shrunk DESeq2 results for contrast matrices...")

shr_tbl <- file.path(base_dir, "results_0.05", "1_DEAnalysis", "tables", "shrunk")
tt_files <- list.files(shr_tbl, pattern = "^tT_.*\\.tsv$", full.names = TRUE)

stat_list <- lapply(tt_files, function(fp) {
  tag <- sub("^tT_", "", tools::file_path_sans_ext(basename(fp)))
  df <- tryCatch(read.delim(fp, check.names = FALSE), error = function(e) NULL)
  if (is.null(df)) return(NULL)
  if (!all(c("gene_id", "log2FoldChange", "lfcSE") %in% colnames(df))) return(NULL)
  df %>%
    transmute(
      gene_id = strip_version(gene_id),
      symbol = id2sym[gene_id],
      stat = log2FoldChange / lfcSE
    ) %>%
    filter(!is.na(stat), is.finite(stat)) %>%
    mutate(label = ifelse(!is.na(symbol) & symbol != "", symbol, gene_id)) %>%
    dplyr::select(label, stat) %>%
    distinct(label, .keep_all = TRUE) %>%
    rename(!!tag := stat)
})
stat_list <- Filter(Negate(is.null), stat_list)

stat_mat <- Reduce(function(a, b) full_join(a, b, by = "label"), stat_list)
rownames(stat_mat) <- stat_mat$label
stat_mat <- stat_mat %>% dplyr::select(-label) %>% as.matrix()

message("  Contrast stat matrix: ", nrow(stat_mat), " genes × ", ncol(stat_mat), " contrasts")

write.table(
  as.data.frame(stat_mat) %>% rownames_to_column("gene"),
  file.path(tbl_dir, "contrast_stat_matrix.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

## ============================================================
## 4. decoupleR — TF activity inference
##    Mouse DoRothEA A/B/C
##    No OmniPath connection required
## ============================================================
message("Loading mouse DoRothEA TF regulon...")

data(dorothea_mm, package = "dorothea")

net_tf <- dorothea_mm %>%
  filter(confidence %in% c("A", "B", "C")) %>%
  dplyr::select(source = tf, target, mor) %>%
  filter(!is.na(source), !is.na(target), !is.na(mor)) %>%
  distinct()

message(
  "  TF network: ", nrow(net_tf),
  " interactions, ", length(unique(net_tf$source)), " TFs"
)

## ---- 4a. TF activity on VST matrix (sample-level) ----
message("Running decoupleR ULM on VST matrix (sample-level)...")

tf_acts_sample <- decoupleR::run_ulm(
  mat = mat_vst,
  net = net_tf,
  .source = source,
  .target = target,
  .mor = mor,
  minsize = min_targets
)

tf_sample_wide <- tf_acts_sample %>%
  filter(statistic == "ulm") %>%
  dplyr::select(source, condition, score) %>%
  pivot_wider(names_from = condition, values_from = score) %>%
  column_to_rownames("source") %>%
  as.matrix()

write.table(
  as.data.frame(tf_sample_wide) %>% rownames_to_column("TF"),
  file.path(tbl_dir, "TF_activity_sample_ULM.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

## ---- 4b. TF activity on contrast stat matrix ----
message("Running decoupleR ULM on contrast stat matrix...")

tf_acts_contrast <- decoupleR::run_ulm(
  mat = stat_mat,
  net = net_tf,
  .source = source,
  .target = target,
  .mor = mor,
  minsize = min_targets
)

tf_contrast_wide <- tf_acts_contrast %>%
  filter(statistic == "ulm") %>%
  dplyr::select(source, condition, score) %>%
  pivot_wider(names_from = condition, values_from = score) %>%
  column_to_rownames("source") %>%
  as.matrix()

tf_contrast_pval <- tf_acts_contrast %>%
  filter(statistic == "ulm") %>%
  dplyr::select(source, condition, p_value) %>%
  pivot_wider(names_from = condition, values_from = p_value) %>%
  column_to_rownames("source") %>%
  as.matrix()

write.table(
  as.data.frame(tf_contrast_wide) %>% rownames_to_column("TF"),
  file.path(tbl_dir, "TF_activity_contrast_ULM.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

write.table(
  as.data.frame(tf_contrast_pval) %>% rownames_to_column("TF"),
  file.path(tbl_dir, "TF_activity_contrast_pval.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

## ============================================================
## 5. Plotting helpers
## ============================================================

make_star_mat <- function(pval_mat, cutoff = 0.05) {
  star <- matrix("", nrow = nrow(pval_mat), ncol = ncol(pval_mat), dimnames = dimnames(pval_mat))
  star[!is.na(pval_mat) & pval_mat < cutoff] <- "★"
  star
}

save_heatmap_ph <- function(mat, star_mat = NULL, title, out_base,
                            cluster_rows = TRUE, cluster_cols = FALSE,
                            row_fontsize = 8, col_fontsize = 8,
                            annotation_col = NA, annotation_colors = NA) {
  max_abs <- max(abs(mat), na.rm = TRUE)
  if (!is.finite(max_abs) || max_abs == 0) max_abs <- 1
  breaks <- seq(-max_abs, max_abs, length.out = 101)
  pal <- colorRampPalette(c("#2166AC", "white", "#B2182B"))(100)
  n_rows <- nrow(mat); n_cols <- ncol(mat)
  fig_h <- max(4, n_rows * 0.32 + 2.5)
  fig_w <- max(6, n_cols * 0.5 + 4)
  display_num <- if (!is.null(star_mat)) star_mat else FALSE
  ph <- pheatmap::pheatmap(
    mat, color = pal, breaks = breaks, main = title,
    annotation_col = annotation_col, annotation_colors = annotation_colors,
    cluster_rows = cluster_rows && n_rows >= 2, cluster_cols = cluster_cols,
    clustering_method = "ward.D2", treeheight_row = 25, treeheight_col = 0,
    fontsize_row = row_fontsize, fontsize_col = col_fontsize, fontsize = 8,
    angle_col = 45, na_col = "grey92", display_numbers = display_num,
    number_color = "black", fontsize_number = 8,
    legend_breaks = c(-max_abs, -max_abs/2, 0, max_abs/2, max_abs),
    legend_labels = round(c(-max_abs, -max_abs/2, 0, max_abs/2, max_abs), 2),
    border_color = NA, cellwidth = 24, cellheight = 16, silent = TRUE
  )
  ggsave(paste0(out_base, ".pdf"), ph, width = fig_w, height = fig_h, limitsize = FALSE)
  ggsave(paste0(out_base, ".tiff"), ph, width = fig_w, height = fig_h,
         limitsize = FALSE, dpi = 300, device = "tiff", compression = "lzw")
  invisible(ph)
}

## ============================================================
## 6. TF activity heatmaps
## ============================================================
message(">>> Plotting TF activity...")

## ---- 6a. Sample-level heatmap: top variable TFs ----
ann_col_sample <- metadata %>%
  dplyr::select(drug, Timepoint, OHT) %>%
  mutate(Timepoint = as.character(Timepoint), OHT = as.character(OHT)) %>%
  as.data.frame()
rownames(ann_col_sample) <- metadata$SampleID

ann_colors_sample <- list(
  drug = c(Budesonide = "#D62728", Prednisolone = "#1F77B4", DMSO = "#7F7F7F"),
  OHT = c(OFF = "#D1E5F0", ON = "#2166AC"),
  Timepoint = c("0" = "#F7F7F7", "2" = "#D1E5F0", "4" = "#92C5DE",
                "8" = "#4393C3", "16" = "#2166AC")
)

top_tfs_var <- names(sort(apply(tf_sample_wide, 1, var, na.rm = TRUE),
                          decreasing = TRUE))[1:min(50, nrow(tf_sample_wide))]
mat_tf_top <- tf_sample_wide[top_tfs_var, rownames(ann_col_sample), drop = FALSE]

save_heatmap_ph(
  mat = mat_tf_top,
  title = "TF activity (ULM) — top 50 variable TFs",
  out_base = file.path(fig_dir, "TF_sample_heatmap_top50"),
  cluster_cols = TRUE,
  annotation_col = ann_col_sample,
  annotation_colors = ann_colors_sample
)
message("  Saved: TF_sample_heatmap_top50")

## ---- 6b. Contrast-level HM2 heatmaps ----
get_contrast_cols <- function(score_mat, pval_mat, col_tags, col_labels,
                              sig_tag, cutoff = 0.05) {
  
  ## Check that the requested significance contrast exists
  if (!sig_tag %in% colnames(pval_mat)) {
    message("  Significance contrast not found: ", sig_tag)
    return(NULL)
  }
  
  ## Keep only contrasts that are present
  present <- col_tags[col_tags %in% colnames(score_mat)]
  labels <- col_labels[col_tags %in% colnames(score_mat)]
  
  if (length(present) == 0) return(NULL)
  
  ## Build heatmap matrices
  mat <- score_mat[, present, drop = FALSE]
  pmat <- pval_mat[, present, drop = FALSE]
  
  colnames(mat) <- colnames(pmat) <- labels
  
  ## ------------------------------------------------------------
  ## IMPORTANT:
  ## Keep ONLY TFs significant at 16h Drug + OHT vs OHT
  ## ------------------------------------------------------------
  sig_rows <- !is.na(pval_mat[, sig_tag]) &
    pval_mat[, sig_tag] < cutoff
  
  mat <- mat[sig_rows, , drop = FALSE]
  pmat <- pmat[sig_rows, , drop = FALSE]
  
  if (nrow(mat) == 0) return(NULL)
  
  ## Rank retained TFs by strongest activity across all displayed
  ## contrasts
  top_n <- min(40, nrow(mat))
  
  top_idx <- order(
    apply(abs(mat), 1, max, na.rm = TRUE),
    decreasing = TRUE
  )[seq_len(top_n)]
  
  list(
    mat = mat[top_idx, , drop = FALSE],
    star = make_star_mat(
      pmat[top_idx, , drop = FALSE],
      cutoff = cutoff
    )
  )
}

for (drug in drugs) {
  
  message("  TF HM2 heatmap — ", drug)
  
  hm2_tags <- c(
    paste0(drug, "_", timepoints),
    paste0(drug, "_", timepoints, "_OHT")
  )
  
  hm2_labels <- c(
    paste0(drug, " T", timepoints, " vs DMSO"),
    paste0(drug, "+OHT T", timepoints, " vs OHT")
  )
  
  ## Significance is determined ONLY from:
  ## 16h Drug + OHT vs OHT
  sig_tag <- paste0(drug, "_16_OHT")
  
  res <- get_contrast_cols(
    score_mat = tf_contrast_wide,
    pval_mat = tf_contrast_pval,
    col_tags = hm2_tags,
    col_labels = hm2_labels,
    sig_tag = sig_tag,
    cutoff = padj_cutoff
  )
  
  if (is.null(res)) {
    message("    No TFs significant at 16h Drug + OHT vs OHT, skipped")
    next
  }
  
  out_base <- file.path(
    fig_dir,
    paste0("TF_", drug, "_HM2_sig16hOHT")
  )
  
  save_heatmap_ph(
    mat = res$mat,
    star_mat = res$star,
    title = paste0(
      "TF activity — ", drug,
      " — TFs significant at 16h Drug + OHT vs OHT"
    ),
    out_base = out_base
  )
  
  write.table(
    as.data.frame(res$mat) %>%
      rownames_to_column("TF"),
    file.path(
      tbl_dir,
      paste0("TF_", drug, "_HM2_sig16hOHT_scores.tsv")
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}

## ============================================================
## 7. TF activity scatterplots
## X = Drug vs Untreated
## Y = Drug + OHT vs OHT
## ============================================================
message(">>> Plotting TF scatterplots...")

for (drug in drugs) {
  scatter_list <- lapply(timepoints, function(tp) {
    x_tag <- paste0(drug, "_", tp)
    y_tag <- paste0(drug, "_", tp, "_OHT")
    if (!x_tag %in% colnames(tf_contrast_wide) || !y_tag %in% colnames(tf_contrast_wide)) {
      message("    Missing contrasts for ", drug, " ", tp, "h")
      return(NULL)
    }
    data.frame(
      TF = rownames(tf_contrast_wide),
      x_score = tf_contrast_wide[, x_tag],
      y_score = tf_contrast_wide[, y_tag],
      x_pval = tf_contrast_pval[, x_tag],
      y_pval = tf_contrast_pval[, y_tag],
      Timepoint = tp,
      stringsAsFactors = FALSE
    ) %>%
      filter(is.finite(x_score), is.finite(y_score)) %>%
      mutate(
        sig_x = !is.na(x_pval) & x_pval < padj_cutoff,
        sig_y = !is.na(y_pval) & y_pval < padj_cutoff,
        sig_group = case_when(
          sig_x & sig_y ~ "Significant in both",
          sig_x & !sig_y ~ "Drug vs Untreated only",
          !sig_x & sig_y ~ "Drug + OHT vs OHT only",
          TRUE ~ "Not significant"
        ),
        importance = abs(x_score) + abs(y_score)
      )
  }) %>% bind_rows()
  
  if (nrow(scatter_list) == 0) next
  
  scatter_list$Timepoint <- factor(
    scatter_list$Timepoint,
    levels = timepoints,
    labels = paste0(timepoints, " h")
  )
  
  ## Label only the strongest significant TFs in each timepoint
  ## Label TFs with priority on Drug + OHT vs OHT response
  label_df <- scatter_list %>%
    filter(sig_group != "Not significant") %>%
    mutate(
      label_priority = case_when(
        sig_group == "Drug + OHT vs OHT only" ~ 1,
        sig_group == "Significant in both" ~ 2,
        sig_group == "Drug vs Untreated only" ~ 3
      ),
      label_strength = abs(y_score)
    ) %>%
    group_by(Timepoint) %>%
    arrange(label_priority, desc(label_strength)) %>%
    slice_head(n = 8) %>%
    ungroup()
  
  p <- ggplot(scatter_list, aes(x = x_score, y = y_score)) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3, color = "grey65") +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3, color = "grey65") +
    geom_point(aes(color = sig_group), size = 1.8, alpha = 0.75) +
    geom_text_repel(
      data = label_df,
      aes(label = TF, color = sig_group),
      size = 2.7,
      fontface = "plain",
      box.padding = 0.5,
      point.padding = 0.25,
      force = 1.5,
      max.overlaps = Inf,
      min.segment.length = 0,
      segment.size = 0.25,
      show.legend = FALSE
    ) +
    facet_wrap(~Timepoint, nrow = 1) +
    scale_color_manual(values = c(
      "Not significant" = "grey78",
      "Drug vs Untreated only" = "#2166AC",
      "Drug + OHT vs OHT only" = "#D62728",
      "Significant in both" = "#7B3294"
    )) +
    labs(
      title = paste0("TF activity — ", drug),
      x = paste0(drug, " vs Untreated — ULM score"),
      y = paste0(drug, " + OHT vs OHT — ULM score"),
      color = NULL
    ) +
    theme_bw(10) +
    theme(
      panel.grid = element_blank(),
      strip.background = element_rect(fill = "grey95", color = "grey70"),
      strip.text = element_text(face = "bold"),
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.position = "bottom",
      legend.text = element_text(size = 9)
    )
  
  print(p)
  
  out_base <- file.path(fig_dir, paste0("TF_scatter_HM2_", drug))
  
  ggsave(paste0(out_base, ".pdf"), p, width = 15, height = 5, limitsize = FALSE)
  ggsave(paste0(out_base, ".tiff"), p, width = 15, height = 5, dpi = 300,
         device = "tiff", compression = "lzw", limitsize = FALSE)
  
  write.table(scatter_list, paste0(out_base, "_data.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  
  message("  Saved: TF_scatter_HM2_", drug)
}


## ---- Separate tables for the three significance categories ----
sig_categories <- c("Drug vs Untreated only", "Drug + OHT vs OHT only", "Significant in both")

for (cat in sig_categories) {
  df_cat <- scatter_list %>% dplyr::filter(sig_group == cat)
  if (nrow(df_cat) == 0) next
  
  cat_tag <- gsub("[^A-Za-z0-9]+", "_", cat)  # safe filename
  write.table(
    df_cat,
    file.path(tbl_dir, paste0("TF_scatter_HM2_", drug, "_", cat_tag, ".tsv")),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  message("    Saved: TF_scatter_HM2_", drug, "_", cat_tag, ".tsv (", nrow(df_cat), " TF×timepoint rows)")
}

## ---- 4c. Nfe2l2 target genes driving TF activity ----
message("Extracting genes associated with Nfe2l2 activity...")

target_tf <- "Nrf1"

## Get Nfe2l2 regulon
nfe2l2_net <- net_tf %>%
  filter(source == target_tf) %>%
  distinct(target, mor)

message("  Klf9 targets in DoRothEA: ", nrow(nfe2l2_net))

## Keep only genes present in the contrast matrix
nfe2l2_targets <- nfe2l2_net %>%
  filter(target %in% rownames(stat_mat))

## Build gene-level table for every contrast
nfe2l2_gene_table <- lapply(colnames(stat_mat), function(con) {

  data.frame(
    TF = target_tf,
    gene = nfe2l2_targets$target,
    mor = nfe2l2_targets$mor,
    stat = stat_mat[nfe2l2_targets$target, con],
    contrast = con,
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      gene_effect = stat * mor,
      predicted_direction = case_when(
        gene_effect > 0 ~ "Activates Nrf1",
        gene_effect < 0 ~ "Inhibits Nrf1",
        TRUE ~ "No effect"
      )
    )
}) %>%
  bind_rows()

## Save complete Nfe2l2 target table
write.table(
  nfe2l2_gene_table,
  file.path(tbl_dir, "Nrf1_target_genes_all_contrasts.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

## ------------------------------------------------------------
## Significant genes driving Nfe2l2 activity
## ------------------------------------------------------------

## Get gene-level DESeq2 padj for each contrast
gene_padj_list <- lapply(tt_files, function(fp) {

  tag <- sub("^tT_", "", tools::file_path_sans_ext(basename(fp)))

  df <- tryCatch(
    read.delim(fp, check.names = FALSE),
    error = function(e) NULL
  )

  if (is.null(df)) return(NULL)

  df %>%
    transmute(
      gene_id = strip_version(gene_id),
      symbol = id2sym[gene_id],
      padj = padj
    ) %>%
    mutate(
      gene = ifelse(
        !is.na(symbol) & symbol != "",
        symbol,
        gene_id
      ),
      contrast = tag
    ) %>%
    dplyr::select(gene, contrast, padj) %>%
    distinct(gene, contrast, .keep_all = TRUE)
}) %>%
  Filter(Negate(is.null), .) %>%
  bind_rows()

## Add gene-level significance
nfe2l2_gene_table_sig <- nfe2l2_gene_table %>%
  left_join(
    gene_padj_list,
    by = c("gene", "contrast")
  ) %>%
  mutate(
    significant_gene = !is.na(padj) & padj < padj_cutoff
  )

## Save all targets with gene-level significance
write.table(
  nfe2l2_gene_table_sig,
  file.path(tbl_dir, "Nrf1_target_genes_with_DEG_significance.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

## ------------------------------------------------------------
## Only significant genes that activate/inhibit Nfe2l2
## ------------------------------------------------------------

nfe2l2_sig_drivers <- nfe2l2_gene_table_sig %>%
  filter(significant_gene) %>%
  arrange(contrast, desc(abs(gene_effect)))

write.table(
  nfe2l2_sig_drivers,
  file.path(tbl_dir, "Nfe2l2_significant_driver_genes.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message(
  "  Significant Nfe2l2 driver genes: ",
  nrow(nfe2l2_sig_drivers)
)
## ============================================================
## 8. Reproducibility
## ============================================================

writeLines(
  capture.output(sessionInfo()),
  file.path(res_dir, "sessionInfo.txt")
)

message(">>> All done.")
message("Tables : ", tbl_dir)
message("Figures: ", fig_dir)

## ============================================================
## 8. PROGENy pathway activity plots
## ============================================================
message(">>> Plotting PROGENy pathway activity...")

## ---- 8a. Sample-level heatmap (all pathways) ----
mat_prog_sample <- prog_sample_wide[, rownames(ann_col_sample), drop = FALSE]

save_heatmap_ph(
  mat            = mat_prog_sample,
  title          = "PROGENy pathway activity (MLM) — all samples",
  out_base       = file.path(fig_dir, "PROGENy_sample_heatmap"),
  cluster_cols   = TRUE,
  annotation_col = ann_col_sample,
  annotation_colors = ann_colors_sample,
  row_fontsize   = 9
)
message("  Saved: PROGENy_sample_heatmap")

## ---- 8b. Contrast-level heatmaps: per drug, HM1/2/3 ----
for (drug in drugs) {
  message("  PROGENy heatmaps — ", drug)
  
  hm1_tags   <- c("DMSO_OHT_vs_DMSO",
                  paste0(drug, "_", timepoints),
                  paste0(drug, "_", timepoints, "_OHT_vs_DMSO"))
  hm1_labels <- c("OHT vs DMSO",
                  paste0(drug, " T", timepoints, " vs DMSO"),
                  paste0(drug, "+OHT T", timepoints, " vs DMSO"))
  
  hm2_tags   <- c(paste0(drug, "_", timepoints),
                  paste0(drug, "_", timepoints, "_OHT"))
  hm2_labels <- c(paste0(drug, " T", timepoints, " vs DMSO"),
                  paste0(drug, "+OHT T", timepoints, " vs OHT"))
  
  hm3_tags   <- c("DMSO_OHT_vs_DMSO",
                  paste0(drug, "_", timepoints, "_OHT__", drug, "_", timepoints))
  hm3_labels <- c("OHT vs DMSO",
                  paste0(drug, "+OHT T", timepoints, " vs ", drug))
  
  for (hm in list(
    list(tags = hm1_tags, labels = hm1_labels, name = "HM1"),
    list(tags = hm2_tags, labels = hm2_labels, name = "HM2"),
    list(tags = hm3_tags, labels = hm3_labels, name = "HM3")
  )) {
    ## For PROGENy all pathways are shown (only ~14), no significance filter needed
    present <- hm$tags[hm$tags %in% colnames(prog_contrast_wide)]
    labels  <- hm$labels[hm$tags %in% colnames(prog_contrast_wide)]
    if (length(present) == 0) next
    
    mat_hm  <- prog_contrast_wide[, present, drop = FALSE]
    pmat_hm <- prog_contrast_pval[, present, drop = FALSE]
    colnames(mat_hm) <- colnames(pmat_hm) <- labels
    
    out_base <- file.path(fig_dir, paste0("PROGENy_", drug, "_", hm$name))
    save_heatmap_ph(
      mat      = mat_hm,
      star_mat = make_star_mat(pmat_hm),
      title    = paste0("PROGENy — ", drug, " — ", hm$name),
      out_base = out_base,
      row_fontsize = 9
    )
    write.table(
      as.data.frame(mat_hm) %>% rownames_to_column("pathway"),
      file.path(tbl_dir, paste0("PROGENy_", drug, "_", hm$name, "_scores.tsv")),
      sep = "\t", quote = FALSE, row.names = FALSE
    )
    message("    Saved: PROGENy_", drug, "_", hm$name)
  }
}

## ---- 8c. PROGENy bar plots per contrast ----
message("  PROGENy bar plots per contrast...")

prog_contrast_long <- prog_acts_contrast %>%
  filter(statistic == "mlm") %>%
  mutate(
    neg_log10p  = -log10(p_value + 1e-300),
    significant = p_value < padj_cutoff,
    dir         = ifelse(score > 0, "Activated", "Repressed")
  )

for (contrast_tag in unique(prog_contrast_long$condition)) {
  df_c <- prog_contrast_long %>%
    filter(condition == contrast_tag) %>%
    mutate(source = fct_reorder(source, score))
  
  p <- ggplot(df_c, aes(x = score, y = source, fill = dir)) +
    geom_col() +
    geom_point(
      data = dplyr::filter(df_c, significant),
      aes(x = score + sign(score) * 0.05, y = source),
      shape = 8, size = 2, color = "black",
      inherit.aes = FALSE
    ) +
    geom_text(aes(x = 0, label = ifelse(significant,
                                        formatC(p_value, digits = 2, format = "e"),
                                        "")),
              size = 2.6, color = "black", hjust = 0.5) +
    scale_fill_manual(values = c(Activated = "#D62728", Repressed = "#1F77B4")) +
    labs(title = paste0("PROGENy: ", contrast_tag),
         x = "MLM score", y = NULL, fill = NULL) +
    theme_bw(9) +
    theme(panel.grid = element_blank(),
          plot.title = element_text(face = "bold", hjust = 0.5, size = 9))
  
  out <- file.path(fig_dir, paste0("PROGENy_bar_", contrast_tag))
  ggsave(paste0(out, ".pdf"),  p, width = 5, height = 4)
  ggsave(paste0(out, ".tiff"), p, width = 5, height = 4,
         dpi = 300, device = "tiff", compression = "lzw")
}

## ============================================================
## 9. Joint TF + PROGENy summary: apoptosis-relevant regulators
##    Focus: TP53, E2F1, MYC, BCL2, CASP3-adjacent TFs
##    + p53, TNF, Apoptosis, PI3K PROGENy pathways
## ============================================================
message(">>> Apoptosis-focused summary plots...")

apoptosis_tfs <- c("TP53", "E2F1", "E2F2", "E2F3", "MYC", "MYCN",
                   "FOXO3", "FOXO1", "NFκB1", "NFKB1", "RELA",
                   "BCL2", "BAX", "SP1", "YBX1", "HIF1A",
                   "STAT3", "STAT1", "IRF1", "NR3C1")

## Filter to TFs actually present in our results
apop_tfs_present <- intersect(apoptosis_tfs, rownames(tf_contrast_wide))

if (length(apop_tfs_present) > 0) {
  ## All contrasts for apoptosis TFs
  mat_apop_tf <- tf_contrast_wide[apop_tfs_present, , drop = FALSE]
  pmat_apop_tf <- tf_contrast_pval[apop_tfs_present, , drop = FALSE]
  
  ## Save full table
  write.table(
    as.data.frame(mat_apop_tf) %>% rownames_to_column("TF"),
    file.path(tbl_dir, "TF_apoptosis_relevant_scores.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  
  ## Per drug structured heatmap
  for (drug in drugs) {
    all_tags   <- c("DMSO_OHT_vs_DMSO",
                    paste0(drug, "_", timepoints),
                    paste0(drug, "_", timepoints, "_OHT_vs_DMSO"),
                    paste0(drug, "_", timepoints, "_OHT"),
                    paste0(drug, "_", timepoints, "_OHT__", drug, "_", timepoints))
    all_labels <- c("OHT vs DMSO",
                    paste0(drug, " T", timepoints, " vs DMSO"),
                    paste0(drug, "+OHT T", timepoints, " vs DMSO"),
                    paste0(drug, "+OHT T", timepoints, " vs OHT"),
                    paste0(drug, "+OHT T", timepoints, " vs ", drug))
    
    present <- all_tags[all_tags %in% colnames(mat_apop_tf)]
    labels  <- all_labels[all_tags %in% colnames(mat_apop_tf)]
    if (length(present) == 0) next
    
    mat_sub  <- mat_apop_tf[, present, drop = FALSE]
    pmat_sub <- pmat_apop_tf[, present, drop = FALSE]
    colnames(mat_sub) <- colnames(pmat_sub) <- labels
    
    out_base <- file.path(fig_dir, paste0("TF_apoptosis_", drug, "_allcontrasts"))
    save_heatmap_ph(
      mat      = mat_sub,
      star_mat = make_star_mat(pmat_sub),
      title    = paste0("Apoptosis-relevant TF activity — ", drug),
      out_base = out_base,
      row_fontsize = 9
    )
    message("  Saved: TF_apoptosis_", drug, "_allcontrasts")
  }
}

## PROGENy apoptosis-adjacent pathways
apop_prog_pathways <- c("p53", "TNFa", "Apoptosis", "PI3K", "MAPK",
                        "NFkB", "Hypoxia", "Trail")
apop_prog_present  <- intersect(apop_prog_pathways, rownames(prog_contrast_wide))

if (length(apop_prog_present) > 0) {
  mat_apop_prog  <- prog_contrast_wide[apop_prog_present, , drop = FALSE]
  pmat_apop_prog <- prog_contrast_pval[apop_prog_present, , drop = FALSE]
  
  write.table(
    as.data.frame(mat_apop_prog) %>% rownames_to_column("pathway"),
    file.path(tbl_dir, "PROGENy_apoptosis_relevant_scores.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  
  for (drug in drugs) {
    all_tags   <- c("DMSO_OHT_vs_DMSO",
                    paste0(drug, "_", timepoints),
                    paste0(drug, "_", timepoints, "_OHT_vs_DMSO"),
                    paste0(drug, "_", timepoints, "_OHT"),
                    paste0(drug, "_", timepoints, "_OHT__", drug, "_", timepoints))
    all_labels <- c("OHT vs DMSO",
                    paste0(drug, " T", timepoints, " vs DMSO"),
                    paste0(drug, "+OHT T", timepoints, " vs DMSO"),
                    paste0(drug, "+OHT T", timepoints, " vs OHT"),
                    paste0(drug, "+OHT T", timepoints, " vs ", drug))
    
    present <- all_tags[all_tags %in% colnames(mat_apop_prog)]
    labels  <- all_labels[all_tags %in% colnames(mat_apop_prog)]
    if (length(present) == 0) next
    
    mat_sub  <- mat_apop_prog[, present, drop = FALSE]
    pmat_sub <- pmat_apop_prog[, present, drop = FALSE]
    colnames(mat_sub) <- colnames(pmat_sub) <- labels
    
    out_base <- file.path(fig_dir, paste0("PROGENy_apoptosis_", drug, "_allcontrasts"))
    save_heatmap_ph(
      mat      = mat_sub,
      star_mat = make_star_mat(pmat_sub),
      title    = paste0("Apoptosis-relevant PROGENy — ", drug),
      out_base = out_base,
      row_fontsize = 9
    )
    message("  Saved: PROGENy_apoptosis_", drug, "_allcontrasts")
  }
}

## ============================================================
## 10. Reproducibility
## ============================================================
writeLines(capture.output(sessionInfo()),
           file.path(res_dir, "sessionInfo.txt"))

message(">>> All done.")
message("Tables : ", tbl_dir)
message("Figures: ", fig_dir)













# Instead of looking at contrasts, score ant TF activity 
# per sample from the VST matrix and plot over time

nr3c1_sample_scores <- tf_acts_sample %>%
  filter(statistic == "ulm", source == "Myc") %>%
  left_join(metadata %>% dplyr::select(SampleID, drug, Timepoint, OHT, group),
            by = c("condition" = "SampleID"))

p <- ggplot(nr3c1_sample_scores,
            aes(x = factor(Timepoint), y = score,
                color = OHT, group = interaction(drug, OHT))) +
  geom_point(size = 2, position = position_dodge(0.3)) +
  geom_line(stat = "summary", fun = mean,
            position = position_dodge(0.3)) +
  facet_wrap(~ drug) +
  scale_color_manual(values = c(OFF = "#1F77B4", ON = "#D62728")) +
  labs(title = "Foxo1 activity over time (absolute, per sample)",
       x = "Timepoint (h)", y = "ULM score") +
  theme_bw(10)
print(p)








library(decoupleR)

url <- "https://rescued.omnipathdb.org/CollecTRI2.tsv.gz"

download.file(
  url,
  destfile = "/hpcnfs/data/BA/glucocorticoid_keep/sinaR/m.rezaei/bax-bac-KO/results/5_decoupler_TF/CollecTRI2.tsv.gz",
  mode = "wb"
)

write.table(
  TFs,
  "CollecTRI_mouse_TFs.txt",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)
