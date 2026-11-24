## ============================================================
## GR Signature Enrichment Analysis
## Fisher's exact test: independent/additive categories vs
## PanCancer GR signature + Choi et al. GR signature
## ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(homologene)
  library(ggplot2)
})

strip_version <- function(x) sub("\\..*$", "", x)

## -------------------- Paths --------------------
base_root <- "/hpcnfs/data/BA/glucocorticoid_keep/sinaR/m.rezaei"
cat_tbl   <- file.path(base_root, "results/4_categories/tables")
fig_out   <- file.path(base_root, "results/8_GR_signature_enrichment/figures")
tbl_out   <- file.path(base_root, "results/8_GR_signature_enrichment/tables")
dir.create(fig_out, recursive = TRUE, showWarnings = FALSE)
dir.create(tbl_out, recursive = TRUE, showWarnings = FALSE)

## -------------------- Load annotation (for symbol mapping) --------------------
annotation <- read.delim(
  file.path(base_root, "annotation.txt"),
  check.names = FALSE
)
sym_col  <- if ("symbol" %in% names(annotation)) "symbol" else "gene_name"
gene2sym <- annotation %>%
  dplyr::select(gene_id, SYMBOL = all_of(sym_col)) %>%
  mutate(gene_id = strip_version(gene_id)) %>%
  distinct(gene_id, .keep_all = TRUE)

## -------------------- Load universe --------------------
## All genes tested in DE analysis
universe_raw <- read_tsv(
  file.path(base_root, "results/1_DEAnalysis/tables/expressed_universe_gene.tsv"),
  show_col_types = FALSE
)
## Detect gene_id column
universe_col <- intersect(
  c("gene_id", "GeneID", "ensembl_gene_id", names(universe_raw)[1]),
  names(universe_raw)
)[1]
universe_ids <- strip_version(universe_raw[[universe_col]])
message("Universe size: ", length(universe_ids))

## ============================================================
## Load GR Signatures
## ============================================================

## ---- 1. PanCancer GR Signature (human symbols -> convert to mouse) ----
pancancer_raw <- read_xlsx(
  "/hpcnfs/data/BA/glucocorticoid_keep/sinaR/44321_2023_BFEMMM202317737_MOESM3_ESM.xlsx",
  col_names = TRUE
)

## First row is the header label "GR Activity Genes"; actual genes start row 2
pancancer_human_symbols <- pancancer_raw[[1]] %>%
  .[-1] %>%
  na.omit() %>%
  as.character() %>%
  trimws()

message("PanCancer human symbols loaded: ", length(pancancer_human_symbols))

## Convert human symbols -> mouse via homologene
## NCBI taxon IDs: 9606 = human, 10090 = mouse
pancancer_mouse_df <- homologene(
  pancancer_human_symbols,
  inTax  = 9606,
  outTax = 10090
)
## Explicitly select "10090" column = mouse gene symbols
## (last column is 10090_ID = Entrez IDs, not symbols)
pancancer_mouse_symbols <- pancancer_mouse_df[["10090"]] %>%
  na.omit() %>%
  unique()

message("PanCancer mouse orthologs found: ", length(pancancer_mouse_symbols),
        " / ", length(pancancer_human_symbols))

## Map mouse symbols -> Ensembl IDs using annotation file, then intersect universe
pancancer_mouse_ids <- gene2sym %>%
  filter(toupper(SYMBOL) %in% toupper(pancancer_mouse_symbols)) %>%
  pull(gene_id) %>%
  intersect(universe_ids)

message("PanCancer Ensembl IDs matched in universe: ", length(pancancer_mouse_ids))

## ---- 2. Choi et al. GR Signature (already mouse Ensembl IDs, use all 46) ----
choi_raw <- read_xlsx(
  "/hpcnfs/data/BA/glucocorticoid_keep/sinaR/thesis/choi_2024.xlsx"
)

choi_ids <- strip_version(choi_raw$ID) %>%
  na.omit() %>%
  unique() %>%
  intersect(universe_ids)

message("Choi et al. Ensembl IDs matched in universe: ", length(choi_ids))

## Bundle signatures
signatures <- list(
  PanCancer_GR = pancancer_mouse_ids,
  Choi_GR      = choi_ids
)

## ============================================================
## Categories to test
## ============================================================
cat_names_test <- c("independent_up", "independent_down",
                    "additive_up",    "additive_down")

drugs <- c(
  "Budesonide_2",   "Budesonide_4",   "Budesonide_8",   "Budesonide_16",
  "Prednisolone_2", "Prednisolone_4", "Prednisolone_8", "Prednisolone_16"
)

## ============================================================
## Fisher's Exact Test — per drug x category x signature
## ============================================================
N <- length(universe_ids)   # background universe size

results_list <- list()

for (drug in drugs) {
  for (cat in cat_names_test) {
    
    fp <- file.path(cat_tbl, paste0(cat, "_", drug, ".tsv"))
    if (!file.exists(fp)) next
    
    df <- read.delim(fp, check.names = FALSE)
    if (!"gene_id" %in% names(df)) next
    
    cat_genes <- strip_version(df$gene_id) %>%
      intersect(universe_ids)
    
    if (length(cat_genes) < 2) next
    
    for (sig_name in names(signatures)) {
      sig_genes <- signatures[[sig_name]]
      
      ## 2x2 contingency table:
      ##                In_Signature   Not_In_Signature
      ## In_Category         a               b
      ## Not_Category        c               d
      
      a <- length(intersect(cat_genes, sig_genes))
      b <- length(cat_genes)  - a
      c <- length(sig_genes)  - a
      d <- N - a - b - c
      if (d < 0) d <- 0
      
      mat <- matrix(c(a, b, c, d), nrow = 2,
                    dimnames = list(
                      c("In_Category",  "Not_Category"),
                      c("In_Signature", "Not_Signature")
                    ))
      
      ft <- fisher.test(mat, alternative = "greater")
      
      overlap_pct <- if (length(cat_genes) > 0)
        round(a / length(cat_genes) * 100, 1) else 0
      
      overlap_syms <- gene2sym %>%
        filter(gene_id %in% intersect(cat_genes, sig_genes)) %>%
        pull(SYMBOL) %>%
        paste(collapse = ", ")
      
      results_list[[length(results_list) + 1]] <- data.frame(
        Drug             = drug,
        Category         = cat,
        Signature        = sig_name,
        N_universe       = N,
        N_signature      = length(sig_genes),
        N_category       = length(cat_genes),
        N_overlap        = a,
        Overlap_pct      = overlap_pct,
        OddsRatio        = ft$estimate,
        P_value          = ft$p.value,
        Overlap_symbols  = overlap_syms,
        stringsAsFactors = FALSE
      )
    }
  }
}

results <- bind_rows(results_list) %>%
  group_by(Signature) %>%
  mutate(P_adj = p.adjust(P_value, method = "BH")) %>%
  ungroup()

write.table(results,
            file.path(tbl_out, "GR_signature_enrichment_results.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

message("Enrichment results saved: ", nrow(results), " tests")

## ============================================================
## Dot Plot — one per signature + one combined
## ============================================================

cat_colors_plot <- c(
  independent_up   = "#F7B6D2",
  independent_down = "#AEC7E8",
  additive_up      = "#74C476",
  additive_down    = "#FDAE6B"
)

plot_df <- results %>%
  mutate(
    Drug        = factor(Drug, levels = drugs),
    Category    = factor(Category, levels = cat_names_test),
    Sig_label   = ifelse(P_adj < 0.05, "*", ""),
    Significant = P_adj < 0.05,
    Count_label = paste0(N_overlap, "/", N_category)
  )

## ---- Per-signature plots ----
for (sig in names(signatures)) {
  
  df_sig <- plot_df %>% filter(Signature == sig)
  
  p <- ggplot(df_sig, aes(
    x     = Drug,
    y     = Category,
    size  = Overlap_pct,
    fill  = Category,
    color = Significant
  )) +
    geom_point(shape = 21, stroke = 1.4) +
    geom_text(aes(label = Sig_label), size = 5, vjust = -1.2, color = "black") +
    scale_size_continuous(
      name   = "Overlap (%)",
      range  = c(2, 14),
      breaks = c(0, 5, 10, 20, 40)
    ) +
    scale_fill_manual(values = cat_colors_plot, name = "Category") +
    scale_color_manual(
      values = c("TRUE" = "#B2182B", "FALSE" = "grey55"),
      name   = "FDR < 0.05",
      labels = c("TRUE" = "Yes", "FALSE" = "No")
    ) +
    theme_bw(base_size = 12) +
    theme(
      axis.text.x      = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y      = element_text(size = 10),
      strip.background = element_rect(fill = "grey90"),
      strip.text       = element_text(face = "bold"),
      legend.position  = "right",
      panel.grid.major = element_line(color = "grey92")
    ) +
    labs(
      title    = paste0("GR Signature Enrichment — ", gsub("_", " ", sig)),
      subtitle = "Dot size = % of category genes overlapping signature | Red border = FDR < 0.05",
      x        = NULL,
      y        = NULL
    )
  
  ggsave(file.path(fig_out, paste0("enrichment_dotplot_", sig, ".pdf")),
         p, width = 14, height = 5, device = "pdf")
  ggsave(file.path(fig_out, paste0("enrichment_dotplot_", sig, ".tiff")),
         p, width = 14, height = 5, device = "tiff",
         dpi = 300, compression = "lzw")
  
  message("  Plot saved: ", sig)
}

## ---- Combined plot (both signatures, faceted by signature row) ----
p_combined <- ggplot(plot_df, aes(
  x     = Drug,
  y     = Category,
  size  = Overlap_pct,
  fill  = Category,
  color = Significant
)) +
  geom_point(shape = 21, stroke = 1.4) +
  geom_text(aes(label = Sig_label), size = 4, vjust = -1.2, color = "black") +
  scale_size_continuous(
    name   = "Overlap (%)",
    range  = c(2, 14),
    breaks = c(0, 5, 10, 20, 40)
  ) +
  scale_fill_manual(values = cat_colors_plot, name = "Category") +
  scale_color_manual(
    values = c("TRUE" = "#B2182B", "FALSE" = "grey55"),
    name   = "FDR < 0.05",
    labels = c("TRUE" = "Yes", "FALSE" = "No")
  ) +
  facet_grid(Signature ~ ., scales = "free_y") +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 9),
    axis.text.y      = element_text(size = 9),
    strip.background = element_rect(fill = "grey90"),
    strip.text       = element_text(face = "bold", size = 10),
    legend.position  = "right",
    panel.grid.major = element_line(color = "grey92")
  ) +
  labs(
    title    = "GR Signature Enrichment — PanCancer & Choi et al.",
    subtitle = "Dot size = % overlap | Red border = FDR < 0.05 | * = significant",
    x        = NULL,
    y        = NULL
  )

ggsave(file.path(fig_out, "enrichment_dotplot_combined.pdf"),
       p_combined, width = 16, height = 8, device = "pdf")
ggsave(file.path(fig_out, "enrichment_dotplot_combined.tiff"),
       p_combined, width = 16, height = 8, device = "tiff",
       dpi = 300, compression = "lzw")

message(">>> Done.")
message(">>> Figures: ", fig_out)
message(">>> Tables:  ", tbl_out)
