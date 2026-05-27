## ============================================================
## 3_enrichment.R — GSEA + ORA (gene_id-based)
## Color rules + ratio labels
## ============================================================


suppressPackageStartupMessages({
  library(tidyverse)
  library(fgsea)
  library(msigdbr)
  library(clusterProfiler)
  library(enrichplot)
  library(org.Mm.eg.db)
})

set.seed(1)
message(">>> 3_enrichment running")

## -------------------- Parameters --------------------
n_cores     <- 4
padj_cutoff <- 0.05

analysis_level <- "gene"          # used only for universe filename; enrichment uses gene_id
species_name   <- "Mus musculus"  # manual

col_up   <- "#D62728"
col_down <- "#1F77B4"

strip_version <- function(x) sub("\\..*$", "", x)

## -------------------- Parallel backend --------------------
BiocParallel::register(BiocParallel::SnowParam(workers = n_cores, type = "SOCK"))
options(mc.cores = n_cores)

## -------------------- Directories --------------------
setwd("/hpcnfs/data/BA/glucocorticoid_keep/sinaR/m.rezaei")
project_root <- getwd()
data_dir     <- file.path(project_root, "data")

de_tbl   <- file.path(project_root, "results", "1_DEAnalysis", "tables", "shrunk")
int_tbl  <- file.path(project_root, "results", "2_interaction", "tables")

res_dir  <- file.path(project_root, "results", "3_enrichment")
tbl_dir  <- file.path(res_dir, "tables")
fig_dir  <- file.path(res_dir, "figures")
dir.create(tbl_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

## -------------------- Parameters --------------------
padj_cutoff <- 0.05

## -------------------- Annotation --------------------
annotation <- read.delim(file.path(project_root, "annotation.txt"), check.names = FALSE)
annotation <- annotation %>%
  dplyr::select(gene_id, symbol) %>%
  distinct(gene_id, .keep_all = TRUE)
id2sym <- setNames(annotation$symbol, annotation$gene_id)

## -------------------- Universe (REQUIRED) --------------------
universe_fp <- file.path(project_root, "results", "1_DEAnalysis", "tables",
                         paste0("expressed_universe_", analysis_level, ".tsv"))
if (!file.exists(universe_fp)) stop("Missing universe: ", universe_fp)

u <- read.delim(universe_fp, check.names = FALSE)
if (ncol(u) < 1) stop("Universe file has no columns: ", universe_fp)

universe_ens <- unique(na.omit(strip_version(u[[1]])))
if (length(universe_ens) < 50) stop("Universe too small (<50): ", universe_fp)

## Universe (ENTREZ) for KEGG ORA
universe_entrez <- mapIds(
  org.Mm.eg.db,
  keys      = universe_ens,
  column    = "ENTREZID",
  keytype   = "ENSEMBL",
  multiVals = "first"
) %>% unname() %>% na.omit() %>% unique()

## -------------------- MSigDB (for GSEA + ORA) --------------------
message("Loading MSigDB...")
msigdb <- list(
  HALLMARK = msigdbr::msigdbr(species = "Mus musculus", category = "H"),
  GOBP     = msigdbr::msigdbr(species = "Mus musculus", category = "C5", subcategory = "BP"),
  C2CP     = msigdbr::msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP")
)

term2gene_list <- lapply(msigdb, function(df) {
  df %>%
    dplyr::transmute(gs_name, ensembl_gene = strip_version(ensembl_gene)) %>%
    dplyr::distinct()
})

term_sizes_u <- lapply(term2gene_list, function(t2g) {
  t2g %>%
    dplyr::filter(ensembl_gene %in% universe_ens) %>%
    dplyr::count(gs_name, name = "setSize")
})

gmt_list <- lapply(term2gene_list, function(t2g) split(t2g$ensembl_gene, t2g$gs_name))

## -------------------- KEGG gene sets for GSEA (ENSEMBL IDs) --------------------
build_kegg_gmt <- function(universe_ens, organism = "mmu") {
  k <- clusterProfiler::download_KEGG(species = organism)
  
  # Join pathway IDs to readable names
  id2name <- k$KEGGPATHID2NAME %>%
    tibble::as_tibble() %>%
    dplyr::transmute(from, pathway_name = to)
  
  t2g_entrez <- k$KEGGPATHID2EXTID %>%
    tibble::as_tibble() %>%
    dplyr::left_join(id2name, by = "from") %>%
    dplyr::transmute(gs_name = ifelse(!is.na(pathway_name), pathway_name, from),
                     ENTREZID = as.character(to)) %>%
    dplyr::filter(!is.na(ENTREZID))
  
  # Keep names on the vector for correct lookup
  entrez_to_ens <- mapIds(
    org.Mm.eg.db,
    keys      = unique(t2g_entrez$ENTREZID),
    keytype   = "ENTREZID",
    column    = "ENSEMBL",
    multiVals = "first"
  )
  
  t2g_ens <- t2g_entrez %>%
    dplyr::mutate(ensembl_gene = strip_version(entrez_to_ens[ENTREZID])) %>%
    dplyr::filter(!is.na(ensembl_gene), ensembl_gene %in% universe_ens) %>%
    dplyr::distinct(gs_name, ensembl_gene)
  
  split(t2g_ens$ensembl_gene, t2g_ens$gs_name)
}

message("Loading KEGG gene sets for GSEA...")
kegg_gmt <- build_kegg_gmt(universe_ens, organism = "mmu")
message("KEGG pathways loaded: ", length(kegg_gmt))

## -------------------- Theme --------------------
theme_enrich <- function() {
  theme_bw(9) +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.text.y = element_text(size = 8)
    )
}

## -------------------- Plots --------------------
plot_gsea_bar <- function(gres_full, out_base, title_text) {
  dfp <- gres_full %>%
    dplyr::filter(!is.na(padj), padj < padj_cutoff)
  if (nrow(dfp) == 0) return(invisible(NULL))
  
  dfp <- dfp %>%
    dplyr::mutate(
      leadingEdgeCount = vapply(leadingEdge, length, integer(1)),
      ratio_label = paste0(leadingEdgeCount, "/", size)
    ) %>%
    dplyr::arrange(dplyr::desc(abs(NES))) %>%
    dplyr::slice_head(n = 20) %>%
    dplyr::mutate(
      pathway = forcats::fct_reorder(pathway, NES),
      dir = ifelse(NES > 0, "Up", "Down")
    )
  
  p <- ggplot(dfp, aes(x = NES, y = pathway, fill = dir)) +
    geom_col() +
    geom_text(aes(x = NES / 2, label = ratio_label),
              color = "black", size = 3.2, fontface = "bold") +
    scale_fill_manual(values = c(Up = col_up, Down = col_down),
                      name = "Direction") +
    labs(title = title_text, x = "NES", y = NULL) +
    theme_enrich()
  
  ggsave(paste0(out_base, ".pdf"),  p, width = 7, height = 5)
  ggsave(paste0(out_base, ".tiff"), p, width = 7, height = 5, dpi = 300,
         device = "tiff", compression = "lzw")
}

plot_ora_bar <- function(df_sig, out_base, fill_col) {
  if (nrow(df_sig) == 0) return(invisible(NULL))
  
  dfp <- df_sig %>%
    dplyr::arrange(p.adjust) %>%
    dplyr::mutate(Description = forcats::fct_reorder(Description, GeneRatio))
  
  p <- ggplot(dfp, aes(x = GeneRatio, y = Description)) +
    geom_col(fill = fill_col) +
    geom_text(aes(x = GeneRatio / 2, label = RatioLabel),
              size = 3.2, color = "black", fontface = "bold") +
    labs(x = "GeneRatio", y = NULL) +
    theme_enrich()
  
  ggsave(paste0(out_base, ".pdf"),  p, width = 7, height = 5)
  ggsave(paste0(out_base, ".tiff"), p, width = 7, height = 5, dpi = 300,
         device = "tiff", compression = "lzw")
}

## -------------------- Helpers --------------------
write_gsea_table <- function(gres, out_tsv) {
  if (is.null(gres) || nrow(gres) == 0) return(invisible(NULL))
  
  out_tbl <- gres %>%
    dplyr::arrange(padj, dplyr::desc(NES)) %>%
    dplyr::transmute(
      pathway = .data$pathway,
      size, ES, NES, pval, padj,
      leadingEdge
    )
  
  if (!is.null(id2sym)) {
    out_tbl <- out_tbl %>%
      dplyr::mutate(
        leadingEdge_symbols = vapply(
          leadingEdge,
          function(ids) paste(na.omit(id2sym[strip_version(ids)]), collapse = ";"),
          character(1)
        )
      )
  }
  
  out_tbl <- out_tbl %>%
    dplyr::mutate(leadingEdge = vapply(leadingEdge, paste, collapse = ";", character(1)))
  
  write.table(out_tbl, out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
  out_tbl
}

## ============================================================
## 1) GSEA (MSigDB + KEGG)
## ============================================================
tt_files <- c(
  list.files(de_tbl,  pattern = "^tT_.*\\.tsv$", full.names = TRUE),
  list.files(int_tbl, pattern = "^tT_.*\\.tsv$", full.names = TRUE)
)
if (length(tt_files) == 0) stop("No tT_*.tsv files found.")

for (fp in tt_files) {
  tag_raw <- tools::file_path_sans_ext(basename(fp))
  tag <- sub("^tT_", "", tag_raw)
  message("GSEA: ", tag)
  
  df <- read.delim(fp, check.names = FALSE)
  if (!("gene_id" %in% names(df))) stop("Missing gene_id in: ", basename(fp))
  
  rank_col <- if ("stat" %in% names(df)) "stat" else if ("log2FoldChange" %in% names(df)) "log2FoldChange" else NA
  if (is.na(rank_col)) stop("Need stat or log2FoldChange in: ", basename(fp))
  
  df2 <- df %>%
    dplyr::transmute(gene_id = strip_version(gene_id), rank = .data[[rank_col]]) %>%
    dplyr::filter(!is.na(gene_id), is.finite(rank), gene_id %in% universe_ens) %>%
    dplyr::group_by(gene_id) %>%
    dplyr::summarise(rank = rank[which.max(abs(rank))], .groups = "drop")
  
  ranks <- setNames(df2$rank, df2$gene_id)
  if (length(ranks) < 50) next
  
  ## ---- MSigDB GSEA ----
  for (set_tag in names(gmt_list)) {
    gres <- tryCatch(
      fgsea::fgseaMultilevel(
        pathways = gmt_list[[set_tag]],
        stats    = ranks,
        BPPARAM  = BiocParallel::bpparam()
      ),
      error = function(e) NULL
    )
    if (is.null(gres) || nrow(gres) == 0) next
    
    out_tsv <- file.path(tbl_dir, paste0("GSEA_", set_tag, "_", tag, ".tsv"))
    write_gsea_table(gres, out_tsv)
    
    out_base <- file.path(fig_dir, paste0("GSEA_", set_tag, "_", tag))
    plot_gsea_bar(
      gres_full  = gres %>% dplyr::transmute(pathway, NES, padj, size, leadingEdge),
      out_base   = out_base,
      title_text = paste0(set_tag, " - ", tag)
    )
  }
  
  ## ---- KEGG GSEA ----
  gres_kegg <- tryCatch(
    fgsea::fgseaMultilevel(
      pathways = kegg_gmt,
      stats    = ranks,
      BPPARAM  = BiocParallel::bpparam()
    ),
    error = function(e) NULL
  )
  
  if (!is.null(gres_kegg) && nrow(gres_kegg) > 0) {
    out_tsv <- file.path(tbl_dir, paste0("GSEA_KEGG_", tag, ".tsv"))
    write_gsea_table(gres_kegg, out_tsv)
    
    out_base <- file.path(fig_dir, paste0("GSEA_KEGG_", tag))
    plot_gsea_bar(
      gres_full  = gres_kegg %>% dplyr::transmute(pathway, NES, padj, size, leadingEdge),
      out_base   = out_base,
      title_text = paste0("KEGG - ", tag)
    )
  }
}

## ============================================================
## 2) ORA (MSigDB) — x=GeneRatio, label=Count/setSize
## ============================================================
up_files <- c(
  list.files(de_tbl,  pattern = "^up_.*\\.tsv$", full.names = TRUE),
  list.files(int_tbl, pattern = "^up_.*\\.tsv$", full.names = TRUE)
)
down_files <- c(
  list.files(de_tbl,  pattern = "^down_.*\\.tsv$", full.names = TRUE),
  list.files(int_tbl, pattern = "^down_.*\\.tsv$", full.names = TRUE)
)
ora_files <- c(up_files, down_files)
if (length(ora_files) == 0) stop("No up_/down_ tables found.")

for (fp in ora_files) {
  tag <- tools::file_path_sans_ext(basename(fp))
  message("ORA: ", tag)
  
  df <- read.delim(fp, check.names = FALSE)
  if (!("gene_id" %in% names(df))) stop("Missing gene_id in: ", basename(fp))
  
  genes_ens <- unique(na.omit(strip_version(df$gene_id)))
  if (length(genes_ens) < 10) next
  
  fill_col <- if (grepl("^up_", tag)) col_up else col_down
  inputN <- length(genes_ens)
  
  for (set_name in names(term2gene_list)) {
    t2g <- term2gene_list[[set_name]]
    
    res <- tryCatch(
      clusterProfiler::enricher(
        gene      = genes_ens,
        universe  = universe_ens,
        TERM2GENE = t2g,
        pAdjustMethod = "BH",
        qvalueCutoff  = padj_cutoff
      ),
      error = function(e) NULL
    )
    if (is.null(res)) next
    
    df_out <- as.data.frame(res)
    if (nrow(df_out) == 0) next
    
    df_out <- df_out %>%
      dplyr::left_join(term_sizes_u[[set_name]], by = c("ID" = "gs_name")) %>%
      dplyr::mutate(
        setSize    = tidyr::replace_na(setSize, 0L),
        GeneRatio  = Count / pmax(inputN, 1L),
        RatioLabel = paste0(Count, "/", setSize)
      )
    
    if (!is.null(id2sym)) {
      df_out$SYMBOLS <- vapply(
        strsplit(df_out$geneID, "/"),
        function(ids) paste(na.omit(id2sym[strip_version(ids)]), collapse = ";"),
        character(1)
      )
    }
    
    out_tsv <- file.path(tbl_dir, paste0("ORA_", set_name, "_", tag, ".tsv"))
    write.table(df_out, out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
    
    df_sig <- df_out %>% dplyr::filter(p.adjust < padj_cutoff)
    if (nrow(df_sig) == 0) next
    
    out_base <- file.path(fig_dir, paste0("ORA_", set_name, "_", tag))
    plot_ora_bar(df_sig, out_base, fill_col = fill_col)
  }
}

## ============================================================
## 3) KEGG ORA (enrichKEGG) — x=GeneRatio, label=Count/setSize
## ============================================================
for (fp in ora_files) {
  tag <- tools::file_path_sans_ext(basename(fp))
  message("KEGG ORA: ", tag)
  
  df <- read.delim(fp, check.names = FALSE)
  if (!("gene_id" %in% names(df))) stop("Missing gene_id in: ", basename(fp))
  
  genes_ens <- unique(na.omit(strip_version(df$gene_id)))
  if (length(genes_ens) < 10) next
  
  fill_col <- if (grepl("^up_", tag)) col_up else col_down
  
  genes_entrez <- mapIds(
    org.Mm.eg.db,
    keys      = genes_ens,
    column    = "ENTREZID",
    keytype   = "ENSEMBL",
    multiVals = "first"
  ) %>% unname() %>% na.omit() %>% unique()
  
  if (length(genes_entrez) < 10 || length(universe_entrez) < 50) next
  
  kegg <- tryCatch(
    clusterProfiler::enrichKEGG(
      gene          = genes_entrez,
      universe      = universe_entrez,
      organism      = "mmu",
      pAdjustMethod = "BH",
      qvalueCutoff  = padj_cutoff
    ),
    error = function(e) NULL
  )
  if (is.null(kegg)) next
  
  df_out <- as.data.frame(kegg)
  if (nrow(df_out) == 0) next
  
  inputN <- length(genes_entrez)
  
  df_out <- df_out %>%
    dplyr::mutate(
      GeneRatio  = Count / pmax(inputN, 1L),
      setSize    = as.integer(sub("/.*$", "", BgRatio)),
      RatioLabel = paste0(Count, "/", setSize)
    )
  
  out_tsv <- file.path(tbl_dir, paste0("ORA_KEGG_", tag, ".tsv"))
  write.table(df_out, out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
  
  df_sig <- df_out %>% dplyr::filter(p.adjust < padj_cutoff)
  if (nrow(df_sig) == 0) next
  
  out_base <- file.path(fig_dir, paste0("ORA_KEGG_", tag))
  plot_ora_bar(df_sig, out_base, fill_col = fill_col)
}



## ===================================
## 4) Heatmap from results of GSEA
##    3 structured heatmaps per drug per gene-set category
## ===================================

message(">>> Building GSEA heatmaps (structured, 3 per drug)...")
library(pheatmap)
library(grid)
library(dplyr)
library(stringr)
library(tibble)

## ---- Drugs and timepoints ----
drugs      <- c("Budesonide", "Prednisolone")
timepoints <- c("2", "4", "8", "16")

## ---- Helper: load one GSEA file and return pathway/NES/padj ----
load_gsea <- function(tag, tbl_dir, categories) {
  out <- list()
  for (cat in categories) {
    fp <- file.path(tbl_dir, paste0("GSEA_", cat, "_", tag, ".tsv"))
    if (!file.exists(fp)) next
    df <- tryCatch(read.delim(fp, check.names = FALSE), error = function(e) NULL)
    if (is.null(df)) next
    if (!all(c("pathway", "NES", "padj") %in% colnames(df))) next
    out[[cat]] <- df %>% dplyr::select(pathway, NES, padj)
  }
  out
}

## ---- Helper: clean pathway names ----
clean_pathway_names <- function(nms) {
  nms <- gsub("^HALLMARK_|^GOBP_|^KEGG_|^REACTOME_|^WP_", "", nms)
  nms <- gsub("_", " ", nms)
  nms <- tolower(nms)
  nms <- paste0(toupper(substr(nms, 1, 1)), substr(nms, 2, nchar(nms)))
  stringr::str_wrap(nms, width = 45)
}

## ---- Helper: build NES + star matrices from a named list of single-column data ----
##   col_data: named list, each element = list(df=data.frame(pathway,NES,padj), label=string)
build_matrices <- function(col_data, padj_cutoff) {
  col_labels <- sapply(col_data, `[[`, "label")
  dfs        <- lapply(col_data, `[[`, "df")
  
  sig_terms <- unique(unlist(lapply(dfs, function(df) {
    if (is.null(df)) return(character(0))
    df$pathway[!is.na(df$padj) & df$padj < padj_cutoff]
  })))
  
  if (length(sig_terms) == 0) return(NULL)
  
  mat      <- matrix(NA_real_, nrow = length(sig_terms), ncol = length(col_data),
                     dimnames = list(sig_terms, col_labels))
  star_mat <- matrix("",      nrow = length(sig_terms), ncol = length(col_data),
                     dimnames = list(sig_terms, col_labels))
  
  for (i in seq_along(col_data)) {
    df  <- dfs[[i]]
    lbl <- col_labels[i]
    if (is.null(df)) next
    hit <- df$pathway %in% sig_terms
    mat[df$pathway[hit], lbl]      <- df$NES[hit]
    sig_hit <- hit & !is.na(df$padj) & df$padj < padj_cutoff
    star_mat[df$pathway[sig_hit], lbl] <- "★"
  }
  
  ## Keep only rows with ≥1 significant result
  has_star <- apply(star_mat, 1, function(x) any(x == "★"))
  mat      <- mat[has_star, , drop = FALSE]
  star_mat <- star_mat[has_star, , drop = FALSE]
  if (nrow(mat) == 0) return(NULL)
  
  ## Top 40 by max |NES|
  top_n   <- min(40, nrow(mat))
  top_idx <- order(apply(abs(mat), 1, max, na.rm = TRUE), decreasing = TRUE)[seq_len(top_n)]
  mat      <- mat[top_idx, , drop = FALSE]
  star_mat <- star_mat[top_idx, , drop = FALSE]
  
  ## Clean row names
  clean <- clean_pathway_names(rownames(mat))
  rownames(mat)      <- clean
  rownames(star_mat) <- clean
  
  list(mat = mat, star_mat = star_mat)
}

## ---- Helper: save a pheatmap ----
save_heatmap <- function(mat, star_mat, out_base, title_text,
                         max_abs = NULL, tbl_dir = NULL, tbl_name = NULL) {
  if (is.null(max_abs) || !is.finite(max_abs) || max_abs == 0)
    max_abs <- max(abs(mat), na.rm = TRUE)
  if (!is.finite(max_abs) || max_abs == 0) max_abs <- 1
  
  breaks <- seq(-max_abs, max_abs, length.out = 101)
  pal    <- colorRampPalette(c("#2166AC", "white", "#B2182B"))(100)
  
  n_rows <- nrow(mat)
  n_cols <- ncol(mat)
  fig_h  <- max(5, n_rows * 0.28 + 3)
  fig_w  <- max(8, n_cols * 0.45 + 5)
  
  ph <- pheatmap::pheatmap(
    mat,
    color             = pal,
    breaks            = breaks,
    main              = title_text,
    cluster_rows      = n_rows >= 2,
    cluster_cols      = FALSE,
    treeheight_row    = 30,
    treeheight_col    = 0,
    clustering_method = "ward.D2",
    fontsize_row      = 7,
    fontsize_col      = 8,
    fontsize          = 8,
    angle_col         = 45,
    na_col            = "grey92",
    display_numbers   = star_mat,
    number_color      = "black",
    fontsize_number   = 9,
    legend_breaks     = c(-max_abs, -max_abs/2, 0, max_abs/2, max_abs),
    legend_labels     = round(c(-max_abs, -max_abs/2, 0, max_abs/2, max_abs), 2),
    border_color      = NA,
    cellwidth         = 22,
    cellheight        = 14,
    silent            = TRUE
  )
  
  ggsave(paste0(out_base, ".pdf"),  ph, width = fig_w, height = fig_h, limitsize = FALSE)
  ggsave(paste0(out_base, ".tiff"), ph, width = fig_w, height = fig_h, limitsize = FALSE,
         dpi = 300, device = "tiff", compression = "lzw")
  
  if (!is.null(tbl_dir) && !is.null(tbl_name)) {
    write.table(
      as.data.frame(mat) %>% tibble::rownames_to_column("pathway"),
      file.path(tbl_dir, paste0(tbl_name, ".tsv")),
      sep = "\t", quote = FALSE, row.names = FALSE
    )
  }
  
  message("    Saved: ", basename(out_base),
          "  (", n_rows, " pathways × ", n_cols, " contrasts)")
}

## ---- Detect available categories from existing GSEA files ----
gsea_files <- list.files(tbl_dir, pattern = "^GSEA_.*\\.tsv$", full.names = TRUE)
if (length(gsea_files) == 0) {
  message("No GSEA files found.")
} else {
  
  ## Parse category from filenames: GSEA_<CATEGORY>_<tag>.tsv
  categories <- unique(sub("^GSEA_([^_]+)_.*\\.tsv$", "\\1", basename(gsea_files)))
  
  ## ---- OHT vs DMSO (shared reference column, one for all drugs/timepoints) ----
  ##   Tag: DMSO_OHT_vs_DMSO
  oht_tag   <- "DMSO_OHT_vs_DMSO"
  oht_label <- "OHT vs DMSO"
  
  ## ---- Loop: drug → category → 3 heatmaps ----
  for (drug in drugs) {
    message("  [Drug] ", drug)
    
    ## Build per-timepoint tag lookup
    ## drug vs DMSO:        e.g. Budesonide_8
    ## drug+OHT vs DMSO:    e.g. Budesonide_8_OHT_vs_DMSO
    ## drug+OHT vs OHT:     e.g. Budesonide_8_OHT
    ## drug+OHT vs drug:    e.g. Budesonide_8_OHT__Budesonide_8
    
    for (cat in categories) {
      message("    [Category] ", cat)
      
      ## -- Load OHT vs DMSO (shared) --
      oht_data <- load_gsea(oht_tag, tbl_dir, cat)[[cat]]  # may be NULL
      
      ## -- Load per-timepoint contrasts --
      tp_drug_dmso     <- list()  # drug vs DMSO
      tp_drugOHT_dmso  <- list()  # drug+OHT vs DMSO
      tp_drugOHT_oht   <- list()  # drug+OHT vs OHT
      tp_drugOHT_drug  <- list()  # drug+OHT vs drug
      
      for (tp in timepoints) {
        tag_drug_dmso    <- paste0(drug, "_", tp)
        tag_drugOHT_dmso <- paste0(drug, "_", tp, "_OHT_vs_DMSO")
        tag_drugOHT_oht  <- paste0(drug, "_", tp, "_OHT")
        tag_drugOHT_drug <- paste0(drug, "_", tp, "_OHT__", drug, "_", tp)
        
        tp_drug_dmso[[tp]]    <- load_gsea(tag_drug_dmso,    tbl_dir, cat)[[cat]]
        tp_drugOHT_dmso[[tp]] <- load_gsea(tag_drugOHT_dmso, tbl_dir, cat)[[cat]]
        tp_drugOHT_oht[[tp]]  <- load_gsea(tag_drugOHT_oht,  tbl_dir, cat)[[cat]]
        tp_drugOHT_drug[[tp]] <- load_gsea(tag_drugOHT_drug,  tbl_dir, cat)[[cat]]
      }
      
      ## ============================================================
      ## HEATMAP 1:
      ##   OHT vs DMSO (0h) | drug_T vs DMSO (×4 tp) | drug+OHT_T vs DMSO (×4 tp)
      ## ============================================================
      col_data_hm1 <- list()
      
      ## Column 1: OHT vs DMSO
      col_data_hm1[["oht"]] <- list(df = oht_data, label = oht_label)
      
      ## Columns 2-5: drug vs DMSO per timepoint
      for (tp in timepoints) {
        key <- paste0("drug_dmso_", tp)
        col_data_hm1[[key]] <- list(
          df    = tp_drug_dmso[[tp]],
          label = paste0(drug, " T", tp, " vs DMSO")
        )
      }
      
      ## Columns 6-9: drug+OHT vs DMSO per timepoint
      for (tp in timepoints) {
        key <- paste0("drugOHT_dmso_", tp)
        col_data_hm1[[key]] <- list(
          df    = tp_drugOHT_dmso[[tp]],
          label = paste0(drug, "+OHT T", tp, " vs DMSO")
        )
      }
      
      res1 <- build_matrices(col_data_hm1, padj_cutoff)
      if (!is.null(res1)) {
        out_base <- file.path(fig_dir, paste0("GSEA_", cat, "_", drug, "_HM1"))
        save_heatmap(
          mat       = res1$mat,
          star_mat  = res1$star_mat,
          out_base  = out_base,
          title_text = paste0(cat, "  —  ", drug, "  —  HM1: OHT | Drug vs DMSO | Drug+OHT vs DMSO"),
          tbl_dir   = tbl_dir,
          tbl_name  = paste0("GSEA_", cat, "_", drug, "_HM1_NES_matrix")
        )
      } else {
        message("      HM1: no significant pathways, skipped")
      }
      
      ## ============================================================
      ## HEATMAP 2:
      ##   drug_T vs DMSO (×4 tp) | drug+OHT_T vs OHT (×4 tp)
      ## ============================================================
      col_data_hm2 <- list()
      
      ## Columns 1-4: drug vs DMSO per timepoint
      for (tp in timepoints) {
        key <- paste0("drug_dmso_", tp)
        col_data_hm2[[key]] <- list(
          df    = tp_drug_dmso[[tp]],
          label = paste0(drug, " T", tp, " vs DMSO")
        )
      }
      
      ## Columns 5-8: drug+OHT vs OHT per timepoint
      for (tp in timepoints) {
        key <- paste0("drugOHT_oht_", tp)
        col_data_hm2[[key]] <- list(
          df    = tp_drugOHT_oht[[tp]],
          label = paste0(drug, "+OHT T", tp, " vs OHT")
        )
      }
      
      res2 <- build_matrices(col_data_hm2, padj_cutoff)
      if (!is.null(res2)) {
        out_base <- file.path(fig_dir, paste0("GSEA_", cat, "_", drug, "_HM2"))
        save_heatmap(
          mat        = res2$mat,
          star_mat   = res2$star_mat,
          out_base   = out_base,
          title_text = paste0(cat, "  —  ", drug, "  —  HM2: Drug vs DMSO | Drug+OHT vs OHT"),
          tbl_dir    = tbl_dir,
          tbl_name   = paste0("GSEA_", cat, "_", drug, "_HM2_NES_matrix")
        )
      } else {
        message("      HM2: no significant pathways, skipped")
      }
      
      ## ============================================================
      ## HEATMAP 3:
      ##   OHT vs DMSO (0h) | drug+OHT_T vs drug_T (×4 tp)
      ## ============================================================
      col_data_hm3 <- list()
      
      ## Column 1: OHT vs DMSO
      col_data_hm3[["oht"]] <- list(df = oht_data, label = oht_label)
      
      ## Columns 2-5: drug+OHT vs drug per timepoint
      for (tp in timepoints) {
        key <- paste0("drugOHT_drug_", tp)
        col_data_hm3[[key]] <- list(
          df    = tp_drugOHT_drug[[tp]],
          label = paste0(drug, "+OHT T", tp, " vs ", drug)
        )
      }
      
      res3 <- build_matrices(col_data_hm3, padj_cutoff)
      if (!is.null(res3)) {
        out_base <- file.path(fig_dir, paste0("GSEA_", cat, "_", drug, "_HM3"))
        save_heatmap(
          mat        = res3$mat,
          star_mat   = res3$star_mat,
          out_base   = out_base,
          title_text = paste0(cat, "  —  ", drug, "  —  HM3: OHT | Drug+OHT vs Drug"),
          tbl_dir    = tbl_dir,
          tbl_name   = paste0("GSEA_", cat, "_", drug, "_HM3_NES_matrix")
        )
      } else {
        message("      HM3: no significant pathways, skipped")
      }
      
    }  # end category loop
  }  # end drug loop
}

message(">>> All GSEA heatmaps done.")

## -------------------- Reproducibility --------------------
writeLines(capture.output(sessionInfo()), file.path(res_dir, "sessionInfo.txt"))

message(">>> done")
message("Tables: ", tbl_dir)
message("Figures: ", fig_dir)

