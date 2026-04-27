setwd("/n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/nf_pipeline/thaDol_mapping_test")



############################################################
## GO for genes overlapping phyloacc + expression
## - Foreground: genes in BOTH datasets
## - Background: OGs expressed in Brain/Heart/Kidney/Liver/Muscle
############################################################


  library(DESeq2)          # for colData(), counts(), etc.
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(dplyr)
  library(readr)
  library(ggplot2)

## -------------------------
## 1. File paths
## -------------------------

## Overlap gene symbols (1177 genes; 1 per line)
overlap_gene_file <- "/n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/cactus-snakemake/thamnophilus-all-species-cactus/thamnophilus-all-species-cactus_output/phyloacc_test_900k_output/overlapping_gene_symbols_ALLtissues_vs_phyloacc.txt"

## R objects produced by 01_DE_and_GO.R
dds_rds        <- "/n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/nf_pipeline/thaDol_mapping_test/obj_dds_DESeq2_OG_all_tissues.rds"
og_rep_all_rds <- "/n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/nf_pipeline/thaDol_mapping_test/obj_og_rep_ALL_OGs_repGeneName.rds"

## Output files
out_prefix   <- "/n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/cactus-snakemake/thamnophilus-all-species-cactus/thamnophilus-all-species-cactus_output/phyloacc_test_900k_output/GO_overlap_phyloacc_expression"
out_table    <- paste0(out_prefix, "_BP_Hs_enrichment.csv")
out_dotplot  <- paste0(out_prefix, "_BP_Hs_dotplot.png")

## -------------------------
## 2. Load foreground (overlap genes)
## -------------------------

if (!file.exists(overlap_gene_file)) {
  stop("Overlap gene file not found: ", overlap_gene_file)
}

fg_symbols <- readLines(overlap_gene_file)
fg_symbols <- fg_symbols[nzchar(fg_symbols)]
fg_symbols <- unique(fg_symbols)

cat("Foreground (overlap) symbols:\n")
cat("  N =", length(fg_symbols), "\n\n")

if (length(fg_symbols) < 5) {
  stop("Too few overlapping genes (", length(fg_symbols), ") for GO.")
}

## -------------------------
## 3. Load dds + global OG->gene_name map
## -------------------------

if (!file.exists(dds_rds)) {
  stop("dds RDS not found: ", dds_rds)
}
if (!file.exists(og_rep_all_rds)) {
  stop("og_rep_all RDS not found: ", og_rep_all_rds)
}

dds        <- readRDS(dds_rds)
og_rep_all <- readRDS(og_rep_all_rds)

## -------------------------
## 4. Build background universe
##    OGs expressed in Brain/Heart/Kidney/Liver/Muscle
## -------------------------

tissues_use <- c("Brain", "Heart", "Kidney", "Liver", "Muscle")

if (!"Tissue" %in% colnames(colData(dds))) {
  stop("dds must have a 'Tissue' column in colData.")
}

cols_tiss <- which(colData(dds)$Tissue %in% tissues_use)
if (length(cols_tiss) == 0) {
  stop("No samples in dds with Tissue in: ", paste(tissues_use, collapse = ", "))
}

cts_tiss <- counts(dds)[, cols_tiss, drop = FALSE]

## OGs with any nonzero counts in those tissues
bg_ogs <- rownames(cts_tiss)[rowSums(cts_tiss > 0) > 0]
cat("Background OGs (expressed in selected tissues):\n")
cat("  N =", length(bg_ogs), "\n\n")

if (length(bg_ogs) < 5) {
  stop("Too few background OGs (", length(bg_ogs), ") for GO.")
}

## Map OGs -> representative gene_name using og_rep_all
## og_rep_all has columns: OG, gene_name, name_source
bg_symbols <- og_rep_all %>%
  filter(OG %in% bg_ogs,
         !is.na(gene_name),
         gene_name != "") %>%
  pull(gene_name) %>%
  unique()

cat("Background gene_name symbols (mapped from OGs):\n")
cat("  N =", length(bg_symbols), "\n\n")

if (length(bg_symbols) < 5) {
  stop("Too few background gene symbols (", length(bg_symbols), ") for GO.")
}

## -------------------------
## 5. Restrict foreground to symbols present in background (optional but tidy)
## -------------------------

fg_symbols_filt <- intersect(fg_symbols, bg_symbols)
cat("Foreground symbols that are also in background:\n")
cat("  N =", length(fg_symbols_filt), "\n\n")

if (length(fg_symbols_filt) < 5) {
  stop("Too few overlapping genes after restricting to background universe (", 
       length(fg_symbols_filt), ").")
}

## -------------------------
## 6. Map foreground and background to human Entrez IDs
## -------------------------

fg_conv <- tryCatch(
  bitr(fg_symbols_filt,
       fromType = "SYMBOL",
       toType   = "ENTREZID",
       OrgDb    = org.Hs.eg.db),
  error = function(e) {
    stop("bitr failed for FOREGROUND: ", e$message)
  }
)

if (is.null(fg_conv) || nrow(fg_conv) == 0) {
  stop("No foreground symbols could be converted to Entrez IDs.")
}

bg_conv <- tryCatch(
  bitr(bg_symbols,
       fromType = "SYMBOL",
       toType   = "ENTREZID",
       OrgDb    = org.Hs.eg.db),
  error = function(e) {
    stop("bitr failed for BACKGROUND: ", e$message)
  }
)

if (is.null(bg_conv) || nrow(bg_conv) == 0) {
  stop("No background symbols could be converted to Entrez IDs.")
}

fg_entrez <- unique(fg_conv$ENTREZID)
bg_entrez <- unique(bg_conv$ENTREZID)

cat("Foreground Entrez IDs: N =", length(fg_entrez), "\n")
cat("Background Entrez IDs: N =", length(bg_entrez), "\n\n")

if (length(fg_entrez) < 5) {
  stop("Too few foreground Entrez IDs (", length(fg_entrez), ") for GO.")
}
if (length(bg_entrez) < 5) {
  stop("Too few background Entrez IDs (", length(bg_entrez), ") for GO.")
}

## -------------------------
## 7. Run enrichGO (BP, human)
## -------------------------

ego <- tryCatch(
  enrichGO(
    gene          = fg_entrez,
    universe      = bg_entrez,
    OrgDb         = org.Hs.eg.db,
    keyType       = "ENTREZID",
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.05,
    readable      = TRUE
  ),
  error = function(e) {
    stop("enrichGO failed: ", e$message)
  }
)

ego_df <- as.data.frame(ego)

if (nrow(ego_df) == 0) {
  warning("No GO BP terms enriched at padj <= 0.05.")
} else {
  write.csv(ego_df, out_table, row.names = FALSE)
  cat("Wrote GO BP enrichment table to:\n  ", out_table, "\n\n")
  
  ## Dotplot
  p_dot <- dotplot(ego, showCategory = 20) +
    ggtitle("GO BP – Overlap of phyloacc & expression genes") +
    theme_bw() +
    theme(
      plot.title   = element_text(size = 11),
      axis.text.x  = element_text(size = 8),
      axis.text.y  = element_text(size = 8),
      axis.title.x = element_text(size = 9),
      axis.title.y = element_text(size = 9)
    )
  
  ggsave(out_dotplot, p_dot, width = 8, height = 6, dpi = 300)
  cat("Wrote GO BP dotplot to:\n  ", out_dotplot, "\n\n")
}

cat("Done.\n")


############################################################
## Assuming you already have:
##   - ego: result from enrichGO(...)
##   - library(clusterProfiler); library(ggplot2)
############################################################

## 1. Basic dotplot (shows in RStudio Plots pane)
dotplot(ego, showCategory = 20) +
  ggtitle("GO BP – Overlap of phyloacc & expression genes") +
  theme_bw() +
  theme(
    plot.title   = element_text(size = 11),
    axis.text.x  = element_text(size = 8),
    axis.text.y  = element_text(size = 8),
    axis.title.x = element_text(size = 9),
    axis.title.y = element_text(size = 9)
  )
