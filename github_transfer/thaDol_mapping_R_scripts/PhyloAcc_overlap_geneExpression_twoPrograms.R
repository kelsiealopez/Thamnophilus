# want to specifically look at gene expression genes that have been detected by two or more programs and overlap with phyloacc 



#!/usr/bin/env Rscript

############################################################
## GO for intersect(phyloacc, expression) genes
## that are supported by ≥ 2 methods in:
##   IntegratedFlags_ALL_tissues_ConvDE_Adaptive_WGCNA_top3_PhyDGET.csv
##
## Methods/columns:
##   ConvDryUp, ConvDryDown, Adaptive, WGCNA_top3, PhyDGET_BFgt1.5
##
## Foreground:
##   gene_name in overlap list AND row has sum(method_flags) >= 2
##   (then collapsed to unique gene_name across tissues/OGs)
##
## Background:
##   OGs expressed in Brain/Heart/Kidney/Liver/Muscle (from dds),
##   mapped via og_rep_all to human symbols.
############################################################


  library(DESeq2)          # for colData(), counts()
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(dplyr)
  library(readr)
  library(ggplot2)


## -------------------------
## 1. File paths
## -------------------------

## Overlap gene symbols between phyloacc & expression (1177 genes; 1 per line)
overlap_gene_file <- "/n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/cactus-snakemake/thamnophilus-all-species-cactus/thamnophilus-all-species-cactus_output/phyloacc_test_900k_output/overlapping_gene_symbols_ALLtissues_vs_phyloacc.txt"

## Integrated flags file (ALL tissues)
integrated_flags_csv <- "/n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/nf_pipeline/thaDol_mapping_test/IntegratedFlags_ALL_tissues_ConvDE_Adaptive_WGCNA_top3_PhyDGET.csv"

## R objects produced by 01_DE_and_GO.R
dds_rds        <- "/n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/nf_pipeline/thaDol_mapping_test/obj_dds_DESeq2_OG_all_tissues.rds"
og_rep_all_rds <- "/n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/nf_pipeline/thaDol_mapping_test/obj_og_rep_ALL_OGs_repGeneName.rds"

## Output files
out_prefix   <- "/n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/cactus-snakemake/thamnophilus-all-species-cactus/thamnophilus-all-species-cactus_output/phyloacc_test_900k_output/GO_overlap_phyloacc_expression_2methods"
out_table    <- paste0(out_prefix, "_BP_Hs_enrichment.csv")
out_dotplot  <- paste0(out_prefix, "_BP_Hs_dotplot.png")
out_fg_genes <- paste0(out_prefix, "_foreground_genes.txt")

## -------------------------
## 2. Load overlap gene symbols
## -------------------------

if (!file.exists(overlap_gene_file)) {
  stop("Overlap gene file not found: ", overlap_gene_file)
}

overlap_symbols <- readLines(overlap_gene_file)
overlap_symbols <- overlap_symbols[nzchar(overlap_symbols)]
overlap_symbols <- unique(overlap_symbols)

cat("Total overlap symbols (phyloacc ∩ expression):\n")
cat("  N =", length(overlap_symbols), "\n\n")

if (length(overlap_symbols) < 5) {
  stop("Too few overlapping genes (", length(overlap_symbols), ") for GO.")
}

## -------------------------
## 3. Load IntegratedFlags and select rows with ≥ 2 methods = 1
## -------------------------

if (!file.exists(integrated_flags_csv)) {
  stop("Integrated flags CSV not found: ", integrated_flags_csv)
}

int_flags <- read.csv(integrated_flags_csv, stringsAsFactors = FALSE)

required_cols <- c("Tissue", "OG", "gene_name", "name_source",
                   "ConvDryUp", "ConvDryDown", "Adaptive",
                   "WGCNA_top3", "PhyDGET_BFgt1.5")
if (!all(required_cols %in% colnames(int_flags))) {
  stop("Integrated flags CSV is missing required columns.\nNeeded: ",
       paste(required_cols, collapse = ", "))
}

## Ensure flag columns are numeric (0/1)
flag_cols <- c("ConvDryUp", "ConvDryDown", "Adaptive",
               "WGCNA_top3", "PhyDGET_BFgt1.5")

int_flags[flag_cols] <- lapply(int_flags[flag_cols], as.numeric)

## Number of methods "on" per row
int_flags$method_sum <- rowSums(int_flags[flag_cols], na.rm = TRUE)

## Keep rows where (a) gene_name in overlap list, (b) at least 2 methods = 1
int_overlap_2meth <- int_flags %>%
  filter(
    !is.na(gene_name),
    gene_name != "",
    gene_name %in% overlap_symbols,
    method_sum >= 2
  )

cat("Rows in IntegratedFlags where gene in overlap AND ≥2 methods = 1:\n")
cat("  N rows =", nrow(int_overlap_2meth), "\n\n")

if (nrow(int_overlap_2meth) == 0) {
  stop("No rows meet the criteria (overlap & ≥2 methods).")
}

## Foreground: unique gene_name from these rows
fg_symbols <- unique(int_overlap_2meth$gene_name)

cat("Foreground symbols (overlap & ≥2 methods, collapsed by gene_name):\n")
cat("  N =", length(fg_symbols), "\n\n")

## Save foreground gene list (for inspection)
writeLines(fg_symbols, con = out_fg_genes)
cat("Wrote foreground gene list to:\n  ", out_fg_genes, "\n\n")

if (length(fg_symbols) < 5) {
  stop("Too few foreground genes after ≥2 methods filter (", length(fg_symbols), ") for GO.")
}

## -------------------------
## 4. Load dds + global OG->gene_name map
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
## 5. Build background universe
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
## 6. Restrict foreground to symbols present in background
## -------------------------

fg_symbols_filt <- intersect(fg_symbols, bg_symbols)
cat("Foreground symbols that are also in background universe:\n")
cat("  N =", length(fg_symbols_filt), "\n\n")

if (length(fg_symbols_filt) < 5) {
  stop("Too few foreground genes after restricting to background universe (", 
       length(fg_symbols_filt), ").")
}

## -------------------------
## 7. Map foreground and background to human Entrez IDs
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
## 8. Run enrichGO (BP, human)
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
}

## -------------------------
## 9. Dotplot (RStudio + PNG)
## -------------------------

if (nrow(ego_df) > 0) {
  p_dot <- dotplot(ego, showCategory = 20) +
    ggtitle("GO BP – phyloacc ∩ expression (≥2 methods)") +
    theme_bw() +
    theme(
      plot.title   = element_text(size = 11),
      axis.text.x  = element_text(size = 8),
      axis.text.y  = element_text(size = 8),
      axis.title.x = element_text(size = 9),
      axis.title.y = element_text(size = 9)
    )
  
  ## Show in RStudio (interactive) if run there
  if (interactive()) {
    print(p_dot)
  }
  
  ## Also save to file
  ggsave(out_dotplot, p_dot, width = 8, height = 6, dpi = 300)
  cat("Wrote GO BP dotplot to:\n  ", out_dotplot, "\n\n")
}

cat("Done.\n")
