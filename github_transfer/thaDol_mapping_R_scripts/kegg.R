###############################
## KEGG ORA for Gallus genes ##
###############################

## 0. Setup ------------------------------------------------------

## If needed, install packages once (uncomment and run)
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(c("clusterProfiler", "org.Gg.eg.db", "enrichplot"))

library(clusterProfiler)
library(org.Gg.eg.db)   # Gallus gallus annotation
library(enrichplot)
library(ggplot2)

## 1. Working directory and input file ---------------------------

## Set this to your phyloacc output directory
setwd("/n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/cactus-snakemake/thamnophilus-all-species-cactus/thamnophilus-all-species-cactus_output/phyloacc_test_900k_output")

## File with your gene symbols (one per line)
gene_file <- "results/top10pct_nearest_gene_symbols_10kb.txt"

if (!file.exists(gene_file)) {
  stop("Gene file not found: ", gene_file)
}

## 2. Read gene symbols and convert to Entrez IDs ---------------

genes_symbol <- scan(gene_file, what = "character", quiet = TRUE)
genes_symbol <- unique(genes_symbol[nzchar(genes_symbol)])

cat("Read", length(genes_symbol), "unique gene symbols\n")

if (length(genes_symbol) < 5) {
  stop("Too few genes (", length(genes_symbol), ") for KEGG enrichment.")
}

## Map SYMBOL -> ENTREZID using chicken OrgDb
conv <- bitr(
  genes_symbol,
  fromType = "SYMBOL",
  toType   = "ENTREZID",
  OrgDb    = org.Gg.eg.db
)

if (is.null(conv) || nrow(conv) == 0) {
  stop("No symbols could be converted to Entrez IDs.")
}

gene_entrez <- unique(conv$ENTREZID)
cat("Successfully mapped to", length(gene_entrez), "Entrez IDs\n")

if (length(gene_entrez) < 5) {
  stop("Too few mapped Entrez IDs (", length(gene_entrez), ") for KEGG.")
}

## 3. (Optional) background / universe ---------------------------
## If you have a background gene list, map it to Entrez and
## set universe = bg_entrez in enrichKEGG(). Otherwise, KEGG’s
## internal background will be used.

## Example (commented out):
# bg_file <- "results/background_gene_symbols.txt"
# bg_symbols <- scan(bg_file, what = "character", quiet = TRUE)
# bg_conv <- bitr(bg_symbols, fromType = "SYMBOL",
#                 toType   = "ENTREZID",
#                 OrgDb    = org.Gg.eg.db)
# bg_entrez <- unique(bg_conv$ENTREZID)

## 4. Run KEGG over-representation analysis ----------------------

ekegg <- enrichKEGG(
  gene         = gene_entrez,
  organism     = "gga",          # Gallus gallus
  keyType      = "ncbi-geneid",
  # universe   = bg_entrez,      # uncomment + define if you have a background
  pvalueCutoff = 0.05,
  pAdjustMethod= "BH",
  minGSSize    = 10,
  maxGSSize    = 500
)

ekegg_df <- as.data.frame(ekegg)

cat("Number of enriched KEGG pathways (p.adjust <= 0.05):",
    sum(ekegg_df$p.adjust <= 0.05, na.rm = TRUE), "\n")

## 5. Save results -----------------------------------------------

out_table  <- "KEGG_enrichment_top10pct_nearest_gene_symbols_10kb_gga.tsv"
out_plot   <- "KEGG_enrichment_top10pct_nearest_gene_symbols_10kb_gga_dotplot.png"

write.table(
  ekegg_df,
  file      = out_table,
  sep       = "\t",
  quote     = FALSE,
  row.names = FALSE
)

cat("Wrote KEGG enrichment table to:\n  ", out_table, "\n")

## 6. Plot a dotplot of top KEGG pathways ------------------------

if (nrow(ekegg_df) > 0) {
  p <- dotplot(ekegg, showCategory = 20) +
    ggtitle("KEGG enrichment – top10pct_nearest_gene_symbols_10kb (Gallus gallus)") +
    theme_bw()
  
  ggsave(out_plot, p, width = 8, height = 6, dpi = 300)
  cat("Wrote KEGG dotplot to:\n  ", out_plot, "\n")
} else {
  cat("No KEGG pathways enriched at the specified cutoff.\n")
}

cat("Done.\n")
