############################################################
## Phydget_Kidney_GO.R
##
## Analyze PhyDGET Kidney results:
##  - Identify OGs with bestmodel in five_*, four_*, three_*
##  - Identify OGs with bestmodel in five_*, four_*, three_*, two_*, one_*
##    AND BF(bestmodel) > 1.5
##  - Map OGs -> representative gene symbols
##  - Build a Kidney-specific background universe
##  - Run human GO BP enrichment for both sets
##  - Save everything to RDS + write gene & GO tables
############################################################

## --------------------------
## 0. Setup
## --------------------------

## Set working directory to the same place as 01_DE_and_GO.R
setwd("/n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/nf_pipeline/thaDol_mapping_test")

library(DESeq2)
library(SummarizedExperiment)
library(dplyr)
library(readr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)

## Input PhyDGET result file (Kidney)
phydget_file <- "phydget/Phydget_Kidney_dry_vs_mesic_results_all_possible_dry_comparisons_more_time.txt"

## Pre-existing objects produced by 01_DE_and_GO.R
## (these should already exist in this directory)
dds_path        <- "obj_dds_DESeq2_OG_all_tissues.rds"
og_rep_all_path <- "obj_og_rep_ALL_OGs_repGeneName.rds"

stopifnot(file.exists(phydget_file),
          file.exists(dds_path),
          file.exists(og_rep_all_path))

## Load DESeq2 object and global OG->gene_name mapping
dds        <- readRDS(dds_path)        # DESeqDataSet with all tissues
og_rep_all <- readRDS(og_rep_all_path) # data.frame with columns: OG, gene_name, name_source

############################################################
## 1. Read PhyDGET Kidney results
############################################################

phyd <- read.delim(phydget_file, header = TRUE, check.names = FALSE)

## Check we have key columns
stopifnot("gene" %in% colnames(phyd),
          "bestmodel" %in% colnames(phyd))

## Convenience: map BF column name per row
## Example: if bestmodel == "five_1", BF column is "BF.five_1"
bf_colname <- paste0("BF.", phyd$bestmodel)

## Check that required BF columns exist where needed
unique_missing_bf <- setdiff(unique(bf_colname), colnames(phyd))
if (length(unique_missing_bf) > 0) {
  warning("Some BF.* columns referenced by bestmodel are missing from the PhyDGET file:\n",
          paste(unique_missing_bf, collapse = ", "))
}

## Get numeric BF for each row where possible
bf_col_index <- match(bf_colname, colnames(phyd))
## bf_vals[i] = BF corresponding to bestmodel in row i (or NA if not found)
bf_vals <- rep(NA_real_, nrow(phyd))
valid_idx <- which(!is.na(bf_col_index))
bf_vals[valid_idx] <- phyd[cbind(valid_idx, bf_col_index[valid_idx])]

############################################################
## 2. Define foreground OG sets
############################################################

## 2a) Set A: bestmodel in five_*, four_*, three_* (any BF)
idx_setA <- grepl("^(five|four|three)_", phyd$bestmodel)
ogs_setA <- phyd$gene[idx_setA]

## 2b) Set B: bestmodel in five_*, four_*, three_*, two_*, one_* AND BF(bestmodel) > 1.5
idx_model_B <- grepl("^(five|four|three|two|one)_", phyd$bestmodel)
idx_BF_B    <- !is.na(bf_vals) & bf_vals > 1.5
idx_setB    <- idx_model_B & idx_BF_B
ogs_setB    <- phyd$gene[idx_setB]

message("Set A (bestmodel = five_/four_/three_): ", length(ogs_setA), " OGs")
message("Set B (bestmodel = five_/four_/three_/two_/one_ & BF > 1.5): ", length(ogs_setB), " OGs")

## Optional: assemble small data.frames for inspection
phyd_setA <- phyd[idx_setA, c("gene", "bestmodel")]
phyd_setB <- data.frame(
  gene      = phyd$gene[idx_setB],
  bestmodel = phyd$bestmodel[idx_setB],
  BF_best   = bf_vals[idx_setB],
  stringsAsFactors = FALSE
)

############################################################
## 3. Map OGs -> representative gene_name (symbols)
############################################################

## og_rep_all has columns: OG, gene_name, name_source
stopifnot(all(c("OG", "gene_name") %in% colnames(og_rep_all)))

## Helper: map OGs to gene_name, drop NAs/empty
map_ogs_to_symbols <- function(ogs, og_rep_all) {
  og_rep_all %>%
    filter(OG %in% ogs, !is.na(gene_name), gene_name != "") %>%
    distinct(OG, gene_name, .keep_all = TRUE)
}

map_setA <- map_ogs_to_symbols(ogs_setA, og_rep_all)
map_setB <- map_ogs_to_symbols(ogs_setB, og_rep_all)

symbols_setA <- unique(map_setA$gene_name)
symbols_setB <- unique(map_setB$gene_name)

message("Set A mapped to ", length(symbols_setA), " unique gene_name symbols")
message("Set B mapped to ", length(symbols_setB), " unique gene_name symbols")

## Write OG lists and OG->gene_name tables
write.table(ogs_setA,
            file = "phydget_Kidney_OGs_bestmodel_543_all.txt",
            quote = FALSE, row.names = FALSE, col.names = FALSE)

write.table(ogs_setB,
            file = "phydget_Kidney_OGs_bestmodel_54321_BFgt1.5.txt",
            quote = FALSE, row.names = FALSE, col.names = FALSE)

write.csv(map_setA,
          file = "phydget_Kidney_OGs_bestmodel_543_all_with_gene_names.csv",
          row.names = FALSE)
write.csv(map_setB,
          file = "phydget_Kidney_OGs_bestmodel_54321_BFgt1.5_with_gene_names.csv",
          row.names = FALSE)

############################################################
## 4. Build Kidney-specific background universe (ENTREZ)
############################################################

## This matches the logic from run_go_for_tissue_dry in 01_DE_and_GO.R

## 4a) Identify Kidney samples in dds
stopifnot("Tissue" %in% colnames(colData(dds)))
cols_kidney <- which(dds$Tissue == "Kidney")

if (length(cols_kidney) == 0) {
  stop("No Kidney samples found in dds$Tissue")
}

cts_kidney <- counts(dds)[, cols_kidney, drop = FALSE]
bg_ogs     <- rownames(cts_kidney)[rowSums(cts_kidney > 0) > 0]

message("Kidney background OGs with >0 counts: ", length(bg_ogs))

if (length(bg_ogs) < 5) {
  stop("Too few background OGs for Kidney; cannot run GO.")
}

bg_map <- og_rep_all %>%
  filter(OG %in% bg_ogs, !is.na(gene_name), gene_name != "")
bg_symbols <- unique(bg_map$gene_name)

message("Kidney background mapped to ", length(bg_symbols), " gene_name symbols")

if (length(bg_symbols) < 5) {
  stop("Too few background symbols for Kidney; cannot run GO.")
}

## 4b) Convert background symbols to ENTREZ IDs
bg_entrez <- tryCatch(
  bitr(bg_symbols,
       fromType = "SYMBOL",
       toType   = "ENTREZID",
       OrgDb    = org.Hs.eg.db),
  error = function(e) {
    stop("bitr conversion failed for Kidney background: ", e$message)
  }
)

if (is.null(bg_entrez) || nrow(bg_entrez) == 0) {
  stop("No background symbols converted to Entrez for Kidney.")
}

bg_universe <- unique(bg_entrez$ENTREZID)
message("Kidney GO universe size (ENTREZ): ", length(bg_universe))

############################################################
## 5. Run GO BP enrichment for Set A and Set B (human) + dotplots
############################################################

run_go_for_symbol_set <- function(symbols,
                                  bg_universe,
                                  tag,
                                  OrgDb = org.Hs.eg.db,
                                  ont = "BP",
                                  pcut = 0.05) {
  message("\n=== Running GO for ", tag, " ===")
  symbols <- unique(symbols)
  message("Foreground symbols: ", length(symbols))
  if (length(symbols) < 5) {
    warning("Too few foreground symbols for ", tag, " (", length(symbols), "); skipping GO.")
    return(NULL)
  }
  
  fg_conv <- tryCatch(
    bitr(symbols,
         fromType = "SYMBOL",
         toType   = "ENTREZID",
         OrgDb    = OrgDb),
    error = function(e) {
      warning("bitr conversion failed for foreground in ", tag, ": ", e$message)
      return(NULL)
    }
  )
  if (is.null(fg_conv) || nrow(fg_conv) == 0) {
    warning("No foreground symbols converted to Entrez for ", tag)
    return(NULL)
  }
  
  ego <- tryCatch(
    enrichGO(
      gene          = fg_conv$ENTREZID,
      universe      = bg_universe,
      OrgDb         = OrgDb,
      keyType       = "ENTREZID",
      ont           = ont,
      pAdjustMethod = "BH",
      pvalueCutoff  = pcut,
      qvalueCutoff  = pcut,
      readable      = TRUE
    ),
    error = function(e) {
      warning("enrichGO failed for ", tag, ": ", e$message)
      return(NULL)
    }
  )
  
  if (is.null(ego) || nrow(as.data.frame(ego)) == 0) {
    warning("No enriched GO terms for ", tag)
    return(NULL)
  }
  
  ## Write table
  ego_df <- as.data.frame(ego)
  out_tab <- paste0("GO_BP_Hs_", tag, "_enrichment.csv")
  write.csv(ego_df, out_tab, row.names = FALSE)
  message("  Wrote GO table: ", out_tab)
  
  ## Dotplot (top 20 by default)
  p_dot <- dotplot(ego, showCategory = 20) +
    ggtitle(paste("GO BP –", tag)) +
    theme_bw() +
    theme(
      plot.title   = element_text(size = 11),
      axis.text.x  = element_text(size = 8),
      axis.text.y  = element_text(size = 8),
      axis.title.x = element_text(size = 9),
      axis.title.y = element_text(size = 9)
    )
  
  out_dot <- paste0("GO_BP_Hs_", tag, "_dotplot.png")
  ggsave(out_dot, p_dot, width = 8, height = 6, dpi = 300)
  message("  Wrote dotplot: ", out_dot)
  
  ego
}

## GO for Set A
ego_setA <- run_go_for_symbol_set(
  symbols     = symbols_setA,
  bg_universe = bg_universe,
  tag         = "Kidney_Phydget_bestmodel_543_all"
)

## GO for Set B (now includes two_ *)
ego_setB <- run_go_for_symbol_set(
  symbols     = symbols_setB,
  bg_universe = bg_universe,
  tag         = "Kidney_Phydget_bestmodel_54321_BFgt1.5"
)

############################################################
## 6. Save everything to an RDS
############################################################

phydget_kidney_results <- list(
  phydget_file    = phydget_file,
  phyd_raw        = phyd,
  bf_vals         = bf_vals,
  ## foreground definitions
  setA = list(
    description       = "bestmodel in five_*, four_*, three_*",
    ogs               = ogs_setA,
    phydget_rows      = phyd_setA,
    og_to_symbol      = map_setA,
    fg_symbols        = symbols_setA,
    go_result         = ego_setA
  ),
  setB = list(
    description       = "bestmodel in five_*, four_*, three_*, two_*, one_* and BF(bestmodel) > 1.5",
    ogs               = ogs_setB,
    phydget_rows      = phyd_setB,
    og_to_symbol      = map_setB,
    fg_symbols        = symbols_setB,
    go_result         = ego_setB
  ),
  ## background universe
  kidney_background = list(
    bg_ogs       = bg_ogs,
    bg_symbols   = bg_symbols,
    bg_universe  = bg_universe
  )
)

saveRDS(phydget_kidney_results,
        file = "phydget_Kidney_GO_results.rds")

message("\nSaved full PhyDGET Kidney GO results to phydget_Kidney_GO_results.rds")

############################################################
## End of Phydget_Kidney_GO.R
############################################################