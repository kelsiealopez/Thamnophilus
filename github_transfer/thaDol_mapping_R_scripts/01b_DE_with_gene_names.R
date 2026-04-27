############################################################
## 01b_DE_with_gene_names.R
##
## Add best/common gene_name to DESeq2 pairwise DE tables
## produced by 01_DE_and_GO.R
##
## Input:
##   - DE_OG_selected_<spA>_vs_<spB>_<Tissue>_kallisto.csv
##   - og_rep_all (global OG -> gene_name mapping)
##
## Output (per input DE file):
##   - DE_OG_selected_<spA>_vs_<spB>_<Tissue>_kallisto_with_geneNames.csv
##       (all OGs, DE table + OG + gene_name + name_source)
##   - DE_SIG_OG_<spA>_vs_<spB>_<Tissue>_kallisto_with_geneNames.csv
##       (only padj < 0.05 & |log2FC| >= 1, with gene_name)
############################################################

setwd("/n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/nf_pipeline/thaDol_mapping_test")

library(dplyr)
library(readr)
library(stringr)

## 0) Load global OG -> gene_name mapping
og_rep_all <- if (file.exists("obj_og_rep_ALL_OGs_repGeneName.rds")) {
  readRDS("obj_og_rep_ALL_OGs_repGeneName.rds")
} else if (file.exists("OGs_ALL_repGeneName.csv")) {
  read.csv("OGs_ALL_repGeneName.csv", stringsAsFactors = FALSE)
} else {
  stop("Cannot find og_rep_all (obj_og_rep_ALL_OGs_repGeneName.rds or OGs_ALL_repGeneName.csv).")
}

stopifnot(all(c("OG","gene_name","name_source") %in% colnames(og_rep_all)))

## 1) Find all DE result files
de_files <- list.files(
  pattern = "^DE_OG_selected_.*_kallisto\\.csv$",
  full.names = TRUE
)

if (length(de_files) == 0) {
  stop("No DE_OG_selected_*_kallisto.csv files found in working directory.")
}

message("Found ", length(de_files), " DE result files.")

## 2) Process each file
for (f in de_files) {
  message("Processing: ", basename(f))
  
  # read DE table
  de_tab <- read_csv(f, show_col_types = FALSE)
  
  # first column is rownames from as.data.frame(res) in 01 script
  # make them explicit OG IDs
  if (!("" %in% colnames(de_tab))) {
    # read_csv drops "" column name; here rownames are in first column already
    # We can re-read with base read.csv to be sure
    de_tab_base <- read.csv(f, row.names = 1, check.names = FALSE)
    de_tab <- as.data.frame(de_tab_base)
    de_tab$OG <- rownames(de_tab_base)
  } else {
    # if "" exists (unlikely with readr), rename to OG
    names(de_tab)[names(de_tab) == ""] <- "OG"
  }
  
  if (!"OG" %in% colnames(de_tab)) {
    stop("Could not find OG column in ", f)
  }
  
  # join with og_rep_all to get gene_name
  de_annot <- de_tab %>%
    dplyr::left_join(og_rep_all, by = "OG")
  
  # Reorder: OG, gene_name, name_source, then the DE columns
  base_cols <- c("OG", "gene_name", "name_source")
  other_cols <- setdiff(colnames(de_annot), base_cols)
  de_annot <- de_annot[, c(base_cols, other_cols)]
  
  # Write full annotated file
  out_full <- sub("\\.csv$", "_with_geneNames.csv", f)
  write.csv(de_annot, out_full, row.names = FALSE)
  message("  Wrote annotated full DE file: ", basename(out_full))
  
  # Filter significant DE genes: padj < 0.05 & |log2FoldChange| >= 1
  if (!all(c("padj","log2FoldChange") %in% colnames(de_annot))) {
    warning("  Missing padj/log2FoldChange in ", basename(f),
            "; skipping sig-only file.")
    next
  }
  
  de_sig <- de_annot %>%
    dplyr::filter(!is.na(padj),
                  padj < 0.05,
                  !is.na(log2FoldChange),
                  abs(log2FoldChange) >= 1)
  
  if (nrow(de_sig) == 0) {
    message("  No significant DE OGs in ", basename(f),
            " at padj < 0.05 & |log2FC| >= 1.")
    next
  }
  
  # Construct sig-DE output filename
  # Example input: DE_OG_selected_dysSti_vs_sakCri_Brain_kallisto.csv
  # → DE_SIG_OG_dysSti_vs_sakCri_Brain_kallisto_with_geneNames.csv
  base_name <- basename(f)
  core      <- sub("^DE_OG_selected_", "", base_name)
  core      <- sub("\\.csv$", "", core)
  out_sig   <- paste0("DE_SIG_OG_", core, "_with_geneNames.csv")
  
  write.csv(de_sig, out_sig, row.names = FALSE)
  message("  Wrote annotated SIG-only DE file: ", out_sig,
          " (n = ", nrow(de_sig), ")")
}

message("Done annotating DE tables with gene names.")
