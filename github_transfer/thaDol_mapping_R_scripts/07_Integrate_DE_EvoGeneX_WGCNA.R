############################################################
## 07_Integrate_DE_EvoGeneX_WGCNA.R
##
## Integrate:
##  - Convergent DE (Dry vs Mesic) from conv_ogs
##  - EvoGeneX adaptive OGs
##  - WGCNA top 3 Dry-associated modules
##
## For each tissue:
##  - Build a table of OGs with flags:
##       ConvDryUp, ConvDryDown, Adaptive, WGCNA_top3
##  - Attach best/common gene_name (og_rep_all)
##  - Save per-tissue CSV and UpSet plot
## Also save a combined table across all tissues.
############################################################

setwd("/n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/nf_pipeline/thaDol_mapping_test")

library(dplyr)
library(tidyr)
library(readr)
library(ComplexUpset)
library(ggplot2)

`%||%` <- function(x, y) if (is.null(x)) y else x

tissues <- c("Brain","Heart","Liver","Muscle","Kidney")

############################################################
## 0) Load core objects and settings
############################################################

## Global OG -> gene_name map
if (file.exists("obj_og_rep_ALL_OGs_repGeneName.rds")) {
  og_rep_all <- readRDS("obj_og_rep_ALL_OGs_repGeneName.rds")
} else if (file.exists("OGs_ALL_repGeneName.csv")) {
  og_rep_all <- read.csv("OGs_ALL_repGeneName.csv", stringsAsFactors = FALSE)
} else {
  stop("Need og_rep_all (obj_og_rep_ALL_OGs_repGeneName.rds or OGs_ALL_repGeneName.csv).")
}
stopifnot(all(c("OG","gene_name","name_source") %in% colnames(og_rep_all)))

## Convergent DE in dry (from 01_DE_and_GO.R)
if (!file.exists("obj_conv_ogs_byTissue_dryDir.rds")) {
  stop("obj_conv_ogs_byTissue_dryDir.rds not found; run 01_DE_and_GO.R first.")
}
conv_ogs <- readRDS("obj_conv_ogs_byTissue_dryDir.rds")

## WGCNA top 3 Dry-associated modules (from 05_WGCNA_all_tissues_and_GO.R)
top3_all <- read.csv("WGCNA_Top3Modules_Dry_byTissue.csv",
                     stringsAsFactors = FALSE)

############################################################
## 1) Helper: get EvoGeneX adaptive OGs for a tissue
############################################################

get_adaptive_ogs <- function(tiss) {
  evog_file <- paste0("EvoGeneX_", tiss, "_results.csv")
  if (!file.exists(evog_file)) {
    warning("No EvoGeneX results for ", tiss, " at ", evog_file)
    return(character(0))
  }
  evog_res <- read.csv(evog_file, stringsAsFactors = FALSE)
  if (!all(c("OG","class","Tissue") %in% colnames(evog_res))) {
    warning("EvoGeneX results for ", tiss, " missing OG/class/Tissue columns.")
    return(character(0))
  }
  evog_res %>%
    dplyr::filter(Tissue == tiss, class == "adaptive") %>%
    dplyr::pull(OG) %>%
    unique()
}

############################################################
## 2) Build per-tissue integration tables
############################################################

all_tissue_tables <- list()

for (tiss in tissues) {
  message("\n=== Tissue: ", tiss, " ===")
  
  ## 2.1 Convergent DE OGs
  conv_t <- conv_ogs[[tiss]]
  if (is.null(conv_t)) {
    message("  No conv_ogs entry for ", tiss, "; treating as empty.")
    conv_up   <- character(0)
    conv_down <- character(0)
  } else {
    conv_up   <- conv_t$dry_up   %||% character(0)
    conv_down <- conv_t$dry_down %||% character(0)
  }
  
  ## 2.2 Adaptive OGs (EvoGeneX)
  adapt_ogs <- get_adaptive_ogs(tiss)
  message("  # adaptive OGs: ", length(adapt_ogs))
  
  ## 2.3 WGCNA top 3 modules
  top3_tiss <- subset(top3_all, Tissue == tiss)
  if (nrow(top3_tiss) == 0) {
    warning("  No top3 WGCNA modules recorded for tissue ", tiss)
    top_colors <- character(0)
  } else {
    top_colors <- unique(top3_tiss$Color)
  }
  message("  Top3 WGCNA module colors: ", paste(top_colors, collapse = ", "))
  
  ## Load geneInfo to get OGs in top3 modules
  gi_file <- paste0("WGCNA_geneInfo_", tiss, "_Dry_vs_Mesic.csv")
  if (!file.exists(gi_file)) {
    warning("  geneInfo file not found for tissue ", tiss, ": ", gi_file)
    og_wgcna <- character(0)
    mod_colors_by_og <- tibble(OG = character(0), WGCNA_modules = character(0))
  } else {
    gi <- read.csv(gi_file, stringsAsFactors = FALSE)
    if (!all(c("OG","moduleColor") %in% colnames(gi))) {
      warning("  geneInfo for ", tiss, " missing OG/moduleColor.")
      og_wgcna <- character(0)
      mod_colors_by_og <- tibble(OG = character(0), WGCNA_modules = character(0))
    } else {
      gi_top <- gi %>%
        dplyr::filter(moduleColor %in% top_colors)
      og_wgcna <- unique(gi_top$OG)
      # which top-3 colors each OG belongs to
      mod_colors_by_og <- gi_top %>%
        dplyr::group_by(OG) %>%
        dplyr::summarise(
          WGCNA_modules = paste(unique(moduleColor), collapse = ";"),
          .groups = "drop"
        )
    }
  }
  message("  # OGs in top3 WGCNA modules: ", length(og_wgcna))
  
  ## 2.4 Union of all OGs that appear in any of these categories for this tissue
  og_union <- unique(c(conv_up, conv_down, adapt_ogs, og_wgcna))
  if (length(og_union) == 0) {
    warning("  No OGs to annotate for tissue ", tiss, "; skipping.")
    next
  }
  
  ## 2.5 Build flag table
  df_tiss <- tibble(
    Tissue        = tiss,
    OG            = og_union,
    ConvDryUp     = as.integer(OG %in% conv_up),
    ConvDryDown   = as.integer(OG %in% conv_down),
    Adaptive      = as.integer(OG %in% adapt_ogs),
    WGCNA_top3    = as.integer(OG %in% og_wgcna)
  )
  
  ## 2.6 Attach WGCNA module color(s)
  if (length(og_wgcna) > 0 && nrow(mod_colors_by_og) > 0) {
    df_tiss <- df_tiss %>%
      dplyr::left_join(mod_colors_by_og, by = "OG")
  } else {
    df_tiss$WGCNA_modules <- NA_character_
  }
  
  ## 2.7 Attach best/common gene_name from og_rep_all
  df_tiss <- df_tiss %>%
    dplyr::left_join(og_rep_all, by = "OG") %>%
    dplyr::relocate(gene_name, name_source, .after = OG)
  
  ## 2.8 Save per-tissue CSV
  out_csv <- paste0("IntegratedFlags_", tiss,
                    "_ConvDE_Adaptive_WGCNA_top3.csv")
  write.csv(df_tiss, out_csv, row.names = FALSE)
  message("  Wrote integrated table: ", out_csv,
          "  (n = ", nrow(df_tiss), " OGs)")
  
  ## 2.9 UpSet plot for this tissue
  # Use only OGs that are in at least one category
  df_sets <- df_tiss %>%
    dplyr::mutate(any_set = ConvDryUp + ConvDryDown + Adaptive + WGCNA_top3) %>%
    dplyr::filter(any_set > 0) %>%
    dplyr::select(-any_set)
  
  if (nrow(df_sets) > 0) {
    set_cols <- c("ConvDryUp","ConvDryDown","Adaptive","WGCNA_top3")
    
    pdf_file <- paste0("UpSet_Integrated_", tiss,
                       "_ConvDE_Adaptive_WGCNA_top3.pdf")
    pdf(pdf_file, width = 7, height = 5)
    print(
      upset(
        df_sets,
        set_cols,
        name = "OGs",
        min_size = 1
      ) +
        ggtitle(paste("Intersection of categories –", tiss))
    )
    dev.off()
    message("  Wrote UpSet plot: ", pdf_file)
  }
  
  all_tissue_tables[[tiss]] <- df_tiss
}

############################################################
## 3) Combined table across all tissues
############################################################

if (length(all_tissue_tables) > 0) {
  integrated_all <- bind_rows(all_tissue_tables)
  out_all <- "IntegratedFlags_ALL_tissues_ConvDE_Adaptive_WGCNA_top3.csv"
  write.csv(integrated_all, out_all, row.names = FALSE)
  message("\nWrote combined integrated table: ", out_all)
} else {
  warning("No tissue tables were generated; nothing to combine.")
}

############################################################
## End of 07_Integrate_DE_EvoGeneX_WGCNA.R
############################################################



integrated_all <- read.csv("IntegratedFlags_ALL_tissues_ConvDE_Adaptive_WGCNA_top3.csv")

# Example: adaptive + convergently DE + in WGCNA top3, per tissue
hits <- integrated_all %>%
  filter(Adaptive == 1, (ConvDryUp == 1 | ConvDryDown == 1), WGCNA_top3 == 1) %>%
  arrange(Tissue, gene_name)

head(hits)



# Save to file
write.csv(hits,
          "Hits_Adaptive_ConvDry_WGCNA_top3_with_geneNames.csv",
          row.names = FALSE)





hits <- read.csv("Hits_Adaptive_ConvDry_WGCNA_top3_with_geneNames.csv",
                 stringsAsFactors = FALSE)

kidney_genes <- hits %>%
  dplyr::filter(Tissue == "Kidney") %>%
  dplyr::pull(gene_name) %>%
  unique()

kidney_genes

write.table(kidney_genes,
            file = "Hits_Adaptive_ConvDry_WGCNA_top3_Kidney_geneNames.txt",
            quote = FALSE, row.names = FALSE, col.names = FALSE)
