## ====================================================================
## ONE-BLOCK SCRIPT: load GO CSVs and plot convergence heatmaps
## ====================================================================

## 0) Set working directory to your project
setwd("/n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/nf_pipeline/thaDol_mapping_test")

## 1) Libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(pheatmap)
})

## 2) Helper: load all GO BP CSVs for one tissue (ALL directions)
##    Files are named like:
##    GO_BP_Hs_<Tissue>_<pair_id>_ALLdir_pair_enrichment.csv
##    e.g. GO_BP_Hs_Muscle_sakLuc_vs_sakCan_ALLdir_pair_enrichment.csv
load_go_csvs_all_pairs <- function(tissue,
                                   prefix = "GO_BP_Hs_",
                                   suffix = "_ALLdir_pair_enrichment.csv") {
  
  pattern <- paste0("^", prefix, tissue, "_.*", suffix, "$")
  files <- list.files(pattern = pattern)
  if (length(files) == 0) {
    stop("No GO CSVs found for tissue '", tissue, "' with pattern ", pattern)
  }
  
  message("Found ", length(files), " GO files for tissue ", tissue, ":")
  print(files)
  
  go_list <- list()
  
  for (f in files) {
    # Extract pair_id from filename:
    # GO_BP_Hs_<Tissue>_<pair_id>_ALLdir_pair_enrichment.csv
    # Remove prefix and suffix
    middle <- sub(paste0("^", prefix, tissue, "_"), "", f)
    pair_id <- sub(suffix, "", middle)
    
    df <- suppressMessages(read_csv(f, show_col_types = FALSE))
    
    ## Expect columns: ID, Description, GeneRatio, BgRatio, pvalue, p.adjust, ...
    ## If they have no header, adjust here, but from your example they do.
    if (!all(c("ID", "Description", "p.adjust") %in% colnames(df))) {
      warning("File ", f, " missing expected columns; skipping.")
      next
    }
    
    go_list[[pair_id]] <- df
  }
  
  if (length(go_list) == 0) {
    stop("No usable GO tables loaded for tissue ", tissue)
  }
  
  go_list
}

## 3) Build GO overlap matrix for one tissue
##    rows: GO IDs, cols: pairs, values: -log10(p.adjust) (0 if not sig, NA if absent)
build_go_overlap_matrix_from_csv <- function(tissue,
                                             alpha = 0.05,
                                             top_n_terms = 50) {
  
  go_list <- load_go_csvs_all_pairs(tissue)
  
  # Collect all GO IDs across all pairs
  all_terms <- unique(unlist(lapply(go_list, function(df) df$ID)))
  if (length(all_terms) == 0) {
    stop("No GO IDs found across all pairs for tissue ", tissue)
  }
  
  pairs <- names(go_list)
  mat <- matrix(NA_real_,
                nrow = length(all_terms),
                ncol = length(pairs),
                dimnames = list(all_terms, pairs))
  
  for (pid in pairs) {
    df <- go_list[[pid]]
    rownames(df) <- df$ID
    common_ids <- intersect(all_terms, df$ID)
    if (length(common_ids) == 0) next
    
    padj <- df[common_ids, "p.adjust", drop = TRUE]
    
    vals <- rep(NA_real_, length(common_ids))
    names(vals) <- common_ids
    
    sig_idx    <- which(!is.na(padj) & padj <= alpha)
    nonsig_idx <- which(!is.na(padj) & padj > alpha)
    
    if (length(sig_idx) > 0) {
      vals[sig_idx] <- -log10(padj[sig_idx])
    }
    if (length(nonsig_idx) > 0) {
      vals[nonsig_idx] <- 0
    }
    
    mat[names(vals), pid] <- vals
  }
  
  # Keep top GO terms by how many pairs they are enriched in
  sig_counts <- rowSums(mat > 0, na.rm = TRUE)
  keep_ids   <- names(sort(sig_counts, decreasing = TRUE))[1:min(top_n_terms,
                                                                 length(sig_counts))]
  mat_filt   <- mat[keep_ids, , drop = FALSE]
  
  # Map GO IDs to Descriptions for nicer row names
  id_desc_df <- do.call(
    rbind,
    lapply(go_list, function(df) df[, c("ID", "Description")])
  ) %>% distinct()
  
  desc_map <- id_desc_df$Description
  names(desc_map) <- id_desc_df$ID
  
  ids <- rownames(mat_filt)
  row_labels <- ifelse(ids %in% names(desc_map), desc_map[ids], ids)
  dup <- duplicated(row_labels)
  if (any(dup)) {
    row_labels[dup] <- paste0(row_labels[dup], " (", ids[dup], ")")
  }
  rownames(mat_filt) <- row_labels
  
  list(
    matrix   = mat_filt,
    raw_list = go_list
  )
}

## 4) Plot heatmap of GO overlap for one tissue
plot_go_overlap_heatmap_from_csv <- function(tissue,
                                             alpha = 0.05,
                                             top_n_terms = 50) {
  res <- build_go_overlap_matrix_from_csv(
    tissue      = tissue,
    alpha       = alpha,
    top_n_terms = top_n_terms
  )
  mat_filt <- res$matrix
  
  has_pos <- any(mat_filt > 0, na.rm = TRUE)
  
  if (!has_pos) {
    message("No significant GO terms for ", tissue,
            " (all cells are 0 or NA). Plotting presence/absence only.")
    breaks <- c(-0.5, 0.5)
    colors <- c("grey80")
    pheatmap(mat_filt,
             color        = colors,
             breaks       = breaks,
             na_col       = "white",
             cluster_rows = TRUE,
             cluster_cols = TRUE,
             main         = paste("GO BP overlap –", tissue,
                                  "\n0 = overlap but padj > ", alpha,
                                  "; >0 none"),
             fontsize_row = 6,
             fontsize_col = 8,
             border_color = NA)
  } else {
    max_val <- max(mat_filt, na.rm = TRUE)
    if (!is.finite(max_val) || max_val <= 0) max_val <- 1
    
    seq_part <- seq(0.00001, max_val, length.out = 100)
    breaks <- unique(c(0, 0.00001, seq_part))
    colors <- c("grey80", colorRampPalette(c("yellow", "red"))(length(breaks) - 2))
    
    pheatmap(mat_filt,
             color        = colors,
             breaks       = breaks,
             na_col       = "white",
             cluster_rows = TRUE,
             cluster_cols = TRUE,
             main         = paste("GO BP overlap –", tissue,
                                  "\n0 = overlap but padj > ", alpha,
                                  "; color = -log10(padj)"),
             fontsize_row = 6,
             fontsize_col = 8,
             border_color = NA)
  }
}

## 5) EXAMPLE: draw Muscle heatmap for ALL-direction GO
##    This will show things like:
##    - whether 'response to lipid', 'regulation of insulin secretion', etc.,
##      are enriched in both sakLuc_vs_sakCan and other pairs.
plot_go_overlap_heatmap_from_csv(
  tissue      = "Muscle",
  alpha       = 0.05,
  top_n_terms = 50
)

## You can similarly do:
plot_go_overlap_heatmap_from_csv(
  tissue      = "Liver",
  alpha       = 0.05,
  top_n_terms = 50
)




plot_go_overlap_heatmap_from_csv(
  tissue      = "Kidney",
  alpha       = 0.05,
  top_n_terms = 50
)


plot_go_overlap_heatmap_from_csv(
  tissue      = "Brain",
  alpha       = 0.05,
  top_n_terms = 50
)


plot_go_overlap_heatmap_from_csv(
  tissue      = "Heart",
  alpha       = 0.05,
  top_n_terms = 50
)





# try with less strict 0.1

plot_go_overlap_heatmap_from_csv(
  tissue      = "Muscle",
  alpha       = 0.1,
  top_n_terms = 50
)

## You can similarly do:
plot_go_overlap_heatmap_from_csv(
  tissue      = "Liver",
  alpha       = 0.1,
  top_n_terms = 50
)




plot_go_overlap_heatmap_from_csv(
  tissue      = "Kidney",
  alpha       = 0.1,
  top_n_terms = 50
)


plot_go_overlap_heatmap_from_csv(
  tissue      = "Brain",
  alpha       = 0.1,
  top_n_terms = 50
)


plot_go_overlap_heatmap_from_csv(
  tissue      = "Heart",
  alpha       = 0.05,
  top_n_terms = 50
)



