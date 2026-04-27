setwd("/n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/nf_pipeline/thaDol_mapping_test")

library(dplyr)
library(readr)
library(pheatmap)
library(RColorBrewer)

############################################################
## 0) Load vst expression + sample metadata
############################################################

vsd <- readRDS("obj_vsd_OG_all_tissues.rds")
vsd_mat <- assay(vsd)                      # OG x sample

coldata <- readRDS("obj_coldata_samples.rds")  # Sample, Species, Tissue

# Match order
stopifnot(all(colnames(vsd_mat) == coldata$Sample))

############################################################
## 1) Load WGCNA geneInfo + top3 modules
############################################################

top3_all <- read.csv("WGCNA_Top3Modules_Dry_byTissue.csv",
                     stringsAsFactors = FALSE)

tissues_to_run <- c("Brain","Heart","Liver","Muscle","Kidney")

# Define dry vs mesic species (for annotation)
dry_species <- c("thaDol", "sakCri", "sakCan", "thaBer", "thaPel", "thaTor")

############################################################
## 2) Helper: make a heatmap for one tissue x module
############################################################

plot_module_heatmap <- function(tiss, module_color,
                                vsd_mat, coldata,
                                max_genes = 100,
                                save_png = TRUE) {
  message("Heatmap for tissue ", tiss, " – module ", module_color)
  
  # 2.1 Load geneInfo for this tissue
  gi_file <- paste0("WGCNA_geneInfo_", tiss, "_Dry_vs_Mesic.csv")
  if (!file.exists(gi_file)) {
    warning("  geneInfo file not found for tissue ", tiss, ": ", gi_file)
    return(invisible(NULL))
  }
  gi <- read.csv(gi_file, stringsAsFactors = FALSE)
  
  # 2.2 Genes in requested module, with names
  gi_mod <- gi %>%
    dplyr::filter(moduleColor == module_color) %>%
    dplyr::mutate(has_name = !is.na(gene_name) & gene_name != "")
  
  if (nrow(gi_mod) == 0) {
    warning("  No genes in module ", module_color, " for tissue ", tiss)
    return(invisible(NULL))
  }
  
  # Optional: sort by |GS_Dry| so top rows are most Dry-associated
  gi_mod <- gi_mod %>%
    dplyr::arrange(desc(abs(GS_Dry)))
  
  # limit to max_genes for readability
  if (nrow(gi_mod) > max_genes) {
    gi_mod <- gi_mod[seq_len(max_genes), ]
  }
  
  og_ids <- gi_mod$OG
  
  # 2.3 Expression matrix for this tissue and these OGs
  samples_tiss <- coldata$Sample[coldata$Tissue == tiss]
  if (length(samples_tiss) == 0) {
    warning("  No samples for tissue ", tiss, " in coldata.")
    return(invisible(NULL))
  }
  
  samples_tiss <- intersect(samples_tiss, colnames(vsd_mat))
  if (length(samples_tiss) == 0) {
    warning("  No matching samples for tissue ", tiss, " in vsd_mat.")
    return(invisible(NULL))
  }
  
  og_ids <- intersect(og_ids, rownames(vsd_mat))
  if (length(og_ids) == 0) {
    warning("  No OGs from module ", module_color, " found in vsd_mat for ", tiss)
    return(invisible(NULL))
  }
  
  expr_sub <- vsd_mat[og_ids, samples_tiss, drop = FALSE]
  
  # 2.4 Row labels: use gene_name if available, otherwise OG
  gene_names_vec <- gi_mod$gene_name[match(og_ids, gi_mod$OG)]
  row_labels <- ifelse(!is.na(gene_names_vec) & gene_names_vec != "",
                       gene_names_vec, og_ids)
  dup <- duplicated(row_labels)
  row_labels[dup] <- paste0(row_labels[dup], " (", og_ids[dup], ")")
  rownames(expr_sub) <- row_labels
  
  # 2.5 Column annotation: Tissue, Species, Dry vs Mesic
  ann_col <- coldata %>%
    dplyr::filter(Sample %in% samples_tiss) %>%
    dplyr::select(Sample, Species, Tissue)
  ann_col <- ann_col[match(colnames(expr_sub), ann_col$Sample), ]
  
  ann_col$DryMesic <- ifelse(ann_col$Species %in% dry_species, "Dry", "Mesic")
  rownames(ann_col) <- ann_col$Sample
  ann_col$Sample <- NULL
  
  ann_colors <- list(
    Tissue = setNames(RColorBrewer::brewer.pal(length(tissues_to_run), "Set3"),
                      tissues_to_run),
    DryMesic = c(Dry = "#E74504", Mesic = "#585757")
  )
  
  # 2.6 Plot to RStudio plot pane
  res <- pheatmap(expr_sub,
                  scale         = "row",
                  clustering_method = "complete",
                  clustering_distance_rows = "euclidean",
                  clustering_distance_cols = "euclidean",
                  annotation_col = ann_col,
                  annotation_colors = ann_colors,
                  show_rownames  = TRUE,
                  show_colnames  = FALSE,
                  fontsize_row   = 6,
                  fontsize_col   = 8,
                  color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
                  main  = paste0("Tissue: ", tiss,
                                 " – module ", module_color,
                                 " (top genes, vst)"))
  
  # 2.7 Optionally save to PNG
  if (save_png) {
    out_png <- paste0("Heatmap_WGCNA_", tiss, "_", module_color, "_vst.png")
    png(filename = out_png, width = 7, height = 8, units = "in", res = 300)
    pheatmap(expr_sub,
             scale         = "row",
             clustering_method = "complete",
             clustering_distance_rows = "euclidean",
             clustering_distance_cols = "euclidean",
             annotation_col = ann_col,
             annotation_colors = ann_colors,
             show_rownames  = TRUE,
             show_colnames  = FALSE,
             fontsize_row   = 6,
             fontsize_col   = 8,
             color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
             main  = paste0("Tissue: ", tiss,
                            " – module ", module_color,
                            " (top genes, vst)"))
    dev.off()
    message("  Wrote heatmap: ", out_png)
  }
  
  invisible(res)
}

############################################################
## 3) Loop over tissues and top 3 modules
############################################################

for (tiss in tissues_to_run) {
  top3_tiss <- subset(top3_all, Tissue == tiss)
  if (nrow(top3_tiss) == 0) {
    warning("No top3 modules recorded for tissue ", tiss)
    next
  }
  
  for (i in seq_len(nrow(top3_tiss))) {
    color <- top3_tiss$Color[i]  # module color name, e.g. "turquoise"
    plot_module_heatmap(tiss, color, vsd_mat, coldata, max_genes = 100)
  }
}

plot_module_heatmap("Brain", "yellow", vsd_mat, coldata, max_genes = 100, save_png = FALSE)
