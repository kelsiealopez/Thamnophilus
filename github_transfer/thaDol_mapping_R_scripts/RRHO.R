setwd("/n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/nf_pipeline/thaDol_mapping_test")

############################################################
## RRHO analysis for convergent Dry vs Mesic expression
## across species pairs and tissues
##
## Assumes you have DE tables:
##   DE_OG_selected_<spA>_vs_<spB>_<Tissue>_kallisto.csv
## produced by your DE script.
############################################################

library(RRHO)         # Bioconductor package
library(dplyr)
library(tidyr)
library(ggplot2)

##----------------------------------------------------------
## 0) Basic config
##----------------------------------------------------------

############################################################
## 1) Pairs, tissues, output directory
############################################################

pairs_of_interest <- tribble(
  ~pair_id,               ~spA,     ~spB,
  "dysSti_vs_sakCri",     "dysSti", "sakCri",
  "sakLuc_vs_sakCan",     "sakLuc", "sakCan",
  "thaAmb_vs_thaPel",     "thaAmb", "thaPel",
  "thaBer_vs_thaAtr",     "thaAtr", "thaBer",
  "thaRuf_vs_thaTor",     "thaRuf", "thaTor"
)

tissues_to_run <- c("Brain", "Heart", "Liver", "Muscle", "Kidney")

outdir_rrho <- "RRHO_results"
if (!dir.exists(outdir_rrho)) dir.create(outdir_rrho, recursive = TRUE)

############################################################
## 2) Builder: rank vector with deduplicated OGs
############################################################

build_rrho_rank_vector <- function(tiss, spA, spB,
                                   prefix = "DE_OG_selected",
                                   suffix = "_kallisto.csv") {
  fname <- paste0(prefix, "_", spA, "_vs_", spB, "_", tiss, suffix)
  if (!file.exists(fname)) {
    warning("DE file not found for ", tiss, " ", spA, " vs ", spB, ": ", fname)
    return(NULL)
  }
  
  de_tab <- read.csv(fname, row.names = 1, check.names = FALSE)
  de_tab$OG <- rownames(de_tab)
  
  if (!all(c("pvalue", "log2FoldChange") %in% colnames(de_tab))) {
    stop("DE table ", fname, " missing pvalue or log2FoldChange columns.")
  }
  
  de_sub <- de_tab %>%
    dplyr::filter(!is.na(pvalue), !is.na(log2FoldChange))
  
  if (nrow(de_sub) == 0) {
    warning("No usable rows in ", fname)
    return(NULL)
  }
  
  # signed score
  eps <- 1e-300
  score <- sign(de_sub$log2FoldChange) * -log10(de_sub$pvalue + eps)
  og_id <- de_sub$OG
  
  # collapse duplicates by max |score|
  score_df <- data.frame(OG = og_id, score = score, stringsAsFactors = FALSE)
  score_collapsed <- tapply(score_df$score, score_df$OG,
                            function(x) x[which.max(abs(x))])
  
  score_vec <- as.numeric(score_collapsed)
  names(score_vec) <- names(score_collapsed)
  
  # sort
  score_vec <- sort(score_vec, decreasing = TRUE)
  score_vec
}

############################################################
## 3) Build rrho_input with deduplicated OGs
############################################################

rrho_input <- list()

for (tiss in tissues_to_run) {
  message("\n=== Building RRHO rank vectors for tissue: ", tiss, " ===")
  rrho_input[[tiss]] <- list()
  
  for (i in seq_len(nrow(pairs_of_interest))) {
    pid <- pairs_of_interest$pair_id[i]
    spA <- pairs_of_interest$spA[i]
    spB <- pairs_of_interest$spB[i]
    
    v <- build_rrho_rank_vector(tiss, spA, spB)
    if (is.null(v)) {
      message("  Skipping ", pid, " in ", tiss, " (no vector).")
      next
    }
    
    rrho_input[[tiss]][[pid]] <- v
    message("  Built rank vector for ", pid, " (", length(v), " unique OGs).")
  }
}

############################################################
## 4) RRHO helper: 2-column data.frames, intersect gene sets
############################################################

run_rrho_pair <- function(vec1, vec2,
                          stepsize = 100,
                          alternative = "enrichment",
                          seed = 123) {
  set.seed(seed)
  
  if (is.null(names(vec1)) || is.null(names(vec2))) {
    stop("Both input vectors must be named (gene IDs as names).")
  }
  
  # Only use genes present in both lists
  common_genes <- intersect(names(vec1), names(vec2))
  if (length(common_genes) < 10) {
    stop("Too few common genes between the two vectors.")
  }
  
  v1 <- vec1[common_genes]
  v2 <- vec2[common_genes]
  
  # 2-column data.frames: gene ID + score
  df1 <- data.frame(gene = common_genes,
                    score = as.numeric(v1),
                    stringsAsFactors = FALSE)
  df2 <- data.frame(gene = common_genes,
                    score = as.numeric(v2),
                    stringsAsFactors = FALSE)
  
  # RRHO expects first column = ID, second column = score
  RRHO(df1, df2,
       stepsize    = stepsize,
       alternative = alternative,
       log10.ind   = TRUE,
       plot        = FALSE)
}

############################################################
## 5) Run empirical RRHO
############################################################

rrho_empirical <- list()

for (tiss in tissues_to_run) {
  vecs <- rrho_input[[tiss]]
  if (is.null(vecs) || length(vecs) < 2) {
    message("Not enough pairs for RRHO in tissue ", tiss)
    next
  }
  
  pairs_ids <- names(vecs)
  combs     <- combn(pairs_ids, 2, simplify = FALSE)
  
  rrho_empirical[[tiss]] <- list()
  
  for (cmb in combs) {
    pid1 <- cmb[1]
    pid2 <- cmb[2]
    v1   <- vecs[[pid1]]
    v2   <- vecs[[pid2]]
    
    message("\nTissue ", tiss, " – RRHO for ", pid1, " vs ", pid2)
    rr <- run_rrho_pair(v1, v2, stepsize = 100, alternative = "enrichment")
    
    key <- paste(pid1, "vs", pid2, sep = "_")
    rrho_empirical[[tiss]][[key]] <- rr
    
    # Save object
    saveRDS(rr,
            file = file.path(outdir_rrho,
                             paste0("RRHO_", tiss, "_", pid1, "_vs_", pid2, "_empirical.rds")))
    
    # Heatmap
    png(file.path(outdir_rrho,
                  paste0("RRHO_", tiss, "_", pid1, "_vs_", pid2, "_empirical.png")),
        width = 800, height = 800)
    try({
      plot(rr,
           main = paste("RRHO –", tiss, "\n", pid1, "vs", pid2, " (empirical)"))
    }, silent = TRUE)
    dev.off()
  }
}

message("\nEmpirical RRHO complete. Results in ", outdir_rrho, "\n")

############################################################
## 5) RRHO null: permuted log2FC/pvalue within each pair
############################################################

build_rrho_rank_vector_permuted <- function(tiss, spA, spB,
                                            prefix = "DE_OG_selected",
                                            suffix = "_kallisto.csv") {
  fname <- paste0(prefix, "_", spA, "_vs_", spB, "_", tiss, suffix)
  if (!file.exists(fname)) {
    warning("DE file not found for ", tiss, " ", spA, " vs ", spB, ": ", fname)
    return(NULL)
  }
  
  de_tab <- read.csv(fname, row.names = 1, check.names = FALSE)
  de_tab$OG <- rownames(de_tab)
  
  if (!all(c("pvalue", "log2FoldChange") %in% colnames(de_tab))) {
    stop("DE table ", fname, " missing pvalue or log2FoldChange columns.")
  }
  
  de_sub <- de_tab %>%
    filter(!is.na(pvalue), !is.na(log2FoldChange))
  
  if (nrow(de_sub) == 0) {
    warning("No usable rows in ", fname)
    return(NULL)
  }
  
  # Permute the (log2FC, pvalue) association across genes
  idx_perm <- sample(seq_len(nrow(de_sub)))
  de_perm  <- de_sub
  de_perm$log2FoldChange <- de_sub$log2FoldChange[idx_perm]
  de_perm$pvalue         <- de_sub$pvalue[idx_perm]
  
  eps   <- 1e-300
  score <- sign(de_perm$log2FoldChange) * -log10(de_perm$pvalue + eps)
  names(score) <- de_perm$OG
  score <- sort(score, decreasing = TRUE)
  
  score
}

# rrho_input_null[[tissue]][[pair_id]] = permuted rank vector
rrho_input_null <- list()

for (tiss in tissues_to_run) {
  message("\n=== Building RRHO *null* rank vectors for tissue: ", tiss, " ===")
  rrho_input_null[[tiss]] <- list()
  
  for (i in seq_len(nrow(pairs_of_interest))) {
    pid <- pairs_of_interest$pair_id[i]
    spA <- pairs_of_interest$spA[i]
    spB <- pairs_of_interest$spB[i]
    
    v <- build_rrho_rank_vector_permuted(tiss, spA, spB)
    if (is.null(v)) {
      message("  Skipping ", pid, " in ", tiss, " (no null vector).")
      next
    }
    
    rrho_input_null[[tiss]][[pid]] <- v
    message("  Built permuted rank vector for ", pid, " (", length(v), " genes).")
  }
}

# RRHO for null
rrho_null <- list()

for (tiss in tissues_to_run) {
  vecs_null <- rrho_input_null[[tiss]]
  if (is.null(vecs_null) || length(vecs_null) < 2) {
    message("Not enough pairs for null RRHO in tissue ", tiss)
    next
  }
  
  pairs_ids <- names(vecs_null)
  combs     <- combn(pairs_ids, 2, simplify = FALSE)
  
  rrho_null[[tiss]] <- list()
  
  for (cmb in combs) {
    pid1 <- cmb[1]
    pid2 <- cmb[2]
    v1   <- vecs_null[[pid1]]
    v2   <- vecs_null[[pid2]]
    
    message("\nTissue ", tiss, " – RRHO NULL for ", pid1, " vs ", pid2)
    rr <- run_rrho_pair(v1, v2, stepsize = 100, alternative = "enrichment")
    
    key <- paste(pid1, "vs", pid2, sep = "_")
    rrho_null[[tiss]][[key]] <- rr
    
    # Save RRHO object
    saveRDS(rr,
            file = file.path(outdir_rrho,
                             paste0("RRHO_", tiss, "_", pid1, "_vs_", pid2, "_null.rds")))
    
    # Plot heatmap
    png(file.path(outdir_rrho,
                  paste0("RRHO_", tiss, "_", pid1, "_vs_", pid2, "_null.png")),
        width = 800, height = 800)
    try({
      plot(rr,
           main = paste("RRHO –", tiss, "\n", pid1, "vs", pid2, "(null permuted)"))
    }, silent = TRUE)
    dev.off()
  }
}

############################################################
## 6) Quick helper: extract strongly concordant genes
############################################################

extract_rrho_concordant_genes <- function(rrho_obj, cutoff.quantile = 0.99) {
  hm <- rrho_obj$hypermat
  
  # Threshold
  thr <- quantile(hm, probs = cutoff.quantile, na.rm = TRUE)
  
  idx <- which(hm >= thr, arr.ind = TRUE)
  if (nrow(idx) == 0) {
    message("No cells above the chosen quantile; returning empty set.")
    return(character(0))
  }
  
  ranks1 <- rrho_obj$dim1$dim[idx[, 1]]
  ranks2 <- rrho_obj$dim2$dim[idx[, 2]]
  
  # Remove NAs
  ranks1 <- ranks1[is.finite(ranks1)]
  ranks2 <- ranks2[is.finite(ranks2)]
  if (length(ranks1) == 0 || length(ranks2) == 0) {
    message("No finite ranks; returning empty set.")
    return(character(0))
  }
  
  genes1 <- names(rrho_obj$list1)[seq_len(max(ranks1))]
  genes2 <- names(rrho_obj$list2)[seq_len(max(ranks2))]
  
  intersect(genes1, genes2)
}

rr_example <- readRDS("RRHO_results/RRHO_Brain_dysSti_vs_sakCri_vs_sakLuc_vs_sakCan_empirical.rds")
concordant_genes <- extract_rrho_concordant_genes(rr_example, cutoff.quantile = 0.95)
length(concordant_genes)
head(concordant_genes)












plot_rrho_heatmap <- function(rr_obj,
                              main = "RRHO",
                              col = colorRampPalette(c("white", "yellow", "red"))(256),
                              zlim = NULL) {
  hm <- rr_obj$hypermat
  
  if (all(is.na(hm))) {
    warning("hypermat is all NA; nothing to plot.")
    return(invisible(NULL))
  }
  
  # Flip vertically so origin matches the usual RRHO orientation
  hm_plot <- t(hm[nrow(hm):1, , drop = FALSE])
  
  if (is.null(zlim)) {
    zlim <- range(hm, finite = TRUE)
  }
  
  image(hm_plot,
        col = col,
        zlim = zlim,
        axes = FALSE,
        xlab = "List 1 rank (increasing → right)",
        ylab = "List 2 rank (increasing → up)",
        main = main)
  box()
}


library(RRHO)

v1 <- rrho_input[["Brain"]][["dysSti_vs_sakCri"]]
v2 <- rrho_input[["Brain"]][["sakLuc_vs_sakCan"]]
rr_test <- run_rrho_pair(v1, v2, stepsize = 100, alternative = "enrichment")

# Plot to screen
plot_rrho_heatmap(rr_test, main = "Test RRHO – Brain dysSti_vs_sakCri vs sakLuc_vs_sakCan")









# plot all !

# make sure plot_rrho_heatmap() is defined and RRHO is loaded
# source("your_script_with_rrho_input_and_rrho_empirical.R") if needed

plot_all_empirical_rrho <- function(rrho_empirical_list) {
  for (tiss in names(rrho_empirical_list)) {
    cat("\n=== Tissue:", tiss, "===\n")
    
    rr_list <- rrho_empirical_list[[tiss]]
    if (is.null(rr_list) || length(rr_list) == 0) next
    
    for (key in names(rr_list)) {
      rr <- rr_list[[key]]
      
      # Make a nice title, e.g. "dysSti_vs_sakCri vs sakLuc_vs_sakCan"
      title_main <- paste("RRHO –", tiss, "\n", gsub("_vs_", " vs ", key))
      
      # open new plot
      dev.new()  # or windows()/x11() depending on your system, but dev.new() is portable
      plot_rrho_heatmap(rr, main = title_main)
      
      cat("Plotted:", tiss, "-", key, "\n")
      readline("Press [Enter] for next plot...")
      dev.off()
    }
  }
}

# run it
plot_all_empirical_rrho(rrho_empirical)






#
#
#
#

# custom colors
rrho_colors <- c(
  "#EEE3D4",
  "#DF9E7B",
  "#DB7658",
  "#D85E4B",
  "#CD454D",
  "#BC3955",
  "#A6335E",
  "#912E60",
  "#752A60",
  "#62265D",
  "#48214E",
  "#341C42",
  "#0D0E21"
)

rrho_palette <- colorRampPalette(rrho_colors)


plot_rrho_heatmap <- function(rr_obj,
                              main = "RRHO",
                              ncol = 256,
                              zlim = NULL) {
  hm <- rr_obj$hypermat
  
  if (all(is.na(hm))) {
    warning("hypermat is all NA; nothing to plot.")
    return(invisible(NULL))
  }
  
  # flip vertically so origin matches usual RRHO orientation
  hm_plot <- t(hm[nrow(hm):1, , drop = FALSE])
  
  # RRHO values are -log10(p), so they should be ≥ 0
  if (is.null(zlim)) {
    zlim <- c(0, max(hm, na.rm = TRUE))
  }
  
  image(hm_plot,
        col = rrho_palette(ncol),
        zlim = zlim,
        axes = FALSE,
        xlab = "List 1 rank (increasing → right)",
        ylab = "List 2 rank (increasing → up)",
        main = main)
  box()
}


plot_rrho_heatmap(rr_test, main = "Test RRHO – Brain dysSti_vs_sakCri vs sakLuc_vs_sakCan")




plot_all_empirical_rrho <- function(rrho_empirical_list) {
  for (tiss in names(rrho_empirical_list)) {
    cat("\n=== Tissue:", tiss, "===\n")
    
    rr_list <- rrho_empirical_list[[tiss]]
    if (is.null(rr_list) || length(rr_list) == 0) next
    
    for (key in names(rr_list)) {
      rr <- rr_list[[key]]
      
      # Title like "RRHO – Brain\n dysSti_vs_sakCri vs sakLuc_vs_sakCan"
      title_main <- paste("RRHO –", tiss, "\n", gsub("_vs_", " vs ", key))
      
      # New device for each plot (you can skip dev.new() if you’re fine reusing the same window)
      dev.new()
      plot_rrho_heatmap(rr, main = title_main)
      
      cat("Plotted:", tiss, "-", key, "\n")
      readline("Press [Enter] for next plot...")
      dev.off()
    }
  }
}

# Run it
plot_all_empirical_rrho(rrho_empirical)




plot_empirical_for_tissue <- function(rrho_empirical_list, tissue, nrow = 2, ncol = 2) {
  rr_list <- rrho_empirical_list[[tissue]]
  if (is.null(rr_list) || length(rr_list) == 0) {
    message("No RRHO empirical objects for tissue ", tissue)
    return(NULL)
  }
  
  comps <- names(rr_list)
  n <- length(comps)
  
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  
  par(mfrow = c(nrow, ncol))
  
  for (i in seq_along(comps)) {
    key <- comps[i]
    rr  <- rr_list[[key]]
    title_main <- paste("RRHO –", tissue, "\n", gsub("_vs_", " vs ", key))
    plot_rrho_heatmap(rr, main = title_main)
    
    if (i %% (nrow * ncol) == 0 && i != n) {
      readline("Press [Enter] for next page of plots...")
      par(mfrow = c(nrow, ncol))
    }
  }
}

# Example: Brain only
plot_empirical_for_tissue(rrho_empirical, "Brain", nrow = 2, ncol = 2)





tab <- read.csv("DE_OG_selected_thaBer_vs_thaAtr_Brain_kallisto.csv",
                row.names = 1, check.names = FALSE)

# Check range
summary(tab$log2FoldChange)[1:6]

# Maybe look at a few extreme genes
head(tab[order(-tab$log2FoldChange), ])  # top positive LFC
head(tab[order(tab$log2FoldChange), ])   # top negative LFC
