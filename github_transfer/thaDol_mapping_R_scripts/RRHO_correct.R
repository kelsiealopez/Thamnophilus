setwd("/n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/nf_pipeline/thaDol_mapping_test")

############################################################
## RRHO analysis for convergent Dry vs Mesic expression
##
## Requires DE tables named:
##   DE_OG_selected_<spA>_vs_<spB>_<Tissue>_kallisto.csv
############################################################

library(RRHO)
library(dplyr)
library(tidyr)
library(ggplot2)

##----------------------------------------------------------
## 1) Pairs, tissues, output directory
##----------------------------------------------------------

## IMPORTANT:
## spA, spB MUST match the order in the DE filenames.
## dry_species = which of those species is the dry (cavity) species.
## We will flip log2FC where needed so that positive = up in dry_species.
pairs_of_interest <- tribble(
  ~pair_id,           ~spA,     ~spB,     ~dry_species,
  # humid vs dry
  "dysSti_vs_sakCri", "dysSti", "sakCri", "sakCri",  # dry = sakCri
  "sakLuc_vs_sakCan", "sakLuc", "sakCan", "sakCan",  # dry = sakCan
  "thaAmb_vs_thaPel", "thaAmb", "thaPel", "thaPel",  # dry = thaPel
  # filename is DE_OG_selected_thaBer_vs_thaAtr_...
  "thaBer_vs_thaAtr", "thaBer", "thaAtr", "thaBer",  # dry = thaBer
  "thaRuf_vs_thaTor", "thaRuf", "thaTor", "thaTor"   # dry = thaTor
)

tissues_to_run <- c("Brain", "Heart", "Liver", "Muscle", "Kidney")

outdir_rrho <- "RRHO_results"
if (!dir.exists(outdir_rrho)) dir.create(outdir_rrho, recursive = TRUE)

############################################################
## 2) Helper: build rank vector with deduplicated OGs
##    enforcing positive score = up in dry species
############################################################

## ASSUMPTION (confirm from DE code):
##   log2FoldChange in DE_OG_selected_<spA>_vs_<spB>_* is (spB - spA).
## Then:
##   - if dry_species == spB: positive LFC already = up in dry
##   - if dry_species == spA: flip sign so positive = up in dry
build_rrho_rank_vector <- function(tiss, spA, spB, dry_species,
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
  
  # adjust log2FC so positive always = up in dry species
  log2fc <- de_sub$log2FoldChange
  if (dry_species == spA) {
    # current LFC is spB - spA; flip so it becomes spA - spB (dry - humid)
    log2fc <- -log2fc
  }
  # if dry_species == spB, leave log2fc as is
  
  # signed score: direction (dry up vs humid up) * significance
  eps   <- 1e-300
  score <- sign(log2fc) * -log10(de_sub$pvalue + eps)
  og_id <- de_sub$OG
  
  # collapse duplicates by max |score|
  score_df <- data.frame(OG = og_id, score = score, stringsAsFactors = FALSE)
  score_collapsed <- tapply(score_df$score, score_df$OG,
                            function(x) x[which.max(abs(x))])
  
  score_vec <- as.numeric(score_collapsed)
  names(score_vec) <- names(score_collapsed)
  
  # sort: largest positive first (strong up in dry)
  sort(score_vec, decreasing = TRUE)
}

############################################################
## 3) Build rrho_input for all tissues/pairs
############################################################

rrho_input <- list()

for (tiss in tissues_to_run) {
  message("\n=== Building RRHO rank vectors for tissue: ", tiss, " ===")
  rrho_input[[tiss]] <- list()
  
  for (i in seq_len(nrow(pairs_of_interest))) {
    pid  <- pairs_of_interest$pair_id[i]
    spA  <- pairs_of_interest$spA[i]
    spB  <- pairs_of_interest$spB[i]
    dry  <- pairs_of_interest$dry_species[i]
    
    v <- build_rrho_rank_vector(tiss, spA, spB, dry)
    if (is.null(v)) {
      message("  Skipping ", pid, " in ", tiss, " (no vector).")
      next
    }
    
    rrho_input[[tiss]][[pid]] <- v
    message("  Built rank vector for ", pid, " (", length(v), " unique OGs).")
  }
}

############################################################
## 4) RRHO helper to run pairwise RRHO
############################################################

run_rrho_pair <- function(vec1, vec2,
                          stepsize = 100,
                          alternative = "enrichment",
                          seed = 123) {
  set.seed(seed)
  
  if (is.null(names(vec1)) || is.null(names(vec2))) {
    stop("Both input vectors must be named (gene IDs as names).")
  }
  
  common_genes <- intersect(names(vec1), names(vec2))
  if (length(common_genes) < 10) {
    stop("Too few common genes between the two vectors.")
  }
  
  v1 <- vec1[common_genes]
  v2 <- vec2[common_genes]
  
  df1 <- data.frame(gene = common_genes,
                    score = as.numeric(v1),
                    stringsAsFactors = FALSE)
  df2 <- data.frame(gene = common_genes,
                    score = as.numeric(v2),
                    stringsAsFactors = FALSE)
  
  RRHO(df1, df2,
       stepsize    = stepsize,
       alternative = alternative,
       log10.ind   = TRUE,
       plot        = FALSE)
}

############################################################
## 5) Run empirical RRHO and save objects
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
    
    saveRDS(rr,
            file = file.path(outdir_rrho,
                             paste0("RRHO_", tiss, "_", pid1, "_vs_", pid2, "_empirical.rds")))
  }
}

message("\nEmpirical RRHO complete. Results in ", outdir_rrho, "\n")

############################################################
## 6) Build permuted (null) rank vectors with same dry/humid convention
############################################################

build_rrho_rank_vector_permuted <- function(tiss, spA, spB, dry_species,
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
  
  # Now enforce dry vs humid direction as above
  log2fc <- de_perm$log2FoldChange
  if (dry_species == spA) {
    log2fc <- -log2fc
  }
  
  eps   <- 1e-300
  score <- sign(log2fc) * -log10(de_perm$pvalue + eps)
  names(score) <- de_perm$OG
  
  sort(score, decreasing = TRUE)
}

rrho_input_null <- list()

for (tiss in tissues_to_run) {
  message("\n=== Building RRHO *null* rank vectors for tissue: ", tiss, " ===")
  rrho_input_null[[tiss]] <- list()
  
  for (i in seq_len(nrow(pairs_of_interest))) {
    pid <- pairs_of_interest$pair_id[i]
    spA <- pairs_of_interest$spA[i]
    spB <- pairs_of_interest$spB[i]
    dry <- pairs_of_interest$dry_species[i]
    
    v <- build_rrho_rank_vector_permuted(tiss, spA, spB, dry)
    if (is.null(v)) {
      message("  Skipping ", pid, " in ", tiss, " (no null vector).")
      next
    }
    
    rrho_input_null[[tiss]][[pid]] <- v
    message("  Built permuted rank vector for ", pid, " (", length(v), " genes).")
  }
}

############################################################
## 7) RRHO for null
############################################################

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
    
    saveRDS(rr,
            file = file.path(outdir_rrho,
                             paste0("RRHO_", tiss, "_", pid1, "_vs_", pid2, "_null.rds")))
  }
}

message("\nRRHO (empirical + null) complete. Results in ", outdir_rrho, "\n")

############################################################
## 8) Pretty plotting helpers (one plot at a time)
############################################################

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

## Plot ALL empirical RRHO maps, one at a time
plot_all_empirical_rrho <- function(rrho_empirical_list) {
  for (tiss in names(rrho_empirical_list)) {
    cat("\n=== Tissue:", tiss, "===\n")
    
    rr_list <- rrho_empirical_list[[tiss]]
    if (is.null(rr_list) || length(rr_list) == 0) next
    
    for (key in names(rr_list)) {
      rr <- rr_list[[key]]
      
      title_main <- paste("RRHO –", tiss, "\n", gsub("_vs_", " vs ", key))
      
      dev.new()
      plot_rrho_heatmap(rr, main = title_main)
      
      cat("Plotted empirical:", tiss, "-", key, "\n")
      readline("Press [Enter] for next plot...")
      dev.off()
    }
  }
}

## Plot ALL null RRHO maps, one at a time (optional)
plot_all_null_rrho <- function(rrho_null_list) {
  for (tiss in names(rrho_null_list)) {
    cat("\n=== Tissue (NULL):", tiss, "===\n")
    
    rr_list <- rrho_null_list[[tiss]]
    if (is.null(rr_list) || length(rr_list) == 0) next
    
    for (key in names(rr_list)) {
      rr <- rr_list[[key]]
      
      title_main <- paste("RRHO NULL –", tiss, "\n", gsub("_vs_", " vs ", key))
      
      dev.new()
      plot_rrho_heatmap(rr, main = title_main)
      
      cat("Plotted null:", tiss, "-", key, "\n")
      readline("Press [Enter] for next null plot...")
      dev.off()
    }
  }
}

## To actually plot everything (empirical only):
plot_all_empirical_rrho(rrho_empirical)

## Or to also look at nulls:
plot_all_null_rrho(rrho_null)






#trying different plotting commands to increasse contrast

plot_rrho_heatmap <- function(rr_obj,
                              main = "RRHO",
                              ncol = 256,
                              zlim = NULL,
                              vmax_quantile = 0.99) {
  hm <- rr_obj$hypermat
  
  if (all(is.na(hm))) {
    warning("hypermat is all NA; nothing to plot.")
    return(invisible(NULL))
  }
  
  hm_plot <- t(hm[nrow(hm):1, , drop = FALSE])
  
  if (is.null(zlim)) {
    # lower bound 0; upper bound = high quantile of hm
    upper <- as.numeric(quantile(hm, probs = vmax_quantile, na.rm = TRUE))
    if (!is.finite(upper) || upper <= 0) {
      upper <- max(hm, na.rm = TRUE)
    }
    zlim <- c(0, upper)
  }
  
  # clip values above zlim[2] so they all get max color
  hm_plot_clipped <- pmin(hm_plot, zlim[2])
  
  image(hm_plot_clipped,
        col  = rrho_palette(ncol),
        zlim = zlim,
        axes = FALSE,
        xlab = "List 1 rank (increasing → right)",
        ylab = "List 2 rank (increasing → up)",
        main = main)
  box()
}

plot_all_empirical_rrho(rrho_empirical)











# try rainbow palette

rrho_rainbow_palette <- colorRampPalette(c(
  "#00004E",  # dark blue
  "#0000FF",  # blue
  "#00FFFF",  # cyan
  "#00FF00",  # green
  "#FFFF00",  # yellow
  "#FFA500",  # orange
  "#FF0000",  # red
  "#800000"   # dark red
))

plot_rrho_heatmap <- function(rr_obj,
                              main = "RRHO",
                              ncol = 256,
                              zlim = NULL,
                              vmax_quantile = 0.99,
                              palette = rrho_palette) {
  hm <- rr_obj$hypermat
  if (all(is.na(hm))) {
    warning("hypermat is all NA; nothing to plot.")
    return(invisible(NULL))
  }
  
  hm_plot <- t(hm[nrow(hm):1, , drop = FALSE])
  
  if (is.null(zlim)) {
    upper <- as.numeric(quantile(hm, probs = vmax_quantile, na.rm = TRUE))
    if (!is.finite(upper) || upper <= 0) {
      upper <- max(hm, na.rm = TRUE)
    }
    zlim <- c(0, upper)
  }
  hm_plot_clipped <- pmin(hm_plot, zlim[2])
  
  image(hm_plot_clipped,
        col  = palette(ncol),
        zlim = zlim,
        axes = FALSE,
        xlab = "List 1 rank (increasing → right)",
        ylab = "List 2 rank (increasing → up)",
        main = main)
  box()
}

# pretty purple palette (your current one)
plot_rrho_heatmap(rr, main = "Empirical", palette = rrho_palette)

# rainbow-style palette
plot_rrho_heatmap(rr, main = "Empirical (rainbow)", palette = rrho_rainbow_palette)







## Plot ALL empirical RRHO maps, one at a time, with chosen palette
plot_all_empirical_rrho <- function(rrho_empirical_list,
                                    palette = rrho_palette,
                                    vmax_quantile = 0.99) {
  for (tiss in names(rrho_empirical_list)) {
    cat("\n=== Tissue:", tiss, "===\n")
    
    rr_list <- rrho_empirical_list[[tiss]]
    if (is.null(rr_list) || length(rr_list) == 0) next
    
    for (key in names(rr_list)) {
      rr <- rr_list[[key]]
      
      title_main <- paste("RRHO –", tiss, "\n", gsub("_vs_", " vs ", key))
      
      dev.new()
      plot_rrho_heatmap(rr,
                        main          = title_main,
                        palette       = palette,
                        vmax_quantile = vmax_quantile)
      
      cat("Plotted empirical:", tiss, "-", key, "\n")
      readline("Press [Enter] for next plot...")
      dev.off()
    }
  }
}


plot_all_empirical_rrho(rrho_empirical, palette = rrho_rainbow_palette, vmax_quantile = 0.99)





#
#
#

# add legend

plot_rrho_heatmap <- function(rr_obj,
                              main = "RRHO",
                              ncol = 256,
                              zlim = NULL,
                              vmax_quantile = 0.99,
                              palette = rrho_palette) {
  hm <- rr_obj$hypermat
  if (all(is.na(hm))) {
    warning("hypermat is all NA; nothing to plot.")
    return(invisible(NULL))
  }
  
  # Flip vertically so origin matches usual RRHO orientation
  hm_plot <- t(hm[nrow(hm):1, , drop = FALSE])
  
  # Determine zlim using a high quantile for contrast
  if (is.null(zlim)) {
    upper <- as.numeric(quantile(hm, probs = vmax_quantile, na.rm = TRUE))
    if (!is.finite(upper) || upper <= 0) {
      upper <- max(hm, na.rm = TRUE)
    }
    zlim <- c(0, upper)
  }
  
  # Clip values above zlim[2] so they all get max color
  hm_plot_clipped <- pmin(hm_plot, zlim[2])
  
  # Split device: 1 = heatmap, 2 = color bar
  layout(matrix(c(1, 2), ncol = 2), widths = c(4, 1))
  
  ## 1) Heatmap
  par(mar = c(5, 5, 4, 1))  # bottom, left, top, right
  image(hm_plot_clipped,
        col  = palette(ncol),
        zlim = zlim,
        axes = FALSE,
        xlab = "List 1 rank (increasing → right)",
        ylab = "List 2 rank (increasing → up)",
        main = main)
  box()
  
  ## 2) Color bar
  par(mar = c(5, 1, 4, 3))
  # simple 0–1 vertical gradient
  grad <- matrix(seq(0, 1, length.out = ncol), nrow = 1)
  image(1, seq(0, 1, length.out = ncol), grad,
        col  = palette(ncol),
        xlab = "",
        ylab = expression(-log[10](p)),
        axes = FALSE)
  # map 0–1 positions to zlim values on the axis
  at_pos  <- seq(0, 1, length.out = 5)
  at_vals <- seq(zlim[1], zlim[2], length.out = 5)
  axis(4, at = at_pos, labels = signif(at_vals, 2))
  box()
  
  # Reset layout
  layout(1)
}



plot_all_empirical_rrho(rrho_empirical,
                        palette       = rrho_rainbow_palette,
                        vmax_quantile = 0.99)




# save as pngs

## Save ALL empirical RRHO maps to PNGs (one file per comparison)
save_all_empirical_rrho_png <- function(rrho_empirical_list,
                                        outdir = "RRHO_results",
                                        palette = rrho_rainbow_palette,
                                        vmax_quantile = 0.99,
                                        prefix = "RRHO") {
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  
  for (tiss in names(rrho_empirical_list)) {
    cat("\n=== Tissue (saving PNGs):", tiss, "===\n")
    
    rr_list <- rrho_empirical_list[[tiss]]
    if (is.null(rr_list) || length(rr_list) == 0) next
    
    for (key in names(rr_list)) {
      rr <- rr_list[[key]]
      
      # key is like "dysSti_vs_sakCri_vs_sakLuc_vs_sakCan"
      # construct a filename; add suffix to indicate rainbow palette
      fname <- file.path(
        outdir,
        paste0(prefix, "_", tiss, "_", key, "_empirical_rainbow.png")
      )
      
      title_main <- paste("RRHO –", tiss, "\n", gsub("_vs_", " vs ", key))
      
      png(fname, width = 900, height = 800)
      plot_rrho_heatmap(rr,
                        main          = title_main,
                        palette       = palette,
                        vmax_quantile = vmax_quantile)
      dev.off()
      
      cat("Saved empirical PNG:", fname, "\n")
    }
  }
}

save_all_empirical_rrho_png(rrho_empirical,
                            outdir        = outdir_rrho,
                            palette       = rrho_rainbow_palette,
                            vmax_quantile = 0.99)











# also wnat to save the null plots !

## Save ALL null RRHO maps to PNGs (one file per comparison)
save_all_null_rrho_png <- function(rrho_null_list,
                                   outdir = "RRHO_results",
                                   palette = rrho_rainbow_palette,
                                   vmax_quantile = 0.99,
                                   prefix = "RRHO") {
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  
  for (tiss in names(rrho_null_list)) {
    cat("\n=== Tissue (saving NULL PNGs):", tiss, "===\n")
    
    rr_list <- rrho_null_list[[tiss]]
    if (is.null(rr_list) || length(rr_list) == 0) next
    
    for (key in names(rr_list)) {
      rr <- rr_list[[key]]
      
      # key is like "dysSti_vs_sakCri_vs_sakLuc_vs_sakCan"
      # construct filename; mark as null + rainbow
      fname <- file.path(
        outdir,
        paste0(prefix, "_", tiss, "_", key, "_null_rainbow.png")
      )
      
      title_main <- paste("RRHO NULL –", tiss, "\n", gsub("_vs_", " vs ", key))
      
      png(fname, width = 900, height = 800)
      plot_rrho_heatmap(rr,
                        main          = title_main,
                        palette       = palette,
                        vmax_quantile = vmax_quantile)
      dev.off()
      
      cat("Saved NULL PNG:", fname, "\n")
    }
  }
}


save_all_null_rrho_png(rrho_null,
                       outdir        = outdir_rrho,
                       palette       = rrho_rainbow_palette,
                       vmax_quantile = 0.99)









#
#
#
#
#
#
#
# reverse the orientation !!!!! bottom , right, etc. 
#
#
#
#
#
#
#
#
#




############################################################
## RRHO analysis for convergent Dry vs Mesic expression
## using DE_OG_selected_<spA>_vs_<spB>_<Tissue>_kallisto.csv
##
## IMPORTANT: In your DE code:
##   results(dds_sub, contrast = c("Species", spA, spB))
## gives log2FoldChange = log2(spA / spB) = spA - spB
##
## Here we explicitly encode which species is dry (cavity)
## and flip signs so that final RRHO score > 0 always means
## "up in dry" and < 0 means "up in humid".
############################################################

setwd("/n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/nf_pipeline/thaDol_mapping_test")

library(RRHO)
library(dplyr)
library(tidyr)
library(ggplot2)

##----------------------------------------------------------
## 1) Pairs, tissues, output directory
##----------------------------------------------------------

## spA, spB must match filenames in DE_OG_selected_<spA>_vs_<spB>_...
## dry_species tells us which of those is dry
pairs_of_interest <- tribble(
  ~pair_id,           ~spA,     ~spB,     ~dry_species,
  # DE files: DE_OG_selected_dysSti_vs_sakCri_...
  "dysSti_vs_sakCri", "dysSti", "sakCri", "sakCri",  # dry = spB
  "sakLuc_vs_sakCan", "sakLuc", "sakCan", "sakCan",  # dry = spB
  "thaAmb_vs_thaPel", "thaAmb", "thaPel", "thaPel",  # dry = spB
  # DE files: DE_OG_selected_thaBer_vs_thaAtr_...
  "thaBer_vs_thaAtr", "thaBer", "thaAtr", "thaBer",  # dry = spA
  "thaRuf_vs_thaTor", "thaRuf", "thaTor", "thaTor"   # dry = spB
)

tissues_to_run <- c("Brain", "Heart", "Liver", "Muscle", "Kidney")

outdir_rrho <- "RRHO_results"
if (!dir.exists(outdir_rrho)) dir.create(outdir_rrho, recursive = TRUE)

############################################################
## 2) Helper: build rank vector with deduplicated OGs
##    enforcing positive score = up in dry species
############################################################

## From DE:
##   log2FoldChange = spA - spB
##
## So:
##   if dry_species == spA:
##       log2FC > 0 -> up in dry (OK)
##   if dry_species == spB:
##       log2FC > 0 -> up in humid (NOT OK)
##       so we flip sign to make positive = up in dry
build_rrho_rank_vector <- function(tiss, spA, spB, dry_species,
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
  
  de_sub <- de_tab %>% dplyr::filter(!is.na(pvalue), !is.na(log2FoldChange))
  if (nrow(de_sub) == 0) {
    warning("No usable rows in ", fname)
    return(NULL)
  }
  
  # adjust log2FC so positive always = up in dry species
  log2fc <- de_sub$log2FoldChange
  if (dry_species == spB) {
    # current LFC = spA - spB; flip so becomes (spB - spA) = dry - humid
    log2fc <- -log2fc
  }
  # if dry_species == spA, leave log2fc as is (already dry - humid)
  
  # signed score: direction (dry up vs humid up) * significance
  eps   <- 1e-300
  score <- sign(log2fc) * -log10(de_sub$pvalue + eps)
  og_id <- de_sub$OG
  
  # collapse duplicates by max |score|
  score_df <- data.frame(OG = og_id, score = score, stringsAsFactors = FALSE)
  score_collapsed <- tapply(score_df$score, score_df$OG,
                            function(x) x[which.max(abs(x))])
  
  score_vec <- as.numeric(score_collapsed)
  names(score_vec) <- names(score_collapsed)
  
  # sort: largest positive first (strong up in dry)
  sort(score_vec, decreasing = TRUE)
}

############################################################
## 3) Build rrho_input for all tissues/pairs
############################################################

rrho_input <- list()

for (tiss in tissues_to_run) {
  message("\n=== Building RRHO rank vectors for tissue: ", tiss, " ===")
  rrho_input[[tiss]] <- list()
  
  for (i in seq_len(nrow(pairs_of_interest))) {
    pid  <- pairs_of_interest$pair_id[i]
    spA  <- pairs_of_interest$spA[i]
    spB  <- pairs_of_interest$spB[i]
    dry  <- pairs_of_interest$dry_species[i]
    
    v <- build_rrho_rank_vector(tiss, spA, spB, dry)
    if (is.null(v)) {
      message("  Skipping ", pid, " in ", tiss, " (no vector).")
      next
    }
    
    rrho_input[[tiss]][[pid]] <- v
    message("  Built rank vector for ", pid, " (", length(v), " unique OGs).")
  }
}

############################################################
## 4) RRHO helper to run pairwise RRHO
############################################################

run_rrho_pair <- function(vec1, vec2,
                          stepsize = 100,
                          alternative = "enrichment",
                          seed = 123) {
  set.seed(seed)
  
  if (is.null(names(vec1)) || is.null(names(vec2))) {
    stop("Both input vectors must be named (gene IDs as names).")
  }
  
  common_genes <- intersect(names(vec1), names(vec2))
  if (length(common_genes) < 10) {
    stop("Too few common genes between the two vectors.")
  }
  
  v1 <- vec1[common_genes]
  v2 <- vec2[common_genes]
  
  df1 <- data.frame(gene = common_genes,
                    score = as.numeric(v1),
                    stringsAsFactors = FALSE)
  df2 <- data.frame(gene = common_genes,
                    score = as.numeric(v2),
                    stringsAsFactors = FALSE)
  
  RRHO(df1, df2,
       stepsize    = stepsize,
       alternative = alternative,
       log10.ind   = TRUE,
       plot        = FALSE)
}

############################################################
## 5) Run empirical RRHO and save objects
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
    
    saveRDS(rr,
            file = file.path(outdir_rrho,
                             paste0("RRHO_", tiss, "_", pid1, "_vs_", pid2, "_empirical.rds")))
  }
}

message("\nEmpirical RRHO complete. Results in ", outdir_rrho, "\n")

############################################################
## 6) Build permuted (null) rank vectors with same dry/humid convention
############################################################

build_rrho_rank_vector_permuted <- function(tiss, spA, spB, dry_species,
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
  
  de_sub <- de_tab %>% dplyr::filter(!is.na(pvalue), !is.na(log2FoldChange))
  if (nrow(de_sub) == 0) {
    warning("No usable rows in ", fname)
    return(NULL)
  }
  
  # permute association across genes
  idx_perm <- sample(seq_len(nrow(de_sub)))
  de_perm  <- de_sub
  de_perm$log2FoldChange <- de_sub$log2FoldChange[idx_perm]
  de_perm$pvalue         <- de_sub$pvalue[idx_perm]
  
  # enforce dry vs humid direction as above
  log2fc <- de_perm$log2FoldChange
  if (dry_species == spB) {
    log2fc <- -log2fc
  }
  
  eps   <- 1e-300
  score <- sign(log2fc) * -log10(de_perm$pvalue + eps)
  names(score) <- de_perm$OG
  
  sort(score, decreasing = TRUE)
}

rrho_input_null <- list()

for (tiss in tissues_to_run) {
  message("\n=== Building RRHO *null* rank vectors for tissue: ", tiss, " ===")
  rrho_input_null[[tiss]] <- list()
  
  for (i in seq_len(nrow(pairs_of_interest))) {
    pid  <- pairs_of_interest$pair_id[i]
    spA  <- pairs_of_interest$spA[i]
    spB  <- pairs_of_interest$spB[i]
    dry  <- pairs_of_interest$dry_species[i]
    
    v <- build_rrho_rank_vector_permuted(tiss, spA, spB, dry)
    if (is.null(v)) {
      message("  Skipping ", pid, " in ", tiss, " (no null vector).")
      next
    }
    
    rrho_input_null[[tiss]][[pid]] <- v
    message("  Built permuted rank vector for ", pid, " (", length(v), " genes).")
  }
}

############################################################
## 7) RRHO for null
############################################################

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
    
    saveRDS(rr,
            file = file.path(outdir_rrho,
                             paste0("RRHO_", tiss, "_", pid1, "_vs_", pid2, "_null.rds")))
  }
}

message("\nRRHO (empirical + null) complete. Results in ", outdir_rrho, "\n")






##----------------------------------------------------------
## Rainbow palette + plotting helpers (screen only)
##----------------------------------------------------------

# rainbow palette: dark blue (low) → dark red (high)
rrho_rainbow_palette <- colorRampPalette(c(
  "#00004E",  # dark blue
  "#0000FF",  # blue
  "#00FFFF",  # cyan
  "#00FF00",  # green
  "#FFFF00",  # yellow
  "#FFA500",  # orange
  "#FF0000",  # red
  "#800000"   # dark red
))

# Single RRHO plot with rainbow palette and color bar
plot_rrho_heatmap <- function(rr_obj,
                              main = "RRHO",
                              ncol = 256,
                              zlim = NULL,
                              vmax_quantile = 0.99,
                              palette = rrho_rainbow_palette) {
  hm <- rr_obj$hypermat
  if (all(is.na(hm))) {
    warning("hypermat is all NA; nothing to plot.")
    return(invisible(NULL))
  }
  
  # flip vertically so origin matches usual RRHO orientation
  hm_plot <- t(hm[nrow(hm):1, , drop = FALSE])
  
  # determine zlim using a high quantile for contrast
  if (is.null(zlim)) {
    upper <- as.numeric(quantile(hm, probs = vmax_quantile, na.rm = TRUE))
    if (!is.finite(upper) || upper <= 0) {
      upper <- max(hm, na.rm = TRUE)
    }
    zlim <- c(0, upper)
  }
  
  # clip values above zlim[2] so they all get max color
  hm_plot_clipped <- pmin(hm_plot, zlim[2])
  
  # split device: 1 = heatmap, 2 = color bar
  layout(matrix(c(1, 2), ncol = 2), widths = c(4, 1))
  
  ## 1) Heatmap
  par(mar = c(5, 5, 4, 1))  # bottom, left, top, right
  image(hm_plot_clipped,
        col  = palette(ncol),
        zlim = zlim,
        axes = FALSE,
        xlab = "List 1 rank (increasing → right)",
        ylab = "List 2 rank (increasing → up)",
        main = main)
  box()
  
  ## 2) Color bar
  par(mar = c(5, 1, 4, 3))
  grad <- matrix(seq(0, 1, length.out = ncol), nrow = 1)
  image(1, seq(0, 1, length.out = ncol), grad,
        col  = palette(ncol),
        xlab = "",
        ylab = expression(-log[10](p)),
        axes = FALSE)
  at_pos  <- seq(0, 1, length.out = 5)
  at_vals <- seq(zlim[1], zlim[2], length.out = 5)
  axis(4, at = at_pos, labels = signif(at_vals, 2))
  box()
  
  # reset layout
  layout(1)
}

##----------------------------------------------------------
## Plot ALL empirical RRHO maps (screen, rainbow, no PNG)
##----------------------------------------------------------

plot_all_empirical_rrho <- function(rrho_empirical_list,
                                    palette = rrho_rainbow_palette,
                                    vmax_quantile = 0.99) {
  for (tiss in names(rrho_empirical_list)) {
    cat("\n=== Tissue:", tiss, "===\n")
    
    rr_list <- rrho_empirical_list[[tiss]]
    if (is.null(rr_list) || length(rr_list) == 0) next
    
    for (key in names(rr_list)) {
      rr <- rr_list[[key]]
      title_main <- paste("RRHO –", tiss, "\n", gsub("_vs_", " vs ", key))
      
      dev.new()
      plot_rrho_heatmap(rr,
                        main          = title_main,
                        palette       = palette,
                        vmax_quantile = vmax_quantile)
      
      cat("Plotted empirical:", tiss, "-", key, "\n")
      readline("Press [Enter] for next plot...")
      dev.off()
    }
  }
}

## Call this to view all empirical comparisons
plot_all_empirical_rrho(rrho_empirical)
