############################################################
## 05_WGCNA_all_tissues_and_GO.R
##
## WGCNA across multiple tissues for Dry vs Mesic
## + GO enrichment on the top 3 Dry-associated modules
##
## Prerequisites:
##   - 01_DE_and_GO.R has been run, producing:
##       obj_dds_DESeq2_OG_all_tissues.rds
##       obj_og_rep_ALL_OGs_repGeneName.rds
############################################################

setwd("/n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/nf_pipeline/thaDol_mapping_test")

library(WGCNA)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

options(stringsAsFactors = FALSE)
allowWGCNAThreads()

############################################################
## 0) Load core objects
############################################################

# DESeq2 object with OGs across all tissues
dds <- readRDS("obj_dds_DESeq2_OG_all_tissues.rds")

# Global OG -> representative gene_name mapping for ALL OGs
if (file.exists("obj_og_rep_ALL_OGs_repGeneName.rds")) {
  og_rep_all <- readRDS("obj_og_rep_ALL_OGs_repGeneName.rds")
} else if (file.exists("OGs_ALL_repGeneName.csv")) {
  og_rep_all <- read.csv("OGs_ALL_repGeneName.csv", stringsAsFactors = FALSE)
} else {
  stop("Need og_rep_all (obj_og_rep_ALL_OGs_repGeneName.rds or OGs_ALL_repGeneName.csv).")
}

############################################################
## 1) Settings
############################################################

# dry vs mesic species
dry_species <- c("thaDol", "sakCri", "sakCan", "thaBer", "thaPel", "thaTor")

# Tissues to analyze
tissues_to_run <- c("Brain", "Heart", "Liver", "Muscle", "Kidney")

# WGCNA parameters
softPower      <- 6
minModuleSize  <- 30
mergeCutHeight <- 0.25

############################################################
## 2) Helper: tissue-specific background universe for GO
############################################################

get_tissue_background_entrez <- function(tiss,
                                         dds,
                                         og_rep_all,
                                         OrgDb = org.Hs.eg.db) {
  cols_tiss <- which(dds$Tissue == tiss)
  if (length(cols_tiss) == 0) {
    warning("No samples for tissue ", tiss, " in dds; cannot build background.")
    return(NULL)
  }
  
  cts_tiss <- counts(dds)[, cols_tiss, drop = FALSE]
  bg_ogs   <- rownames(cts_tiss)[rowSums(cts_tiss > 0) > 0]
  if (length(bg_ogs) < 5) {
    warning("Too few background OGs (", length(bg_ogs),
            ") in tissue ", tiss, "; skipping background.")
    return(NULL)
  }
  
  bg_map <- og_rep_all %>%
    dplyr::filter(OG %in% bg_ogs, !is.na(gene_name), gene_name != "")
  bg_symbols <- unique(bg_map$gene_name)
  
  if (length(bg_symbols) < 5) {
    warning("Too few background symbols (", length(bg_symbols),
            ") in tissue ", tiss, "; skipping background.")
    return(NULL)
  }
  
  bg_conv <- tryCatch(
    bitr(bg_symbols,
         fromType = "SYMBOL",
         toType   = "ENTREZID",
         OrgDb    = OrgDb),
    error = function(e) {
      warning("Background bitr failed for ", tiss, ": ", e$message)
      return(NULL)
    }
  )
  if (is.null(bg_conv) || nrow(bg_conv) == 0) {
    warning("No background symbols converted for tissue ", tiss)
    return(NULL)
  }
  
  unique(bg_conv$ENTREZID)
}

############################################################
## 3) Main WGCNA function for one tissue
############################################################

run_wgcna_for_tissue <- function(tissue_of_interest) {
  message("\n========================")
  message("Running WGCNA for tissue: ", tissue_of_interest)
  message("========================")
  
  ## 3.1 Expression data (normalized counts -> log2)
  norm_counts <- counts(dds, normalized = TRUE)   # OG x sample
  expr_log2   <- log2(norm_counts + 1)
  
  coldata_dds <- as.data.frame(colData(dds))
  coldata_dds$Sample <- rownames(coldata_dds)
  stopifnot(all(colnames(norm_counts) == coldata_dds$Sample))
  
  samples_tiss <- rownames(coldata_dds)[coldata_dds$Tissue == tissue_of_interest]
  if (length(samples_tiss) < 5) {
    warning("Tissue ", tissue_of_interest, " has fewer than 5 samples. Skipping.")
    return(NULL)
  }
  
  expr_tiss <- expr_log2[, samples_tiss, drop = FALSE]  # OG x samples
  cd_tiss   <- coldata_dds[samples_tiss, , drop = FALSE]
  
  # WGCNA expects samples as rows, genes as columns
  datExpr0 <- t(expr_tiss)  # samples x OG
  
  ## 3.2 Filter lowly expressed OGs
  norm_tiss <- norm_counts[, samples_tiss, drop = FALSE]  # OG x samples
  keep_ogs  <- rowSums(norm_tiss >= 15) >= (0.10 * length(samples_tiss))
  datExpr0  <- datExpr0[, keep_ogs, drop = FALSE]
  
  cat("After filtering, ", ncol(datExpr0), " OGs remain in tissue ", tissue_of_interest, "\n")
  
  ## 3.3 Good samples / genes
  gsg <- goodSamplesGenes(datExpr0, verbose = 3)
  if (!gsg$allOK) {
    if (sum(!gsg$goodGenes) > 0)
      printFlush(paste("Removing genes:",
                       paste(colnames(datExpr0)[!gsg$goodGenes], collapse = ", ")))
    if (sum(!gsg$goodSamples) > 0)
      printFlush(paste("Removing samples:",
                       paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")))
    datExpr0 <- datExpr0[gsg$goodSamples, gsg$goodGenes, drop = FALSE]
    cd_tiss  <- cd_tiss[rownames(datExpr0), , drop = FALSE]
  }
  
  ## 3.4 Trait: Dry vs Mesic
  dry_trait <- ifelse(cd_tiss$Species %in% dry_species, 1, 0)
  datTraits <- data.frame(
    Dry = dry_trait,
    row.names = rownames(cd_tiss)
  )
  
  nGenes   <- ncol(datExpr0)
  nSamples <- nrow(datExpr0)
  
  ## 3.5 Network construction
  net <- blockwiseModules(
    datExpr0,
    power               = softPower,
    TOMType             = "signed",
    networkType         = "signed hybrid",
    minModuleSize       = minModuleSize,
    reassignThreshold   = 0,
    mergeCutHeight      = mergeCutHeight,
    numericLabels       = TRUE,
    pamRespectsDendro   = FALSE,
    saveTOMs            = FALSE,
    corType             = "bicor",
    maxBlockSize        = nGenes,
    verbose             = 3
  )
  
  moduleLabels <- net$colors              # numeric labels (0,1,2,...)
  moduleColors <- labels2colors(moduleLabels)  # color names ("turquoise", "blue", ...)
  MEs          <- orderMEs(net$MEs)           # module eigengenes
  geneTree     <- net$dendrograms[[1]]
  
  ## 3.6 Dendrogram + module colors
  png(paste0("WGCNA_", tissue_of_interest, "_dendrogram_modules.png"),
      width = 800, height = 600)
  par(mfrow = c(1,1), mar = c(0,4,2,0))
  plotDendroAndColors(
    geneTree,
    moduleColors[net$blockGenes[[1]]],
    "Module colors",
    dendroLabels = FALSE,
    hang         = 0.03,
    addGuide     = TRUE,
    guideHang    = 0.05,
    main         = paste("Gene dendrogram and module colors –", tissue_of_interest)
  )
  dev.off()
  
  ## 3.7 Module–trait correlations
  modTraitCor <- cor(MEs, datTraits, use = "p")
  modTraitP   <- corPvalueStudent(modTraitCor, nSamples)
  
  textMatrix <- paste(signif(modTraitCor, 2), "\n(",
                      signif(modTraitP,   1), ")", sep = "")
  dim(textMatrix) <- dim(modTraitCor)
  
  png(paste0("WGCNA_", tissue_of_interest, "_ModuleTrait_Dry.png"),
      width = 800, height = 800)
  par(mar = c(6, 8.5, 3, 1))
  labeledHeatmap(
    Matrix      = modTraitCor,
    xLabels     = colnames(datTraits),
    yLabels     = rownames(modTraitCor),
    ySymbols    = rownames(modTraitCor),
    colorLabels = FALSE,
    colors      = blueWhiteRed(50),
    textMatrix  = textMatrix,
    cex.text    = 0.7,
    zlim        = c(-1, 1),
    main        = paste("Module–Dry relationships –", tissue_of_interest)
  )
  dev.off()
  
  ## 3.8 Gene-level measures: MM and GS(Dry)
  modNames <- substring(colnames(MEs), 3)  # "0","1",...
  
  geneModuleMembership <- as.data.frame(cor(datExpr0, MEs, use = "p"))
  MMPvalue             <- as.data.frame(
    corPvalueStudent(as.matrix(geneModuleMembership), nSamples)
  )
  names(geneModuleMembership) <- paste0("MM_", modNames)
  names(MMPvalue)             <- paste0("p.MM_", modNames)
  
  geneTraitSignificance <- as.data.frame(cor(datExpr0, datTraits$Dry, use = "p"))
  GSPvalue              <- as.data.frame(
    corPvalueStudent(as.matrix(geneTraitSignificance), nSamples)
  )
  names(geneTraitSignificance) <- "GS_Dry"
  names(GSPvalue)              <- "p.GS_Dry"
  
  geneInfo <- data.frame(
    OG              = colnames(datExpr0),
    moduleColor     = moduleColors,
    geneTraitSignificance,
    GSPvalue,
    geneModuleMembership,
    MMPvalue,
    stringsAsFactors = FALSE
  )
  
  # Add common gene names if available
  if (exists("og_rep_all")) {
    geneInfo <- geneInfo %>%
      dplyr::left_join(og_rep_all, by = "OG")
  }
  
  out_geneInfo <- paste0("WGCNA_geneInfo_", tissue_of_interest, "_Dry_vs_Mesic.csv")
  write.csv(geneInfo, out_geneInfo, row.names = FALSE)
  message("Wrote: ", out_geneInfo)
  
  ## 3.9 Rank modules by |cor(ME, Dry)| and pick top 3
  absCor <- abs(modTraitCor[, "Dry"])
  top3_idx <- order(absCor, decreasing = TRUE)[1:min(3, length(absCor))]
  top3_MEs <- rownames(modTraitCor)[top3_idx]   # e.g. "ME1", "ME2", ...
  
  # Map ME name to actual module color name (string)
  get_color_for_ME <- function(me) {
    lab <- as.numeric(sub("^ME", "", me))          # numeric label
    cols <- unique(moduleColors[moduleLabels == lab])
    if (length(cols) == 0) NA_character_ else cols[1]
  }
  
  top3_colors <- sapply(top3_MEs, get_color_for_ME)
  
  top3_summary <- data.frame(
    Tissue  = tissue_of_interest,
    ME      = top3_MEs,
    Color   = top3_colors,                        # e.g. "turquoise", "blue", ...
    r_Dry   = modTraitCor[top3_MEs, "Dry"],
    p_Dry   = modTraitP[top3_MEs, "Dry"],
    stringsAsFactors = FALSE
  )
  
  # 3.10 Plots for the best (rank 1) module
  bestME    <- top3_MEs[1]
  bestColor <- top3_colors[1]
  
  cat("Tissue", tissue_of_interest, ": best Dry module =",
      bestME, "(color", bestColor, "), r =",
      round(modTraitCor[bestME, "Dry"], 3), ", p =",
      signif(modTraitP[bestME, "Dry"], 3), "\n")
  
  ME_best <- MEs[, bestME]
  df_me <- data.frame(
    Sample  = rownames(datExpr0),
    ME      = ME_best,
    Species = cd_tiss$Species,
    Dry     = datTraits$Dry
  )
  
  png(paste0("WGCNA_", tissue_of_interest, "_Eigengene_", bestColor, "_byDry.png"),
      width = 600, height = 500)
  boxplot(ME ~ Dry, data = df_me,
          names = c("Mesic", "Dry"),
          main  = paste(tissue_of_interest, bestME, "(", bestColor, ")", "by Dry vs Mesic"),
          ylab  = "Module eigengene")
  dev.off()
  
  # MM vs GS for best module
  bestLabel <- as.numeric(sub("^ME", "", bestME))
  moduleGenes <- (moduleLabels == bestLabel)
  mm_colname  <- paste0("MM_", bestLabel)
  if (mm_colname %in% names(geneModuleMembership) && any(moduleGenes)) {
    png(paste0("WGCNA_", tissue_of_interest, "_MM_vs_GS_", bestColor, ".png"),
        width = 600, height = 500)
    verboseScatterplot(
      x = abs(geneModuleMembership[moduleGenes, mm_colname]),
      y = abs(geneTraitSignificance[moduleGenes, 1]),
      xlab = paste("Module Membership in", bestColor, "module"),
      ylab = "Gene significance for Dry",
      main = paste("MM vs GS_Dry –", tissue_of_interest, "(", bestColor, ")"),
      cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.1,
      col = bestColor
    )
    dev.off()
  }
  
  list(
    tissue       = tissue_of_interest,
    geneInfo     = geneInfo,
    MEs          = MEs,
    modTraitCor  = modTraitCor,
    modTraitP    = modTraitP,
    top3_summary = top3_summary
  )
}

############################################################
## 4) Run WGCNA for all tissues
############################################################

wgcna_results <- lapply(tissues_to_run, run_wgcna_for_tissue)
names(wgcna_results) <- tissues_to_run

# (optional) save results object
saveRDS(wgcna_results, file = "obj_WGCNA_all_tissues_results.rds")

# Summary of top 3 modules per tissue
top3_all <- do.call(rbind, lapply(wgcna_results, function(x) {
  if (is.null(x)) return(NULL)
  x$top3_summary
}))
write.csv(top3_all, "WGCNA_Top3Modules_Dry_byTissue.csv", row.names = FALSE)
print(top3_all)

############################################################
## 5) GO enrichment for top 3 modules per tissue
############################################################

for (tiss in tissues_to_run) {
  res_t <- wgcna_results[[tiss]]
  if (is.null(res_t)) next
  
  geneInfo   <- res_t$geneInfo
  top3       <- res_t$top3_summary
  bg_universe <- get_tissue_background_entrez(tiss, dds, og_rep_all, org.Hs.eg.db)
  if (is.null(bg_universe)) {
    warning("No background universe for ", tiss, "; skipping GO for its modules.")
    next
  }
  
  for (i in seq_len(nrow(top3))) {
    me_name  <- top3$ME[i]        # e.g. "ME1"
    color    <- top3$Color[i]     # e.g. "turquoise"
    r_val    <- top3$r_Dry[i]
    p_val    <- top3$p_Dry[i]
    
    message("Tissue ", tiss, " – module ", me_name,
            " (color ", color, "), r_Dry = ", round(r_val,3),
            ", p_Dry = ", signif(p_val,3))
    
    genes_mod <- geneInfo %>%
      dplyr::filter(moduleColor == color,
                    !is.na(gene_name), gene_name != "")
    
    gene_symbols <- unique(genes_mod$gene_name)
    if (length(gene_symbols) < 5) {
      warning("Too few genes (", length(gene_symbols),
              ") in tissue ", tiss, " module ", color, " for GO; skipping.")
      next
    }
    
    gmap <- tryCatch(
      bitr(gene_symbols,
           fromType = "SYMBOL",
           toType   = "ENTREZID",
           OrgDb    = org.Hs.eg.db),
      error = function(e) {
        warning("bitr failed for ", tiss, " module ", color, ": ", e$message)
        return(NULL)
      }
    )
    if (is.null(gmap) || nrow(gmap) == 0) {
      warning("No symbols converted to Entrez for ", tiss, " module ", color)
      next
    }
    
    ego <- tryCatch(
      enrichGO(
        gene          = gmap$ENTREZID,
        universe      = bg_universe,
        OrgDb         = org.Hs.eg.db,
        keyType       = "ENTREZID",
        ont           = "BP",
        pAdjustMethod = "BH",
        pvalueCutoff  = 0.05,
        qvalueCutoff  = 0.05,
        readable      = TRUE
      ),
      error = function(e) {
        warning("enrichGO failed for ", tiss, " module ", color, ": ", e$message)
        return(NULL)
      }
    )
    
    if (is.null(ego) || nrow(as.data.frame(ego)) == 0) {
      message("No enriched GO terms for ", tiss, " module ", color)
      next
    }
    
    ego_df <- as.data.frame(ego)
    tag <- paste0(tiss, "_", color)
    
    out_tab <- paste0("GO_BP_Hs_WGCNA_", tag, "_enrichment.csv")
    write.csv(ego_df, out_tab, row.names = FALSE)
    message("  Wrote GO table: ", out_tab)
    
    p_dot <- dotplot(ego, showCategory = 20) +
      ggtitle(paste0("GO BP – ", tiss, " WGCNA module ", color,
                     "\n(r_Dry = ", round(r_val,3),
                     ", p = ", signif(p_val,3), ")")) +
      theme_bw() +
      theme(
        plot.title   = element_text(size = 11),
        axis.text.x  = element_text(size = 8),
        axis.text.y  = element_text(size = 8),
        axis.title.x = element_text(size = 9),
        axis.title.y = element_text(size = 9)
      )
    
    out_dot <- paste0("GO_BP_Hs_WGCNA_", tag, "_dotplot.png")
    ggsave(out_dot, p_dot, width = 8, height = 6, dpi = 300)
    message("  Wrote GO dotplot: ", out_dot)
  }
}

############################################################
## End of 05_WGCNA_all_tissues_and_GO.R
############################################################

