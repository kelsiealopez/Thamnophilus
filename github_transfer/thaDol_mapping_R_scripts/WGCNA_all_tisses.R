############################################################
## WGCNA across multiple tissues for Dry vs Mesic
## Tissues: Brain, Heart, Liver, Muscle, Kidney
############################################################

# need to run this code first !DEseq_fastOMA_OGs_KALLISTO_cleanedUp_testMaptoDoliatus_only_fastOMA_OGs.R

library(WGCNA)
options(stringsAsFactors = FALSE)
allowWGCNAThreads()

# dry vs mesic species
dry_species <- c("thaDol", "sakCri", "sakCan", "thaBer", "thaPel", "thaTor")

# Tissues to analyze
tissues_to_run <- c("Brain", "Heart", "Liver", "Muscle", "Kidney")

# Common WGCNA parameters, as in the paper
softPower      <- 6
minModuleSize  <- 30
mergeCutHeight <- 0.25

############################################################
## Helper: run WGCNA for one tissue
############################################################

run_wgcna_for_tissue <- function(tissue_of_interest) {
  message("\n========================")
  message("Running WGCNA for tissue: ", tissue_of_interest)
  message("========================")
  
  # 1. Expression matrices (normalized counts -> log2)
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
  
  # 2. Filter lowly expressed OGs: keep OG if >=15 counts in at least 10% of samples
  norm_tiss <- norm_counts[, samples_tiss, drop = FALSE]  # OG x samples
  keep_ogs  <- rowSums(norm_tiss >= 15) >= (0.10 * length(samples_tiss))
  datExpr0  <- datExpr0[, keep_ogs, drop = FALSE]
  
  cat("After filtering, ", ncol(datExpr0), " OGs remain in tissue ", tissue_of_interest, "\n")
  
  # 3. Check good samples & genes
  gsg <- goodSamplesGenes(datExpr0, verbose = 3)
  if (!gsg$allOK) {
    if (sum(!gsg$goodGenes) > 0)
      printFlush(paste("Removing genes:",
                       paste(colnames(datExpr0)[!gsg$goodGenes], collapse = ", ")))
    if (sum(!gsg$goodSamples) > 0)
      printFlush(paste("Removing samples:",
                       paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")))
    datExpr0 <- datExpr0[gsg$goodSamples, gsg$goodGenes]
    cd_tiss  <- cd_tiss[rownames(datExpr0), , drop = FALSE]
  }
  
  # 4. Trait: Dry vs Mesic (1/0)
  dry_trait <- ifelse(cd_tiss$Species %in% dry_species, 1, 0)
  datTraits <- data.frame(
    Dry = dry_trait,
    row.names = rownames(cd_tiss)
  )
  
  # 5. Network construction
  nGenes   <- ncol(datExpr0)
  nSamples <- nrow(datExpr0)
  
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
  
  moduleLabels <- net$colors
  moduleColors <- labels2colors(moduleLabels)
  MEs          <- orderMEs(net$MEs)
  geneTree     <- net$dendrograms[[1]]
  
  # Save dendrogram + module colors
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
  
  # 6. Module–trait correlations
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
    yLabels     = names(MEs),
    ySymbols    = names(MEs),
    colorLabels = FALSE,
    colors      = blueWhiteRed(50),
    textMatrix  = textMatrix,
    cex.text    = 0.7,
    zlim        = c(-1, 1),
    main        = paste("Module–Dry relationships –", tissue_of_interest)
  )
  dev.off()
  
  # 7. Gene-level measures: MM and GS
  modNames <- substring(names(MEs), 3)  # "0","1",...
  
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
  
  if (exists("og_rep_all")) {
    geneInfo <- geneInfo %>%
      dplyr::left_join(og_rep_all, by = c("OG" = "OG"))
  }
  
  out_geneInfo <- paste0("WGCNA_geneInfo_", tissue_of_interest, "_Dry_vs_Mesic.csv")
  write.csv(geneInfo, out_geneInfo, row.names = FALSE)
  message("Wrote: ", out_geneInfo)
  
  ##########################################################
  ## 8. Pick the most Dry-associated module & make plots
  ##########################################################
  
  absCor <- abs(modTraitCor[,"Dry"])
  bestMEidx <- which.max(absCor)
  bestME    <- names(MEs)[bestMEidx]   # e.g. "ME11"
  bestLabel <- as.numeric(gsub("^ME", "", bestME))
  bestColor <- unique(moduleColors[moduleLabels == bestLabel])
  
  cat("Tissue", tissue_of_interest, ": best Dry module =",
      bestME, "(label", bestLabel, ", color", bestColor, "), r =",
      round(modTraitCor[bestME, "Dry"], 3), ", p =",
      signif(modTraitP[bestME, "Dry"], 3), "\n")
  
  # 8a. Eigengene vs Dry boxplot
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
  
  # 8b. Module membership vs GS(Dry) for the best module
  moduleGenes <- moduleColors == bestColor
  mm_colname  <- paste0("MM_", bestLabel)
  mm_col_idx  <- which(names(geneModuleMembership) == mm_colname)
  
  if (length(mm_col_idx) == 1 && any(moduleGenes)) {
    png(paste0("WGCNA_", tissue_of_interest, "_MM_vs_GS_", bestColor, ".png"),
        width = 600, height = 500)
    verboseScatterplot(
      x = abs(geneModuleMembership[moduleGenes, mm_col_idx]),
      y = abs(geneTraitSignificance[moduleGenes, 1]),
      xlab = paste("Module Membership in", bestColor, "module"),
      ylab = "Gene significance for Dry",
      main = paste("MM vs GS_Dry –", tissue_of_interest, "(", bestColor, ")"),
      cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.1,
      col = bestColor
    )
    dev.off()
  } else {
    warning("Could not make MM vs GS plot for ", tissue_of_interest,
            " best module; mm_col_idx length = ", length(mm_col_idx),
            ", any(moduleGenes) = ", any(moduleGenes))
  }
  
  # Return a small summary list
  list(
    tissue      = tissue_of_interest,
    bestME      = bestME,
    bestLabel   = bestLabel,
    bestColor   = bestColor,
    r_Dry       = modTraitCor[bestME, "Dry"],
    p_Dry       = modTraitP[bestME, "Dry"],
    geneInfo    = geneInfo
  )
}

############################################################
## Run WGCNA for all tissues
############################################################

wgcna_results <- lapply(tissues_to_run, run_wgcna_for_tissue)
names(wgcna_results) <- tissues_to_run

# Inspect summary of best modules per tissue
best_modules_summary <- do.call(rbind, lapply(wgcna_results, function(x) {
  if (is.null(x)) return(NULL)
  data.frame(
    Tissue    = x$tissue,
    ME        = x$bestME,
    Label     = x$bestLabel,
    Color     = x$bestColor,
    r_Dry     = x$r_Dry,
    p_Dry     = x$p_Dry,
    stringsAsFactors = FALSE
  )
}))

print(best_modules_summary)
write.csv(best_modules_summary,
          "WGCNA_BestModules_Dry_byTissue.csv",
          row.names = FALSE)