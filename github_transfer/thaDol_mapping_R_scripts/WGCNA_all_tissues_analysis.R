
setwd("/n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/nf_pipeline/thaDol_mapping_test")

library(dplyr)
library(readr)

# Read the summary of best modules
best_mods <- read_csv("WGCNA_BestModules_Dry_byTissue.csv")

# Tissues you ran
tissues_to_run <- c("Brain", "Heart", "Liver", "Muscle", "Kidney")

# List to store genes per tissue
best_module_genes <- list()

for (tiss in tissues_to_run) {
  message("Processing ", tiss, " ...")
  
  # Read tissue-specific geneInfo
  geneinfo_file <- paste0("WGCNA_geneInfo_", tiss, "_Dry_vs_Mesic.csv")
  gi <- read_csv(geneinfo_file)
  
  # Get the best module color for this tissue
  best_color <- best_mods %>%
    filter(Tissue == tiss) %>%
    pull(Color)
  
  if (length(best_color) != 1) {
    warning("Could not find a unique best color for tissue ", tiss)
    next
  }
  
  # Subset genes in the best module
  gi_best <- gi %>% filter(moduleColor == best_color)
  
  # Store full table and just gene names
  best_module_genes[[tiss]] <- gi_best
  
  # Write out a simple list of gene_name and OG for this module
  out_file <- paste0("BestModuleGenes_", tiss, "_", best_color, ".csv")
  gi_best %>% 
    select(OG, gene_name, name_source, moduleColor, GS_Dry, `p.GS_Dry`) %>%
    write_csv(out_file)
  
  message("  Tissue ", tiss, 
          ": best color = ", best_color, 
          ", n genes = ", nrow(gi_best),
          " → wrote ", out_file)
}

# Optional: Print gene names for each tissue to console
for (tiss in names(best_module_genes)) {
  cat("\n==== ", tiss, " ====\n")
  print(best_module_genes[[tiss]]$gene_name)
}




library(WGCNA)
options(stringsAsFactors = FALSE)

# Suppose you're inside run_wgcna_for_tissue() AFTER net is computed
# and you know `bestColor`

# 1. Get indices for genes in the best module
moduleGenes <- (moduleColors == bestColor)

# 2. Recalculate TOM for genes in that module only (smaller matrix)
dissTOM <- TOMdist(datExpr0[, moduleGenes], TOMType = "signed")
TOM <- 1 - dissTOM

# 3. Intramodular connectivity = rowSums(TOM) or WGCNA::intramodularConnectivity
kIN <- rowSums(TOM)  # simple measure

# 4. Quick “network-like” heatmap
png(paste0("WGCNA_", tissue_of_interest, "_", bestColor, "_TOM_heatmap.png"),
    width = 800, height = 800)
TOMplot(
  dissTOM,
  hclust(as.dist(dissTOM), method = "average"),
  as.character(which(moduleGenes)),  # labels
  main = paste("TOM heatmap –", tissue_of_interest, bestColor, "module")
)
dev.off()





# Only do this for the best module
moduleGenes <- (moduleColors == bestColor)

# Calculate TOM for the whole dataset (may be large) OR for module only
# For module only:
TOM_best <- TOMsimilarityFromExpr(datExpr0[, moduleGenes],
                                  power = softPower,
                                  TOMType = "signed",
                                  networkType = "signed hybrid")

# Choose a threshold for edges (e.g. TOM > 0.1 or 0.2)
tom_thresh <- 0.1

cyt <- exportNetworkToCytoscape(
  TOM_best,
  edgeFile      = paste0("CytoscapeEdges_", tissue_of_interest, "_", bestColor, ".txt"),
  nodeFile      = paste0("CytoscapeNodes_", tissue_of_interest, "_", bestColor, ".txt"),
  weighted      = TRUE,
  threshold     = tom_thresh,
  nodeNames     = colnames(datExpr0)[moduleGenes],
  altNodeNames  = geneInfo$gene_name[moduleGenes],
  nodeAttr      = moduleColors[moduleGenes]
)