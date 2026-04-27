## Run the following code in its entirety before running this because you need the objects

# DEseq_fastOMA_OGs_KALLISTO_cleanedUp_testMaptoDoliatus_only_fastOMA_OGs.R



############################################################
## WGCNA: dry vs mesic species on Thamnophilus OG counts
## Assumes you already have:
##   - dds (DESeqDataSet) with OG rows, samples as columns
##   - coldata or colData(dds) with columns Species and Tissue
##   - og_rep_all (optional) mapping OG -> gene_name
############################################################
#install.packages("WGCNA")
library(WGCNA)



## 1. Install BiocManager if needed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

## 2. Install the missing Bioconductor dependency
BiocManager::install("impute")

## 3. Now install WGCNA (this will reuse the already installed fastcluster, dynamicTreeCut)
install.packages("WGCNA")

## 4. Load it
library(WGCNA)
options(stringsAsFactors = FALSE)
allowWGCNAThreads()

############################################################
## 0. Input objects and basic setup
############################################################

# If you don't already have these in memory, load them:
# load("EvoGeneX_Thamnophilus_allTissues.RData")  # for dds, coldata, etc.
# og_rep_all <- read.csv("OGs_ALL_repGeneName.csv", stringsAsFactors = FALSE)

# Define dry species (must match Species names in coldata/dds)
dry_species <- c("thaDol", "sakCri", "sakCan", "thaBer", "thaPel", "thaTor")

# Pick ONE tissue for now (you can loop over tissues later)
tissue_of_interest <- "Kidney"   # change to "Brain","Heart","Liver","Muscle" etc.

############################################################
## 1. Extract normalized expression matrix for this tissue
############################################################

# Get DESeq2 normalized counts and log2-transform
norm_counts <- counts(dds, normalized = TRUE)
expr_log2   <- log2(norm_counts + 1)

# Build a coldata data.frame from dds
coldata_dds <- as.data.frame(colData(dds))
coldata_dds$Sample <- rownames(coldata_dds)

stopifnot(all(colnames(norm_counts) == coldata_dds$Sample))

# Subset samples to this tissue
samples_tiss <- rownames(coldata_dds)[coldata_dds$Tissue == tissue_of_interest]
expr_tiss    <- expr_log2[, samples_tiss, drop = FALSE]   # OG x samples
cd_tiss      <- coldata_dds[coldata_dds$Tissue == tissue_of_interest, , drop = FALSE]

# WGCNA expects samples as rows, genes (OGs) as columns
datExpr0 <- t(expr_tiss)  # samples x OG

# Optional: rename columns OG -> something friendlier; here we keep OG IDs
# colnames(datExpr0) <- rownames(expr_tiss)

############################################################
## 2. Basic filtering: remove lowly expressed OGs
##    (like in the paper: <15 norm counts in 90% of samples)
############################################################

# Filter done on normalized counts (not log2) for this tissue
norm_tiss <- norm_counts[, samples_tiss, drop = FALSE]   # OG x samples

keep_ogs <- rowSums(norm_tiss >= 15) >= (0.10 * length(samples_tiss))
table(keep_ogs)

datExpr0 <- datExpr0[, keep_ogs, drop = FALSE]   # samples x OG (filtered)

cat("After filtering, ", ncol(datExpr0), " OGs remain in tissue ", tissue_of_interest, "\n")

############################################################
## 3. Check for good samples/genes and remove outliers
############################################################

gsg <- goodSamplesGenes(datExpr0, verbose = 3)
gsg$allOK

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

# Sample clustering to detect outliers
sampleTree <- hclust(dist(datExpr0), method = "average")

par(cex = 0.6, mar = c(0,4,2,0))
plot(sampleTree, main = paste("Sample clustering –", tissue_of_interest),
     xlab = "", sub = "")

## OPTIONAL: if you see clear outliers, set a cut height and remove them.
## Here we assume no outlier removal, but you can uncomment/edit:

# cutHeight <- 100    # CHANGE if you see a horizontal break in the dendrogram
# abline(h = cutHeight, col = "red")
# clust <- cutreeStatic(sampleTree, cutHeight = cutHeight, minSize = 5)
# keepSamples <- (clust == 1)
# datExpr0    <- datExpr0[keepSamples, ]
# cd_tiss     <- cd_tiss[rownames(datExpr0), , drop = FALSE]
# sampleTree  <- hclust(dist(datExpr0), method = "average")

############################################################
## 4. Construct trait data: Dry vs Mesic
############################################################

# For each sample: 1 = dry, 0 = mesic
dry_trait <- ifelse(cd_tiss$Species %in% dry_species, 1, 0)
datTraits <- data.frame(
  Dry = dry_trait,
  row.names = rownames(cd_tiss)
)

############################################################
## 5. Pick soft-thresholding power (or set to 6 like paper)
############################################################

# You can run this to inspect:
powers <- c(1:10, 12, 14, 16, 18, 20)
sft <- pickSoftThreshold(datExpr0,
                         powerVector = powers,
                         verbose = 5,
                         networkType = "signed hybrid",
                         corFnc = "bicor",
                         corOptions = list(use = "pairwise.complete.obs"))

# Plot scale-free topology fit index vs power
par(mfrow = c(1,2), mar = c(4,4,2,1))
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3]) * sft$fitIndices[,2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     type = "n", main = paste("Scale independence -", tissue_of_interest))
text(sft$fitIndices[,1],
     -sign(sft$fitIndices[,3]) * sft$fitIndices[,2],
     labels = sft$fitIndices[,1], col = "red")
abline(h = 0.80, col = "red")

# Plot mean connectivity vs power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity", type = "n",
     main = paste("Mean connectivity -", tissue_of_interest))
text(sft$fitIndices[,1], sft$fitIndices[,5],
     labels = sft$fitIndices[,1], col = "red")

# Choose power (either inspect sft$powerEstimate or force 6 as in the paper)
softPower <- 6
cat("Using softPower =", softPower, "\n")

############################################################
## 6. Network construction + module detection (blockwiseModules)
############################################################

net <- blockwiseModules(
  datExpr0,
  power               = softPower,
  TOMType             = "signed",           # with networkType="signed hybrid" below
  networkType         = "signed hybrid",    # as in paper
  minModuleSize       = 30,
  reassignThreshold   = 0,
  mergeCutHeight      = 0.25,               # merge modules with eigengene cor > 0.75
  numericLabels       = TRUE,
  pamRespectsDendro   = FALSE,
  saveTOMs            = FALSE,
  corType             = "bicor",
  maxBlockSize        = ncol(datExpr0),
  verbose             = 3
)

moduleLabels  <- net$colors
moduleColors  <- labels2colors(moduleLabels)
MEs           <- net$MEs
MEs           <- orderMEs(MEs)
geneTree      <- net$dendrograms[[1]]

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

############################################################
## 7. Relate modules to the Dry trait
############################################################

# Correlate eigengenes with trait
nSamples <- nrow(datExpr0)
modTraitCor <- cor(MEs, datTraits, use = "p")
modTraitP   <- corPvalueStudent(modTraitCor, nSamples)

# Heatmap with correlations and p-values
textMatrix <- paste(signif(modTraitCor, 2), "\n(",
                    signif(modTraitP,   1), ")", sep = "")
dim(textMatrix) <- dim(modTraitCor)

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

############################################################
## 8. Gene-level measures: Module Membership & Gene Significance
############################################################

# Remove "ME" prefix for convenience
modNames <- substring(names(MEs), 3)

# Module membership: correlation of each gene with each module eigengene
geneModuleMembership <- as.data.frame(cor(datExpr0, MEs, use = "p"))
MMPvalue            <- as.data.frame(
  corPvalueStudent(as.matrix(geneModuleMembership), nSamples)
)
names(geneModuleMembership) <- paste0("MM_", modNames)
names(MMPvalue)            <- paste0("p.MM_", modNames)

# Gene significance for Dry (trait)
geneTraitSignificance <- as.data.frame(cor(datExpr0, datTraits$Dry, use = "p"))
GSPvalue             <- as.data.frame(
  corPvalueStudent(as.matrix(geneTraitSignificance), nSamples)
)
names(geneTraitSignificance) <- "GS_Dry"
names(GSPvalue)              <- "p.GS_Dry"

# Combine into one table
geneInfo <- data.frame(
  OG              = colnames(datExpr0),
  moduleColor     = moduleColors,
  geneTraitSignificance,
  GSPvalue,
  geneModuleMembership,
  MMPvalue,
  stringsAsFactors = FALSE
)

# Optionally map OG -> gene_name if og_rep_all available
if (exists("og_rep_all")) {
  geneInfo <- geneInfo %>%
    dplyr::left_join(og_rep_all, by = c("OG" = "OG"))
}

write.csv(
  geneInfo,
  file = paste0("WGCNA_geneInfo_", tissue_of_interest, "_Dry_vs_Mesic.csv"),
  row.names = FALSE
)

############################################################
## 9. Example: scatterplot for a Dry-associated module
############################################################

# 1. Find index of ME with strongest correlation to Dry
absCor <- abs(modTraitCor[,"Dry"])
bestMEidx <- which.max(absCor)      # index in MEs / modNames
bestME    <- names(MEs)[bestMEidx]  # e.g. "MEturquoise" or "ME11"

# 2. Get the numeric module label underlying that ME
bestLabel <- as.numeric(gsub("^ME", "", bestME))  # e.g. 11

# 3. Convert numeric label -> actual color name
bestColor <- unique(moduleColors[moduleLabels == bestLabel])

cat("Most Dry-associated module in", tissue_of_interest,
    "is label", bestLabel, "color", bestColor, "\n")

# 4. Now select genes in that color
moduleGenes <- moduleColors == bestColor

moduleColumn <- bestMEidx  # column in MEs / geneModuleMembership

par(mfrow = c(1,1), mar = c(5,5,3,1))
verboseScatterplot(
  x = abs(geneModuleMembership[moduleGenes, moduleColumn]),
  y = abs(geneTraitSignificance[moduleGenes, 1]),
  xlab = paste("Module Membership in", bestColor, "module"),
  ylab = "Gene significance for Dry",
  main = paste("Module membership vs gene significance\n", tissue_of_interest),
  cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.1, col = bestColor
)





modTraitCor["MEgreenyellow", "Dry"]
modTraitP["MEgreenyellow", "Dry"]


# Genes in greenyellow module
greenGenes <- geneInfo$OG[geneInfo$moduleColor == "greenyellow"]

# Tissue-specific expressed background (you already implemented this)
bg_ogs <- rownames(norm_tiss)[rowSums(norm_tiss > 0) > 0]
# map to gene_name via og_rep_all, then do clusterProfiler::enrichGO


MEg <- MEs[ , "MEgreenyellow"]
df_me <- data.frame(
  Sample  = rownames(datExpr0),
  ME      = MEg,
  Species = cd_tiss$Species,
  Dry     = datTraits$Dry
)

boxplot(ME ~ Dry, data = df_me,
        names = c("Mesic", "Dry"),
        main  = "Kidney MEgreenyellow by Dry vs Mesic")








############################################################
## 1. Inspect ME names and module–trait correlations
############################################################

# See what MEs you actually have
print(names(MEs))
print(rownames(modTraitCor))

# For safety, print the whole correlation rownames
# rownames(modTraitCor)  # should be like "ME0","ME1",...

############################################################
## 2. Map color "greenyellow" -> numeric label -> ME column
############################################################

# moduleLabels: numeric labels (0,1,2,...) from blockwiseModules
# moduleColors: same length, color names
# Find numeric label(s) for greeyellow
gy_labels <- unique(moduleLabels[moduleColors == "greenyellow"])
gy_labels
# Often this will be a single number, e.g. 14

if (length(gy_labels) != 1) {
  stop("Expected exactly one numeric label for greeyellow, found: ",
       paste(gy_labels, collapse = ", "))
}

gy_label <- gy_labels[1]
gy_ME_name <- paste0("ME", gy_label)
gy_ME_name

# Confirm this ME exists
if (!gy_ME_name %in% names(MEs)) {
  stop("Eigengene ", gy_ME_name, " not found in names(MEs). Check labels.")
}

# 2a. Correlation of greeyellow eigengene with Dry
gy_cor <- modTraitCor[gy_ME_name, "Dry"]
gy_p   <- modTraitP[gy_ME_name, "Dry"]
cat("Kidney:", gy_ME_name, "(color = greeyellow) vs Dry: r =",
    round(gy_cor, 3), ", p =", signif(gy_p, 3), "\n")

############################################################
## 3. Eigengene vs Dry boxplot
############################################################

MEg <- MEs[, gy_ME_name]

df_me <- data.frame(
  Sample  = rownames(datExpr0),
  ME      = MEg,
  Species = cd_tiss$Species,
  Dry     = datTraits$Dry
)

boxplot(ME ~ Dry, data = df_me,
        names = c("Mesic", "Dry"),
        main  = paste("Kidney", gy_ME_name, "(greenyellow) by Dry vs Mesic"),
        ylab  = "Module eigengene")

############################################################
## 4. Scatter: module membership vs GS for greeyellow
############################################################

# Column in geneModuleMembership corresponding to this ME
gy_MM_col <- which(names(geneModuleMembership) == paste0("MM_", gy_label))
if (length(gy_MM_col) != 1) {
  stop("Could not uniquely find MM column for label ", gy_label)
}

gy_module_genes <- moduleColors == "greenyellow"

verboseScatterplot(
  x = abs(geneModuleMembership[gy_module_genes, gy_MM_col]),
  y = abs(geneTraitSignificance[gy_module_genes, 1]),
  xlab = "Module Membership in greeyellow module",
  ylab = "Gene significance for Dry",
  main = paste("Kidney greeyellow (", gy_ME_name, ") – MM vs GS_Dry", sep = ""),
  cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.1,
  col = "greenyellow"
)

############################################################
## 5. Export greeyellow module network to Cytoscape
############################################################

# Subset TOM to greeyellow genes
inModule <- which(moduleColors == "greenyellow")
modTOM  <- TOM[inModule, inModule, drop = FALSE]

# Make sure nodeNames length matches nrow(modTOM)
nodeNames_sub <- colnames(datExpr0)[inModule]

if (nrow(modTOM) != length(nodeNames_sub)) {
  stop("Mismatch: nrow(modTOM) =", nrow(modTOM),
       " but length(nodeNames_sub) =", length(nodeNames_sub))
}

# Optional alternative names from og_rep_all
if (exists("og_rep_all")) {
  nm <- nodeNames_sub
  m2 <- og_rep_all[match(nm, og_rep_all$OG), "gene_name"]
  altNames <- ifelse(is.na(m2) | m2 == "", nm, m2)
} else {
  altNames <- nodeNames_sub
}

cyt <- exportNetworkToCytoscape(
  modTOM,
  edgeFile  = paste0("CytoscapeInput-edges_", tissue_of_interest, "_greenyellow.txt"),
  nodeFile  = paste0("CytoscapeInput-nodes_", tissue_of_interest, "_greenyellow.txt"),
  threshold = 0.1,
  weighted  = TRUE,
  nodeNames = nodeNames_sub,
  altNodeNames = altNames,
  nodeAttr = moduleColors[inModule]
)

############################################################
## 6. Double-check gene list for greeyellow (optional)
############################################################

gy_genes <- geneInfo[geneInfo$moduleColor == "greenyellow",
                     c("OG", "gene_name", "GS_Dry", "moduleColor")]
head(gy_genes)
