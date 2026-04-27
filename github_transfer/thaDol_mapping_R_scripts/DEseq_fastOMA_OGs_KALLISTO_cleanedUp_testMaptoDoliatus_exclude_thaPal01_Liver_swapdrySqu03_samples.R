############################################################
## Thamnophilus DE + PCA + Convergence + GO (human only)
##   – using *.kallisto.merged.gene_counts.tsv
############################################################

setwd("/n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/nf_pipeline/thaDol_mapping_test")

library(DESeq2)
library(ggplot2)
library(RColorBrewer)
library(ggrepel)
library(matrixStats)
library(jsonlite)
library(dplyr)
library(tibble)
library(tidyr)
library(ComplexUpset)
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

## ========================================================
## 1) Read gene-level kallisto counts and build DESeq2 object
## ========================================================

## 1a) List your new files
files <- c(
  "drySqu_to_thaDol_kallisto.merged.gene_counts.tsv",
  "dysSti_to_thaDol_kallisto.merged.gene_counts.tsv",
  "sakCan_to_thaDol_kallisto.merged.gene_counts.tsv",
  "sakCri_to_thaDol_kallisto.merged.gene_counts.tsv",
  "sakLuc_to_thaDol_kallisto.merged.gene_counts.tsv",
  "thaAmb_to_thaDol_kallisto.merged.gene_counts.tsv",
  "thaAtr_to_thaDol_kallisto.merged.gene_counts.tsv",
  "thaBer_to_thaDol_kallisto.merged.gene_counts.tsv",
  "thaCae_to_thaDol_kallisto.merged.gene_counts.tsv",
  "thaDol_kallisto.merged.gene_counts.tsv",
  "thaPal_to_thaDol_kallisto.merged.gene_counts.tsv",
  "thaPel_to_thaDol_kallisto.merged.gene_counts.tsv",
  "thaRuf_to_thaDol_kallisto.merged.gene_counts.tsv",
  "thaTor_to_thaDol_kallisto.merged.gene_counts.tsv"
)

## Give them short run names (similar structure to before)
names(files) <- c(
  "drySqu_to_thaDol",
  "dysSti_to_thaDol",
  "sakCan_to_thaDol",
  "sakCri_to_thaDol",
  "sakLuc_to_thaDol",
  "thaAmb_to_thaDol",
  "thaAtr_to_thaDol",
  "thaBer_to_thaDol",
  "thaCae_to_thaDol",
  "thaDol",
  "thaPal_to_thaDol",
  "thaPel_to_thaDol",
  "thaRuf_to_thaDol",
  "thaTor_to_thaDol"
)

## 1b) Read and merge matrices by gene_id (drop gene_name)
gene_counts_list <- list()

for (run in names(files)) {
  path <- files[[run]]
  if (!file.exists(path)) {
    warning("Missing count file for run: ", run, " at ", path)
    next
  }
  cat("Reading:", path, "for run", run, "\n")
  
  df <- read.delim(path, header = TRUE, check.names = FALSE)
  ## Expect first cols: gene_id, gene_name, then samples
  stopifnot(colnames(df)[1] == "gene_id")
  
  # Drop gene_name, keep gene_id + sample columns
  df <- df[, !colnames(df) %in% "gene_name"]
  
  gene_counts_list[[run]] <- df
}
stopifnot(length(gene_counts_list) > 0)

## Merge on gene_id (analogous to OG before)
merged <- Reduce(function(x, y) merge(x, y, by = "gene_id", all = FALSE),
                 gene_counts_list)
cat("Merged matrix has", nrow(merged), "genes and", ncol(merged) - 1, "samples\n")

rownames(merged) <- merged$gene_id
merged$gene_id <- NULL
cts <- merged   # genes × samples

id <- rownames(cts)

gene_root <- ifelse(
  grepl("_stringtie", id),
  sub("_stringtie.*$", "", id),                 # collapse stringtie isoforms
  sub("^([^\\.]+\\.[^\\.]+).*", "\\1", id)      # collapse XM_... etc to first two fields
)

cts_gene <- rowsum(as.matrix(cts), group = gene_root)

dim(cts)       # transcripts × samples
dim(cts_gene)  # collapsed genes × samples

cts <- cts_gene

## -----------------------------------------
## After:  cts <- cts_gene
## -----------------------------------------

## Build coldata from sample names
samples  <- colnames(cts)
tissue   <- sub(".*_", "", samples)           # after last "_"
id_species <- sub("_[^_]*$", "", samples)     # before last "_"
species  <- sub("[0-9]+.*$", "", id_species)  # leading letters

coldata <- data.frame(
  Sample  = samples,
  Species = species,
  Tissue  = tissue,
  stringsAsFactors = FALSE
)

## 1) Fix drySqu03 tissue mislabels
fix_map <- c(
  "drySqu03_Liver" = "Brain",
  "drySqu03_Heart" = "Liver",
  "drySqu03_Muscle" = "Heart",
  "drySqu03_Kidney" = "Muscle",
  "drySqu03_Brain" = "Kidney"
)

coldata$Tissue <- ifelse(
  coldata$Sample %in% names(fix_map),
  fix_map[coldata$Sample],
  coldata$Tissue
)

## 2) Drop thaPal01_Liver completely (from counts and coldata)
drop_samples <- c("thaPal01_Liver")

keep_cols <- !colnames(cts) %in% drop_samples
cts       <- cts[, keep_cols, drop = FALSE]
coldata   <- coldata[keep_cols, , drop = FALSE]

## Fix thaCap -> thaDol in Species only (as before)
coldata$Species[coldata$Species == "thaCap"] <- "thaDol"

stopifnot(all(colnames(cts) == coldata$Sample))

## -----------------------------------------
## Re-run DESeq2 and PCA from here
## -----------------------------------------

## DESeq2 object
cts_int <- round(as.matrix(cts))
storage.mode(cts_int) <- "integer"

dds <- DESeqDataSetFromMatrix(
  countData = cts_int,
  colData   = coldata,
  design    = ~ Species + Tissue
)

## Filter low-count genes (previously OGs)
keep <- rowSums(counts(dds)) >= 3
dds  <- dds[keep,]

dds$Species <- factor(dds$Species)
dds$Tissue  <- factor(dds$Tissue)

dds <- DESeq(dds)

## Variance stabilisation
vsd     <- vst(dds, blind = FALSE)
vsd_mat <- assay(vsd)

## -----------------------------------------
## Global PCA (PC1×PC2, PC2×PC3, PC1×PC3)
## -----------------------------------------

# Top variable genes
ntop <- 500
rv   <- rowVars(vsd_mat)
select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]

pca_res <- prcomp(t(vsd_mat[select, ]), center = TRUE, scale. = FALSE)
pca_var <- 100 * (pca_res$sdev^2) / sum(pca_res$sdev^2)

pca_df  <- as.data.frame(pca_res$x)
pca_df$Sample <- rownames(pca_df)
pca_df <- merge(pca_df, coldata, by = "Sample")

## Rebuild species color palette (as in your script)
all_species <- sort(unique(pca_df$Species))

warm_species <- c(
  "thaCap",  # will be converted to thaDol
  "sakCri",
  "thaPel",
  "thaTor",
  "thaRuf",
  "sakCan",
  "thaBer"
)
cool_species <- c(
  "sakLuc",
  "drySqu",
  "thaAmb",
  "dysSti",
  "thaCae",
  "thaPal",
  "thaAtr"
)

warm_species[warm_species == "thaCap"] <- "thaDol"
warm_species <- intersect(warm_species, all_species)
cool_species <- intersect(cool_species, all_species)

warm_palette <- c(
  "#b2182b", "#e31a1c", "#fb6a4a",
  "#fdae61", "#ffdd55", "#f46d43", "#d73027"
)[seq_along(warm_species)]

cool_palette <- c(
  "#2166ac", "#4393c3", "#1b7837",
  "#5aae61", "#7b3294", "#c2a5cf", "#7bccc4"
)[seq_along(cool_species)]

pal_warmcool <- c(
  setNames(warm_palette, warm_species),
  setNames(cool_palette, cool_species)
)

missing_species <- setdiff(all_species, names(pal_warmcool))
if (length(missing_species) > 0) {
  pal_warmcool <- c(
    pal_warmcool,
    setNames(rep("#999999", length(missing_species)), missing_species)
  )
}

## PC1–PC2
p12 <- ggplot(pca_df, aes(PC1, PC2, color = Species, shape = Tissue)) +
  geom_point(size = 3, alpha = 0.8) +
  xlab(paste0("PC1: ", round(pca_var[1],1), "% variance")) +
  ylab(paste0("PC2: ", round(pca_var[2],1), "% variance")) +
  coord_fixed() +
  theme_bw() +
  scale_color_manual(values = pal_warmcool, name = "Species")
ggsave("PCA_global_PC1_PC2_no_thaPal01_Liver.png", p12, width = 6, height = 5)
p12

## PC2–PC3
p23 <- ggplot(pca_df, aes(PC2, PC3, color = Species, shape = Tissue)) +
  geom_point(size = 3, alpha = 0.8) +
  xlab(paste0("PC2: ", round(pca_var[2],1), "% variance")) +
  ylab(paste0("PC3: ", round(pca_var[3],1), "% variance")) +
  coord_fixed() +
  theme_bw() +
  scale_color_manual(values = pal_warmcool, name = "Species")
ggsave("PCA_global_PC2_PC3_no_thaPal01_Liver.png", p23, width = 6, height = 5)
p23

## PC1–PC3
p13 <- ggplot(pca_df, aes(PC1, PC3, color = Species, shape = Tissue)) +
  geom_point(size = 3, alpha = 0.8) +
  xlab(paste0("PC1: ", round(pca_var[1],1), "% variance")) +
  ylab(paste0("PC3: ", round(pca_var[3],1), "% variance")) +
  coord_fixed() +
  theme_bw() +
  scale_color_manual(values = pal_warmcool, name = "Species")
ggsave("PCA_global_PC1_PC3_no_thaPal01_Liver.png", p13, width = 6, height = 5)
p13