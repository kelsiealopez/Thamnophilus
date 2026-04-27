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