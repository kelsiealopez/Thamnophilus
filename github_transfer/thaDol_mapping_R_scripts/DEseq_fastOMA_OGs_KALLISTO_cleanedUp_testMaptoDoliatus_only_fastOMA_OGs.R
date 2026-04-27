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

## 1a) List your OG_all9 files
files <- c(
  "drySqu_to_thaDol_kallisto.merged.OG_all9_counts.tsv",
  "dysSti_to_thaDol_kallisto.merged.OG_all9_counts.tsv",
  "sakCan_to_thaDol_kallisto.merged.OG_all9_counts.tsv",
  "sakCri_to_thaDol_kallisto.merged.OG_all9_counts.tsv",
  "sakLuc_to_thaDol_kallisto.merged.OG_all9_counts.tsv",
  "thaAmb_to_thaDol_kallisto.merged.OG_all9_counts.tsv",
  "thaAtr_to_thaDol_kallisto.merged.OG_all9_counts.tsv",
  "thaBer_to_thaDol_kallisto.merged.OG_all9_counts.tsv",
  "thaCae_to_thaDol_kallisto.merged.OG_all9_counts.tsv",
  "thaDol_kallisto.merged.OG_all9_counts.tsv",
  "thaPal_to_thaDol_kallisto.merged.OG_all9_counts.tsv",
  "thaPel_to_thaDol_kallisto.merged.OG_all9_counts.tsv",
  "thaRuf_to_thaDol_kallisto.merged.OG_all9_counts.tsv",
  "thaTor_to_thaDol_kallisto.merged.OG_all9_counts.tsv"
)

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

## 1b) Read and merge matrices by OG (rename to gene_id)
gene_counts_list <- list()

for (run in names(files)) {
  path <- files[[run]]
  if (!file.exists(path)) {
    warning("Missing count file for run: ", run, " at ", path)
    next
  }
  cat("Reading:", path, "for run", run, "\n")
  
  df <- read.delim(path, header = TRUE, check.names = FALSE)
  
  ## Check first column is OG and rename to gene_id
  stopifnot(colnames(df)[1] == "OG")
  colnames(df)[1] <- "gene_id"
  
  ## There is no gene_name column in these files, so don't try to drop it
  ## Just keep gene_id + sample columns:
  # df <- df  # nothing to drop
  
  gene_counts_list[[run]] <- df
}

stopifnot(length(gene_counts_list) > 0)

## Merge on gene_id
merged <- Reduce(function(x, y) merge(x, y, by = "gene_id", all = FALSE),
                 gene_counts_list)

cat("Merged matrix has", nrow(merged), "genes and", ncol(merged) - 1, "samples\n")

rownames(merged) <- merged$gene_id
merged$gene_id <- NULL
cts <- merged

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

## =========================================
## EXCLUDE SAMPLES HERE
## =========================================
samples_to_drop <- c("thaPal01_Liver")

keep_cols <- setdiff(colnames(cts), samples_to_drop)
cts <- cts[, keep_cols, drop = FALSE]
## =========================================







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

## Fix thaCap -> thaDol in Species only
coldata$Species[coldata$Species == "thaCap"] <- "thaDol"

stopifnot(all(colnames(cts) == coldata$Sample))

## DESeq2 object
cts_int <- round(as.matrix(cts))
storage.mode(cts_int) <- "integer"

dds <- DESeqDataSetFromMatrix(
  countData = cts_int,
  colData   = coldata,
  design    = ~ Species + Tissue
)

## Filter low-count OGs
keep <- rowSums(counts(dds)) >= 3
dds  <- dds[keep,]

dds$Species <- factor(dds$Species)
dds$Tissue  <- factor(dds$Tissue)

dds <- DESeq(dds)

## Variance stabilisation
vsd     <- vst(dds, blind = FALSE)
vsd_mat <- assay(vsd)

## ========================================================
## 2) PCA plots (global and tissue-specific)
## ========================================================

## 2.1 Global PCA via DESeq2::plotPCA
pcaData <- plotPCA(vsd, intgroup = c("Species", "Tissue"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

# Warm/cool palette by Species
pca_df <- as.data.frame(prcomp(t(vsd_mat[order(rowVars(vsd_mat), decreasing = TRUE)[1:min(500, nrow(vsd_mat))], ]),
                               center = TRUE, scale. = FALSE)$x)
pca_df$Sample <- rownames(pca_df)
pca_df <- merge(pca_df, coldata, by = "Sample")

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
  "#b2182b",  # deep red
  "#e31a1c",  # bright red
  "#fb6a4a",  # orange‑red
  "#fdae61",  # orange
  "#ffdd55",  # warm yellow
  "#f46d43",  # strong orange‑red
  "#d73027"   # deep orange‑red
)[seq_along(warm_species)]

cool_palette <- c(
  "#2166ac",  # blue
  "#4393c3",  # light blue
  "#1b7837",  # green
  "#5aae61",  # light green
  "#7b3294",  # purple
  "#c2a5cf",  # lavender
  "#7bccc4"   # teal
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

## Use variance from prcomp for PC axes
ntop <- 500
rv   <- rowVars(vsd_mat)
select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
pca_res <- prcomp(t(vsd_mat[select, ]), center = TRUE, scale. = FALSE)
pca_var <- 100 * (pca_res$sdev^2) / sum(pca_res$sdev^2)
pca_df  <- as.data.frame(pca_res$x)
pca_df$Sample <- rownames(pca_df)
pca_df <- merge(pca_df, coldata, by = "Sample")

## Global PC1–PC2
p12 <- ggplot(pca_df, aes(PC1, PC2, color = Species, shape = Tissue)) +
  geom_point(size = 3, alpha = 0.8) +
  xlab(paste0("PC1: ", round(pca_var[1],1), "% variance")) +
  ylab(paste0("PC2: ", round(pca_var[2],1), "% variance")) +
  coord_fixed() +
  theme_bw() +
  scale_color_manual(values = pal_warmcool, name = "Species")
ggsave("PCA_OG_all_species_kallisto_PC1_PC2.png", p12, width = 6, height = 5)
p12


## Global PC2–PC3
p23 <- ggplot(pca_df, aes(PC2, PC3, color = Species, shape = Tissue)) +
  geom_point(size = 3, alpha = 0.8) +
  xlab(paste0("PC2: ", round(pca_var[2],1), "% variance")) +
  ylab(paste0("PC3: ", round(pca_var[3],1), "% variance")) +
  coord_fixed() +
  theme_bw() +
  scale_color_manual(values = pal_warmcool, name = "Species")
ggsave("PCA_OG_all_species_kallisto_PC2_PC3.png", p23, width = 6, height = 5)
p23


## Global PC1–PC3
p13 <- ggplot(pca_df, aes(PC1, PC3, color = Species, shape = Tissue)) +
  geom_point(size = 3, alpha = 0.8) +
  xlab(paste0("PC1: ", round(pca_var[1],1), "% variance")) +
  ylab(paste0("PC3: ", round(pca_var[3],1), "% variance")) +
  coord_fixed() +
  theme_bw() +
  scale_color_manual(values = pal_warmcool, name = "Species")
ggsave("PCA_OG_all_species_kallisto_PC1_PC3.png", p13, width = 6, height = 5)
p13

## Tissue-specific PCA (with labels) for PC1–PC3
tissues_to_run <- c("Brain", "Heart", "Kidney", "Liver", "Muscle")

for (tiss in tissues_to_run) {
  message("Making PCA plots for tissue: ", tiss)
  df_t <- subset(pca_df, Tissue == tiss)
  if (nrow(df_t) == 0) next
  
  ## PC1 vs PC2
  p12_t <- ggplot(df_t, aes(PC1, PC2, color = Species)) +
    geom_point(size = 3, alpha = 0.8, shape = 16) +
    geom_text_repel(aes(label = Species),
                    size = 3,
                    show.legend = FALSE,
                    max.overlaps = Inf,
                    box.padding = 0.3,
                    point.padding = 0.1) +
    xlab(paste0("PC1: ", round(pca_var[1],1), "% variance")) +
    ylab(paste0("PC2: ", round(pca_var[2],1), "% variance")) +
    coord_fixed() +
    theme_bw() +
    scale_color_manual(values = pal_warmcool, name = "Species") +
    ggtitle(paste("PCA –", tiss, "(PC1 vs PC2)"))
  ggsave(paste0("PCA_OG_", tiss, "_kallisto_PC1_PC2_labeled.png"),
         p12_t, width = 6, height = 5)
  
  ## PC2 vs PC3
  p23_t <- ggplot(df_t, aes(PC2, PC3, color = Species)) +
    geom_point(size = 3, alpha = 0.8, shape = 16) +
    geom_text_repel(aes(label = Species),
                    size = 3,
                    show.legend = FALSE,
                    max.overlaps = Inf,
                    box.padding = 0.3,
                    point.padding = 0.1) +
    xlab(paste0("PC2: ", round(pca_var[2],1), "% variance")) +
    ylab(paste0("PC3: ", round(pca_var[3],1), "% variance")) +
    coord_fixed() +
    theme_bw() +
    scale_color_manual(values = pal_warmcool, name = "Species") +
    ggtitle(paste("PCA –", tiss, "(PC2 vs PC3)"))
  ggsave(paste0("PCA_OG_", tiss, "_kallisto_PC2_PC3_labeled.png"),
         p23_t, width = 6, height = 5)
  
  ## PC1 vs PC3
  p13_t <- ggplot(df_t, aes(PC1, PC3, color = Species)) +
    geom_point(size = 3, alpha = 0.8, shape = 16) +
    geom_text_repel(aes(label = Species),
                    size = 3,
                    show.legend = FALSE,
                    max.overlaps = Inf,
                    box.padding = 0.3,
                    point.padding = 0.1) +
    xlab(paste0("PC1: ", round(pca_var[1],1), "% variance")) +
    ylab(paste0("PC3: ", round(pca_var[3],1), "% variance")) +
    coord_fixed() +
    theme_bw() +
    scale_color_manual(values = pal_warmcool, name = "Species") +
    ggtitle(paste("PCA –", tiss, "(PC1 vs PC3)"))
  ggsave(paste0("PCA_OG_", tiss, "_kallisto_PC1_PC3_labeled.png"),
         p13_t, width = 6, height = 5)
}

## ========================================================
## 3) Differential expression for 5 pairs in each tissue
## ========================================================

dds$Species <- factor(dds$Species)
dds$Tissue  <- factor(dds$Tissue)

## Helper: DE for one tissue + species pair
run_DE_tissue_pair <- function(dds_full, tissue, spA, spB,
                               out_prefix = "DE_OG") {
  keep   <- dds_full$Tissue == tissue & dds_full$Species %in% c(spA, spB)
  dds_sub <- dds_full[, keep]
  
  if (ncol(dds_sub) < 2 || length(unique(dds_sub$Species)) < 2) {
    warning("Skipping ", tissue, " ", spA, " vs ", spB,
            ": one or both species have no samples in this tissue.")
    return(NULL)
  }
  
  dds_sub$Species <- factor(dds_sub$Species)
  design(dds_sub) <- ~ Species
  dds_sub <- DESeq(dds_sub)
  
  res <- results(dds_sub, contrast = c("Species", spA, spB))
  
  fname <- paste0(out_prefix, "_", spA, "_vs_", spB, "_", tissue, "_kallisto.csv")
  write.csv(as.data.frame(res), fname)
  message("Wrote: ", fname)
  
  return(res)
}

## 5 species pairs of interest
pairs_of_interest <- tribble(
  ~pair_id,               ~spA,     ~spB,
  "dysSti_vs_sakCri",     "dysSti", "sakCri",
  "sakLuc_vs_sakCan",     "sakLuc", "sakCan",
  "thaAmb_vs_thaPel",     "thaAmb", "thaPel",
  "thaBer_vs_thaAtr",     "thaBer", "thaAtr",
  "thaRuf_vs_thaTor",     "thaRuf", "thaTor"
)

tissues_to_run <- c("Brain", "Heart", "Kidney", "Liver", "Muscle")

## de_sig: list[tissue][[pair_id]] = vector of sig OGs
de_sig <- list()

for (tiss in tissues_to_run) {
  message("\n=== Tissue: ", tiss, " ===")
  de_sig[[tiss]] <- list()
  
  for (i in seq_len(nrow(pairs_of_interest))) {
    pair_id <- pairs_of_interest$pair_id[i]
    spA     <- pairs_of_interest$spA[i]
    spB     <- pairs_of_interest$spB[i]
    
    message("  Running DE for ", pair_id, " (", spA, " vs ", spB, ")")
    res <- run_DE_tissue_pair(
      dds_full   = dds,
      tissue     = tiss,
      spA        = spA,
      spB        = spB,
      out_prefix = "DE_OG_selected"
    )
    if (is.null(res)) next
    
    sig <- as.data.frame(res) %>%
      mutate(OG = rownames(res)) %>%
      filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) >= 1)
    
    de_sig[[tiss]][[pair_id]] <- sig$OG
    message("    ", length(sig$OG), " significant OGs (padj < 0.05, |log2FC| >= 1)")
  }
}

## ========================================================
## 4) Build UpSet input (sig DE OGs, with direction)
## ========================================================

## Helper: read DE table back
read_de_table <- function(tissue, pair_id, spA, spB) {
  fname <- paste0("DE_OG_selected_", spA, "_vs_", spB, "_", tissue, "_kallisto.csv")
  if (!file.exists(fname)) return(NULL)
  df <- read.csv(fname, row.names = 1, check.names = FALSE)
  df$OG <- rownames(df)
  df
}

## Build UpSet-style table for one tissue
## Build UpSet-style table for one tissue
build_tissue_upset <- function(tiss, de_sig, pairs_of_interest) {
  if (!tiss %in% names(de_sig)) return(NULL)
  tissue_list <- de_sig[[tiss]]
  pair_ids    <- names(tissue_list)
  if (length(pair_ids) < 1) return(NULL)
  
  all_ogs <- unique(unlist(tissue_list))
  if (length(all_ogs) == 0) return(NULL)
  
  ## Membership matrix (sig yes/no)
  membership_df <- expand.grid(
    OG      = all_ogs,
    pair_id = pair_ids,
    stringsAsFactors = FALSE
  ) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(present = OG %in% tissue_list[[pair_id]]) %>%
    dplyr::ungroup()
  
  membership_wide <- membership_df %>%
    dplyr::mutate(present = as.integer(present)) %>%
    tidyr::pivot_wider(
      names_from  = pair_id,
      values_from = present,
      values_fill = 0
    )
  
  ## Direction matrix (log2FC sign in each pair)
  dir_list <- list()
  for (i in seq_len(nrow(pairs_of_interest))) {
    pid <- pairs_of_interest$pair_id[i]
    if (!pid %in% pair_ids) next
    
    spA <- pairs_of_interest$spA[i]
    spB <- pairs_of_interest$spB[i]
    
    de_tab <- read_de_table(tiss, pid, spA, spB)
    if (is.null(de_tab)) next
    
    sig_ogs <- tissue_list[[pid]]
    if (is.null(sig_ogs) || length(sig_ogs) == 0) next
    
    de_sub <- de_tab %>%
      dplyr::filter(OG %in% sig_ogs) %>%
      dplyr::select(OG, log2FoldChange) %>%
      dplyr::mutate(
        pair_id   = pid,
        direction = dplyr::case_when(
          log2FoldChange > 0 ~ "up_in_spA",
          log2FoldChange < 0 ~ "down_in_spA",
          TRUE               ~ "zero"
        )
      )
    dir_list[[pid]] <- de_sub
  }
  
  if (length(dir_list) == 0) {
    upset_df <- membership_wide
  } else {
    dir_df <- dplyr::bind_rows(dir_list)
    dir_wide <- dir_df %>%
      dplyr::select(OG, pair_id, direction) %>%
      tidyr::pivot_wider(
        names_from  = pair_id,
        values_from = direction,
        values_fill = NA
      )
    
    ## direction columns get "_dir" suffix
    names(dir_wide)[-1] <- paste0(names(dir_wide)[-1], "_dir")
    
    upset_df <- membership_wide %>%
      dplyr::left_join(dir_wide, by = "OG")
  }
  
  upset_df
}

## Rebuild tissue_upset
tissue_upset <- lapply(
  names(de_sig),
  build_tissue_upset,
  de_sig = de_sig,
  pairs_of_interest = pairs_of_interest
)
names(tissue_upset) <- names(de_sig)

## UpSet plots for each tissue (sig DE OGs only)
for (tiss in names(tissue_upset)) {
  df_t <- tissue_upset[[tiss]]
  if (is.null(df_t)) next
  
  set_cols <- names(df_t)[names(df_t) %in% pairs_of_interest$pair_id]
  if (length(set_cols) == 0) next
  
  message("Drawing UpSet (sig DE) for tissue: ", tiss)
  pdf(paste0("UpSet_SIG_DE_OGs_", tiss, ".pdf"), width = 7, height = 5)
  print(upset(
    df_t,
    set_cols,
    name = "Sig DE OGs",
    min_size = 10
  ) + ggtitle(paste("Significant DE OGs –", tiss)))
  dev.off()
}

## ========================================================
## 5) OGs DE in >=2 (Kidney) / >=3 (others) pairs (overlap)
## ========================================================

og_multi_tissue <- list()

for (tiss in names(de_sig)) {
  tissue_list <- de_sig[[tiss]]
  if (is.null(tissue_list) || length(tissue_list) == 0) next
  
  og_counts <- table(unlist(tissue_list))
  og_counts <- sort(og_counts, decreasing = TRUE)
  
  thresh <- if (tiss == "Kidney") 2L else 3L
  og_sel <- names(og_counts[og_counts >= thresh])
  if (length(og_sel) == 0) next
  
  pair_ids <- names(tissue_list)
  og_to_pairs <- lapply(og_sel, function(og) {
    names(Filter(function(v) og %in% v, tissue_list))
  })
  
  df_t <- tibble(
    Tissue  = tiss,
    OG      = og_sel,
    n_pairs = as.integer(og_counts[og_sel]),
    Pairs   = vapply(og_to_pairs, paste, collapse = ";",
                     FUN.VALUE = character(1))
  )
  og_multi_tissue[[tiss]] <- df_t
  
  fn <- paste0("OGs_DE_in_multiple_pairs_", tiss, ".csv")
  write.csv(df_t, fn, row.names = FALSE)
  message("Wrote: ", fn)
}

og_multi_all <- bind_rows(og_multi_tissue)
write.csv(og_multi_all,
          "OGs_DE_in_multiple_pairs_ALL_tissues.csv",
          row.names = FALSE)

## ========================================================
## 6) Map OGs to gene IDs and choose representative gene name
## ========================================================

## Map OG -> gene_id via fastOMA output
og_map_file <- "/n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/fastoma/out_folder_9_species_mamba/OG_gene_map_all9.tsv"
og_map <- read.delim(og_map_file, header = TRUE, stringsAsFactors = FALSE)  # columns: OG, species, gene_id

og_multi_all_mapped <- og_multi_all %>%
  left_join(og_map, by = "OG") %>%
  arrange(Tissue, OG, species)

write.csv(og_multi_all_mapped,
          "OGs_DE_in_multiple_pairs_ALL_tissues_with_geneIDs.csv",
          row.names = FALSE)

## Load isoforms for chicken and zebra finch to help infer gene symbols
gal_iso_path <- "/n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/annotation/02_toga_galGal/sakCri/reference/isoforms.tsv"
tae_iso_path <- "/n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/annotation/02_toga_taeGut/sakCri/reference/isoforms.tsv"

gal_iso <- read.delim(gal_iso_path, header = TRUE, stringsAsFactors = FALSE)
tae_iso <- read.delim(tae_iso_path, header = TRUE, stringsAsFactors = FALSE)

gal_map <- gal_iso %>%
  rename(galGal_geneID = geneID,
         transcriptID  = transcriptID) %>%
  distinct(transcriptID, galGal_geneID, .keep_all = TRUE)

tae_map <- tae_iso %>%
  rename(taeGut_geneID = geneID,
         transcriptID  = transcriptID) %>%
  distinct(transcriptID, taeGut_geneID, .keep_all = TRUE)

clean_transcript_root <- function(x) {
  sub("^([^\\.]+\\.[^\\.]+).*", "\\1", x)
}

og_annot_ref <- og_multi_all_mapped %>%
  mutate(transcript_root = clean_transcript_root(gene_id)) %>%
  left_join(gal_map, by = c("transcript_root" = "transcriptID")) %>%
  left_join(tae_map, by = c("transcript_root" = "transcriptID")) %>%
  arrange(Tissue, OG, species)

write.csv(og_annot_ref,
          "OGs_DE_in_multiple_pairs_ALL_tissues_with_refGeneIDs.csv",
          row.names = FALSE)

#

## Choose representative gene name per (Tissue, OG)
choose_rep_gene <- function(df_group) {
  # 1) prefer taeGut symbol
  tg <- unique(df_group$taeGut_geneID[!is.na(df_group$taeGut_geneID) &
                                        df_group$taeGut_geneID != ""])
  if (length(tg) >= 1) {
    return(data.frame(gene_name = tg[1], name_source = "taeGut_geneID",
                      stringsAsFactors = FALSE))
  }
  # 2) else galGal
  gg <- unique(df_group$galGal_geneID[!is.na(df_group$galGal_geneID) &
                                        df_group$galGal_geneID != ""])
  if (length(gg) >= 1) {
    return(data.frame(gene_name = gg[1], name_source = "galGal_geneID",
                      stringsAsFactors = FALSE))
  }
  # 3) else parse from gene_id (stringtie-style)
  ids    <- unique(df_group$gene_id)
  ids_st <- ids[grepl("stringtie", ids)]
  if (length(ids_st) >= 1) {
    id0  <- ids_st[1]
    base <- sub("_stringtie.*$", "", id0)
    base <- sub("-like$", "", base)
    base <- sub("^([^\\.]+\\.[^\\.]+).*", "\\1", base)
    return(data.frame(gene_name = base, name_source = "parsed_from_gene_id",
                      stringsAsFactors = FALSE))
  }
  # 4) else use transcript_root
  tr <- unique(df_group$transcript_root)
  tr <- tr[!is.na(tr) & tr != ""]
  if (length(tr) >= 1) {
    return(data.frame(gene_name = tr[1], name_source = "transcript_root",
                      stringsAsFactors = FALSE))
  }
  # 5) unknown
  data.frame(gene_name = NA_character_, name_source = "none",
             stringsAsFactors = FALSE)
}

og_meta <- og_annot_ref %>%
  dplyr::select(Tissue, OG, n_pairs, Pairs) %>%
  dplyr::distinct()

rep_names <- og_annot_ref %>%
  dplyr::group_by(Tissue, OG) %>%
  dplyr::do(choose_rep_gene(.)) %>%
  dplyr::ungroup()

og_rep <- og_meta %>%
  dplyr::left_join(rep_names, by = c("Tissue", "OG")) %>%
  dplyr::arrange(Tissue, OG)

write.csv(
  og_rep,
  "OGs_DE_in_multiple_pairs_ALL_tissues_repGeneName.csv",
  row.names = FALSE
)

for (tiss in unique(og_rep$Tissue)) {
  df_t <- og_rep %>%
    dplyr::filter(Tissue == tiss)
  fn <- paste0("OGs_DE_in_multiple_pairs_", tiss, "_repGeneName.csv")
  write.csv(df_t, fn, row.names = FALSE)
}

## ========================================================
## 7) Convergent dry vs mesic directions (OGs)
## ========================================================

pairs_arid <- tribble(
  ~pair_id,             ~spA,     ~spB,     ~dry_species, ~mesic_species,
  "dysSti_vs_sakCri",   "dysSti", "sakCri", "sakCri",     "dysSti",
  "sakLuc_vs_sakCan",   "sakLuc", "sakCan", "sakCan",     "sakLuc",
  "thaAmb_vs_thaPel",   "thaAmb", "thaPel", "thaPel",     "thaAmb",
  "thaBer_vs_thaAtr",   "thaBer", "thaAtr", "thaBer",     "thaAtr",
  "thaRuf_vs_thaTor",   "thaRuf", "thaTor", "thaTor",     "thaRuf"
)

## For each tissue, get direction relative to dry species
get_dry_direction_for_tissue <- function(tiss, de_sig, pairs_arid) {
  if (!tiss %in% names(de_sig)) return(NULL)
  tissue_list <- de_sig[[tiss]]
  
  dir_list <- list()
  
  for (i in seq_len(nrow(pairs_arid))) {
    pid <- pairs_arid$pair_id[i]
    spA <- pairs_arid$spA[i]
    spB <- pairs_arid$spB[i]
    dry <- pairs_arid$dry_species[i]
    
    if (!pid %in% names(tissue_list)) next
    
    sig_ogs <- tissue_list[[pid]]
    if (is.null(sig_ogs) || length(sig_ogs) == 0) next
    
    fname <- paste0("DE_OG_selected_", spA, "_vs_", spB, "_", tiss, "_kallisto.csv")
    if (!file.exists(fname)) next
    
    de_tab <- read.csv(fname, row.names = 1, check.names = FALSE)
    de_tab$OG <- rownames(de_tab)
    
    de_sub <- de_tab[de_tab$OG %in% sig_ogs, c("OG", "log2FoldChange")]
    if (nrow(de_sub) == 0) next
    
    dir_rel <- if (dry == spA) {
      ifelse(de_sub$log2FoldChange > 0, "dry_up",
             ifelse(de_sub$log2FoldChange < 0, "dry_down", "zero"))
    } else if (dry == spB) {
      ifelse(de_sub$log2FoldChange > 0, "dry_down",
             ifelse(de_sub$log2FoldChange < 0, "dry_up", "zero"))
    } else {
      stop("dry_species not equal to spA or spB for pair ", pid)
    }
    
    dir_list[[pid]] <- data.frame(
      Tissue       = tiss,
      OG           = de_sub$OG,
      pair_id      = pid,
      directionDry = dir_rel,
      stringsAsFactors = FALSE
    )
  }
  
  if (length(dir_list) == 0) return(NULL)
  do.call(rbind, dir_list)
}

tissues_to_run <- c("Brain", "Heart", "Kidney", "Liver", "Muscle")
dry_dir_all <- list()

for (tiss in tissues_to_run) {
  message("Computing dry-relative directions for tissue: ", tiss)
  dry_dir_all[[tiss]] <- get_dry_direction_for_tissue(tiss, de_sig = de_sig, pairs_arid = pairs_arid)
}

## Identify convergent dry-up / dry-down OGs per tissue
conv_ogs <- list()

for (tiss in tissues_to_run) {
  df_dir <- dry_dir_all[[tiss]]
  if (is.null(df_dir) || nrow(df_dir) == 0) next
  
  df_ns <- df_dir %>% filter(directionDry %in% c("dry_up", "dry_down"))
  if (nrow(df_ns) == 0) next
  
  counts <- df_ns %>%
    group_by(OG) %>%
    summarise(
      n_dry_up   = sum(directionDry == "dry_up"),
      n_dry_down = sum(directionDry == "dry_down"),
      n_pairs    = n(),
      .groups    = "drop"
    )
  
  thresh <- if (tiss == "Kidney") 2L else 3L
  
  og_up_conv   <- counts$OG[counts$n_dry_up   >= thresh]
  og_down_conv <- counts$OG[counts$n_dry_down >= thresh]
  
  conv_ogs[[tiss]] <- list(
    dry_up   = og_up_conv,
    dry_down = og_down_conv,
    counts   = counts
  )
  
  message("Tissue ", tiss, ": ",
          length(og_up_conv), " OGs convergently UP in dry (>= ", thresh, " pairs); ",
          length(og_down_conv), " OGs convergently DOWN in dry.")
}

## ========================================================
## 8) GO enrichment (human only) for convergent dry genes
## ========================================================

## Helper: GO for one tissue / OG subset (UPdry or DOWNdry, human only)
run_go_for_tissue_dry <- function(tiss,
                                  og_subset,
                                  OrgDb = org.Hs.eg.db,
                                  ont = "BP",
                                  pcut = 0.05,
                                  dir_tag = c("UPdry","DOWNdry")) {
  dir_tag <- match.arg(dir_tag)
  
  fn_rep <- paste0("OGs_DE_in_multiple_pairs_", tiss, "_repGeneName.csv")
  if (!file.exists(fn_rep)) {
    warning("Rep gene file not found for tissue ", tiss, ": ", fn_rep)
    return(NULL)
  }
  df_rep <- read.csv(fn_rep, stringsAsFactors = FALSE)
  if (!all(c("OG","gene_name") %in% names(df_rep))) {
    warning("No OG/gene_name columns in ", fn_rep)
    return(NULL)
  }
  
  if (!is.null(og_subset)) {
    df_rep <- df_rep %>% filter(OG %in% og_subset)
  }
  
  gene_list <- unique(na.omit(df_rep$gene_name))
  message("Tissue ", tiss, " – GO (", dir_tag, ", human) with ",
          length(gene_list), " unique gene_name entries.")
  
  if (length(gene_list) < 5) {
    warning("Too few genes (", length(gene_list),
            ") in tissue ", tiss, " for GO (", dir_tag, "); skipping.")
    return(NULL)
  }
  
  conv <- tryCatch(
    bitr(gene_list,
         fromType = "SYMBOL",
         toType   = "ENTREZID",
         OrgDb    = OrgDb),
    error = function(e) {
      warning("bitr conversion failed for human in ", tiss, " (", dir_tag, "): ", e$message)
      return(NULL)
    }
  )
  if (is.null(conv) || nrow(conv) == 0) {
    warning("No symbols converted to Entrez for human in ", tiss, " (", dir_tag, ")")
    return(NULL)
  }
  
  ego <- tryCatch(
    enrichGO(
      gene          = conv$ENTREZID,
      OrgDb         = OrgDb,
      keyType       = "ENTREZID",
      ont           = ont,
      pAdjustMethod = "BH",
      pvalueCutoff  = pcut,
      qvalueCutoff  = pcut,
      readable      = TRUE
    ),
    error = function(e) {
      warning("enrichGO failed for human in ", tiss, " (", dir_tag, "): ", e$message)
      return(NULL)
    }
  )
  
  if (is.null(ego) || nrow(as.data.frame(ego)) == 0) {
    warning("No enriched GO terms for human in ", tiss, " (", dir_tag, ")")
    return(NULL)
  }
  
  ego_df <- as.data.frame(ego)
  tag <- paste0("DRY_", dir_tag)  # e.g. DRY_UPdry
  
  out_tab <- paste0("GO_BP_Hs_", tiss, "_", tag, "_enrichment.csv")
  write.csv(ego_df, out_tab, row.names = FALSE)
  message("  Wrote GO table: ", out_tab)
  
  p_dot <- dotplot(ego, showCategory = 20) +
    ggtitle(paste("GO BP –", tiss, "(human,", tag, ")"))
  out_dot <- paste0("GO_BP_Hs_", tiss, "_", tag, "_dotplot.png")
  ggsave(out_dot, p_dot, width = 8, height = 6, dpi = 300)
  message("  Wrote dotplot: ", out_dot)
  
  ego
}

## Run GO for convergent dry-up / dry-down in each tissue (human only)
for (tiss in c("Brain", "Heart", "Kidney", "Liver", "Muscle")) {
  og_up   <- conv_ogs[[tiss]]$dry_up
  og_down <- conv_ogs[[tiss]]$dry_down
  
  if (!is.null(og_up) && length(og_up) >= 5) {
    message("\n=== ", tiss, " – convergently UP in dry (human GO) ===")
    run_go_for_tissue_dry(tiss, og_up, dir_tag = "UPdry")
  }
  
  if (!is.null(og_down) && length(og_down) >= 5) {
    message("\n=== ", tiss, " – convergently DOWN in dry (human GO) ===")
    run_go_for_tissue_dry(tiss, og_down, dir_tag = "DOWNdry")
  }
}

############################################################
## 9) Convenience helpers to view plots in RStudio
############################################################

## -------------------------------
## 9.1 PCA plots in the console
## -------------------------------
## Uses: pca_df, pca_var, pal_warmcool (already defined above)

view_pca_global <- function() {
  # Global PC1–PC2
  p12 <- ggplot(pca_df, aes(PC1, PC2, color = Species, shape = Tissue)) +
    geom_point(size = 3, alpha = 0.8) +
    xlab(paste0("PC1: ", round(pca_var[1],1), "% variance")) +
    ylab(paste0("PC2: ", round(pca_var[2],1), "% variance")) +
    coord_fixed() +
    theme_bw() +
    scale_color_manual(values = pal_warmcool, name = "Species") +
    ggtitle("Global PCA (PC1 vs PC2)")
  print(p12)
  
  # Global PC2–PC3
  p23 <- ggplot(pca_df, aes(PC2, PC3, color = Species, shape = Tissue)) +
    geom_point(size = 3, alpha = 0.8) +
    xlab(paste0("PC2: ", round(pca_var[2],1), "% variance")) +
    ylab(paste0("PC3: ", round(pca_var[3],1), "% variance")) +
    coord_fixed() +
    theme_bw() +
    scale_color_manual(values = pal_warmcool, name = "Species") +
    ggtitle("Global PCA (PC2 vs PC3)")
  print(p23)
  
  # Global PC1–PC3
  p13 <- ggplot(pca_df, aes(PC1, PC3, color = Species, shape = Tissue)) +
    geom_point(size = 3, alpha = 0.8) +
    xlab(paste0("PC1: ", round(pca_var[1],1), "% variance")) +
    ylab(paste0("PC3: ", round(pca_var[3],1), "% variance")) +
    coord_fixed() +
    theme_bw() +
    scale_color_manual(values = pal_warmcool, name = "Species") +
    ggtitle("Global PCA (PC1 vs PC3)")
  print(p13)
}

# View all global PCA plots:
# view_pca_global()

view_pca_tissue <- function(tiss) {
  df_t <- subset(pca_df, Tissue == tiss)
  if (nrow(df_t) == 0) {
    warning("No samples for tissue: ", tiss)
    return(invisible(NULL))
  }
  
  message("Showing PCA for tissue: ", tiss)
  
  # PC1 vs PC2
  p12_t <- ggplot(df_t, aes(PC1, PC2, color = Species)) +
    geom_point(size = 3, alpha = 0.8, shape = 16) +
    ggrepel::geom_text_repel(
      aes(label = Species),
      size = 3,
      show.legend = FALSE,
      max.overlaps = Inf,
      box.padding = 0.3,
      point.padding = 0.1
    ) +
    xlab(paste0("PC1: ", round(pca_var[1],1), "% variance")) +
    ylab(paste0("PC2: ", round(pca_var[2],1), "% variance")) +
    coord_fixed() +
    theme_bw() +
    scale_color_manual(values = pal_warmcool, name = "Species") +
    ggtitle(paste("PCA –", tiss, "(PC1 vs PC2)"))
  print(p12_t)
  
  # PC2 vs PC3
  p23_t <- ggplot(df_t, aes(PC2, PC3, color = Species)) +
    geom_point(size = 3, alpha = 0.8, shape = 16) +
    ggrepel::geom_text_repel(
      aes(label = Species),
      size = 3,
      show.legend = FALSE,
      max.overlaps = Inf,
      box.padding = 0.3,
      point.padding = 0.1
    ) +
    xlab(paste0("PC2: ", round(pca_var[2],1), "% variance")) +
    ylab(paste0("PC3: ", round(pca_var[3],1), "% variance")) +
    coord_fixed() +
    theme_bw() +
    scale_color_manual(values = pal_warmcool, name = "Species") +
    ggtitle(paste("PCA –", tiss, "(PC2 vs PC3)"))
  print(p23_t)
  
  # PC1 vs PC3
  p13_t <- ggplot(df_t, aes(PC1, PC3, color = Species)) +
    geom_point(size = 3, alpha = 0.8, shape = 16) +
    ggrepel::geom_text_repel(
      aes(label = Species),
      size = 3,
      show.legend = FALSE,
      max.overlaps = Inf,
      box.padding = 0.3,
      point.padding = 0.1
    ) +
    xlab(paste0("PC1: ", round(pca_var[1],1), "% variance")) +
    ylab(paste0("PC3: ", round(pca_var[3],1), "% variance")) +
    coord_fixed() +
    theme_bw() +
    scale_color_manual(values = pal_warmcool, name = "Species") +
    ggtitle(paste("PCA –", tiss, "(PC1 vs PC3)"))
  print(p13_t)
}

# Examples:
view_pca_tissue("Kidney")
view_pca_tissue("Heart")
view_pca_tissue("Muscle")
view_pca_tissue("Brain")
view_pca_tissue("Liver")


## -------------------------------
## 9.2 UpSet plots in the console
## -------------------------------
## Uses: tissue_upset, pairs_of_interest

view_upset_tissue <- function(tiss, min_size = 10) {
  df_t <- tissue_upset[[tiss]]
  if (is.null(df_t)) {
    warning("No UpSet data for tissue: ", tiss)
    return(invisible(NULL))
  }
  
  set_cols <- names(df_t)[names(df_t) %in% pairs_of_interest$pair_id]
  if (length(set_cols) == 0) {
    warning("No set columns for tissue: ", tiss)
    return(invisible(NULL))
  }
  
  message("Showing UpSet for tissue: ", tiss)
  p <- upset(
    df_t,
    set_cols,
    name = "Sig DE OGs",
    min_size = min_size
  ) + ggtitle(paste("Significant DE OGs –", tiss))
  print(p)
}

# Examples:
# view_upset_tissue("Kidney")
# view_upset_tissue("Brain")


## -------------------------------
## 9.3 GO dotplots (re-run, human only)
## -------------------------------
## Uses: conv_ogs, run_go_for_tissue_dry

view_go_tissue <- function(tiss,
                           dir_tag = c("UPdry","DOWNdry"),
                           ont = "BP",
                           pcut = 0.05) {
  dir_tag <- match.arg(dir_tag)
  
  og_subset <- if (dir_tag == "UPdry") {
    conv_ogs[[tiss]]$dry_up
  } else {
    conv_ogs[[tiss]]$dry_down
  }
  
  if (is.null(og_subset) || length(og_subset) < 5) {
    warning("Not enough convergent OGs for ", tiss, " (", dir_tag, ")")
    return(invisible(NULL))
  }
  
  message("Re-running GO for ", tiss, " (", dir_tag, ", human only)")
  ego <- run_go_for_tissue_dry(
    tiss      = tiss,
    og_subset = og_subset,
    OrgDb     = org.Hs.eg.db,
    ont       = ont,
    pcut      = pcut,
    dir_tag   = dir_tag
  )
  # run_go_for_tissue_dry already prints a dotplot and saves PNG;
  # this function just returns the enrichResult object for inspection.
  invisible(ego)
}

# Examples:
# view_go_tissue("Kidney", "UPdry")
# view_go_tissue("Kidney", "DOWNdry")
# view_go_tissue("Heart", "UPdry")


############################################################
## How to use interactively in RStudio
############################################################

# 1) Global PCA:
view_pca_global()

# 2) PCA for a given tissue:
view_pca_tissue("Kidney")
view_pca_tissue("Heart")

# 3) UpSet plot for a tissue:
view_upset_tissue("Kidney")
view_upset_tissue("Brain")
view_upset_tissue("Muscle")
view_upset_tissue("Heart")
view_upset_tissue("Liver")


# 4) GO dotplots (human) for convergent dry-up / dry-down:
view_go_tissue("Kidney", "UPdry")
view_go_tissue("Kidney", "DOWNdry")
view_go_tissue("Heart",  "UPdry")
view_go_tissue("Heart",  "DOWNdry")
view_go_tissue("Muscle",  "UPdry")
view_go_tissue("Muscle",  "DOWNdry")
view_go_tissue("Brain",  "UPdry")
view_go_tissue("Brain",  "DOWNdry")
view_go_tissue("Liver",  "UPdry")
view_go_tissue("Liver",  "DOWNdry")




## ========================================================
## 3b) Venn / Euler diagrams of DE OG overlap per tissue
## ========================================================

install.packages("eulerr")  # run once if needed
library(eulerr)

# Function to make Venn/Euler for one tissue
plot_venn_for_tissue <- function(tiss, de_sig_t, pairs_of_interest) {
  if (is.null(de_sig_t) || length(de_sig_t) == 0) {
    warning("No DE sets for tissue: ", tiss)
    return(invisible(NULL))
  }
  
  # Keep only pairs that actually have DE OGs
  sets <- de_sig_t
  sets <- sets[!vapply(sets, is.null, logical(1))]
  sets <- sets[vapply(sets, function(x) length(x) > 0, logical(1))]
  
  if (length(sets) < 2) {
    warning("Fewer than 2 non-empty DE sets for tissue: ", tiss)
    return(invisible(NULL))
  }
  
  # eulerr expects a named list of vectors
  # Optionally sort sets by name to keep order consistent
  sets <- sets[sort(names(sets))]
  
  fit <- euler(sets)
  
  p <- plot(
    fit,
    fills = list(fill = RColorBrewer::brewer.pal(max(3, length(sets)), "Set3"), alpha = 0.6),
    labels = list(font = 2),
    edges = list(lwd = 1),
    quantities = list(type = "counts", font = 2),
    main = paste("DE OG overlap –", tiss)
  )
  
  # Show in plot pane
  print(p)
  
  # Save to file
  out_png <- paste0("Venn_DE_OGs_", tiss, "_pairs.png")
  png(out_png, width = 800, height = 800, res = 150)
  plot(
    fit,
    fills = list(fill = RColorBrewer::brewer.pal(max(3, length(sets)), "Set3"), alpha = 0.6),
    labels = list(font = 2),
    edges = list(lwd = 1),
    quantities = list(type = "counts", font = 2),
    main = paste("DE OG overlap –", tiss)
  )
  dev.off()
  message("Wrote Venn/Euler PNG: ", out_png)
  
  invisible(fit)
}

# Make Venn/Euler per tissue
for (tiss in names(de_sig)) {
  message("\n--- Venn/Euler for tissue: ", tiss, " ---")
  plot_venn_for_tissue(tiss, de_sig[[tiss]], pairs_of_interest)
}





###
#
###
#
###
#
###

## Helper: GO for one tissue / OG subset (UPdry or DOWNdry, human only)
run_go_for_tissue_dry <- function(tiss,
                                  og_subset,
                                  OrgDb = org.Hs.eg.db,
                                  ont = "BP",
                                  pcut = 0.05,
                                  dir_tag = c("UPdry","DOWNdry")) {
  dir_tag <- match.arg(dir_tag)
  
  fn_rep <- paste0("OGs_DE_in_multiple_pairs_", tiss, "_repGeneName.csv")
  if (!file.exists(fn_rep)) {
    warning("Rep gene file not found for tissue ", tiss, ": ", fn_rep)
    return(NULL)
  }
  df_rep <- read.csv(fn_rep, stringsAsFactors = FALSE)
  if (!all(c("OG","gene_name") %in% names(df_rep))) {
    warning("No OG/gene_name columns in ", fn_rep)
    return(NULL)
  }
  
  if (!is.null(og_subset)) {
    df_rep <- df_rep %>% dplyr::filter(OG %in% og_subset)
  }
  
  gene_list <- unique(na.omit(df_rep$gene_name))
  message("Tissue ", tiss, " – GO (", dir_tag, ", human) with ",
          length(gene_list), " unique gene_name entries.")
  
  if (length(gene_list) < 5) {
    warning("Too few genes (", length(gene_list),
            ") in tissue ", tiss, " for GO (", dir_tag, "); skipping.")
    return(NULL)
  }
  
  conv <- tryCatch(
    bitr(gene_list,
         fromType = "SYMBOL",
         toType   = "ENTREZID",
         OrgDb    = OrgDb),
    error = function(e) {
      warning("bitr conversion failed for human in ", tiss, " (", dir_tag, "): ", e$message)
      return(NULL)
    }
  )
  if (is.null(conv) || nrow(conv) == 0) {
    warning("No symbols converted to Entrez for human in ", tiss, " (", dir_tag, ")")
    return(NULL)
  }
  
  ego <- tryCatch(
    enrichGO(
      gene          = conv$ENTREZID,
      OrgDb         = OrgDb,
      keyType       = "ENTREZID",
      ont           = ont,
      pAdjustMethod = "BH",
      pvalueCutoff  = pcut,
      qvalueCutoff  = pcut,
      readable      = TRUE
    ),
    error = function(e) {
      warning("enrichGO failed for human in ", tiss, " (", dir_tag, "): ", e$message)
      return(NULL)
    }
  )
  
  if (is.null(ego) || nrow(as.data.frame(ego)) == 0) {
    warning("No enriched GO terms for human in ", tiss, " (", dir_tag, ")")
    return(NULL)
  }
  
  ego_df <- as.data.frame(ego)
  tag <- paste0("DRY_", dir_tag)  # e.g. DRY_UPdry
  
  out_tab <- paste0("GO_BP_Hs_", tiss, "_", tag, "_enrichment.csv")
  write.csv(ego_df, out_tab, row.names = FALSE)
  message("  Wrote GO table: ", out_tab)
  
  # Dotplot with smaller fonts
  p_dot <- dotplot(ego, showCategory = 20) +
    ggtitle(paste("GO BP –", tiss, "(human,", tag, ")")) +
    theme_bw() +
    theme(
      plot.title   = element_text(size = 11),
      axis.text.x  = element_text(size = 8),
      axis.text.y  = element_text(size = 8),
      axis.title.x = element_text(size = 9),
      axis.title.y = element_text(size = 9)
    )
  
  # Show in plot pane
  print(p_dot)
  
  # Save PNG
  out_dot <- paste0("GO_BP_Hs_", tiss, "_", tag, "_dotplot.png")
  ggsave(out_dot, p_dot, width = 8, height = 6, dpi = 300)
  message("  Wrote dotplot: ", out_dot)
  
  ego
}



view_go_tissue("Kidney", "UPdry")
view_go_tissue("Kidney", "DOWNdry")
view_go_tissue("Liver", "UPdry")
view_go_tissue("Liver", "DOWNdry")
view_go_tissue("Muscle", "UPdry")
view_go_tissue("Muscle", "DOWNdry")
view_go_tissue("Brain", "UPdry")
view_go_tissue("Brain", "DOWNdry")
view_go_tissue("Heart", "UPdry")
view_go_tissue("Heart", "DOWNdry")



#
#
#
# make it one big go analysis 

## ========================================================
## 10) Global GO enrichment: all convergent dry OGs combined (human only)
## ========================================================

run_global_go_all_dry <- function(OrgDb = org.Hs.eg.db,
                                  ont   = "BP",
                                  pcut  = 0.05) {
  tissues <- c("Brain", "Heart", "Kidney", "Liver", "Muscle")
  
  # 1) Collect all convergent OGs across tissues & directions
  all_conv_ogs <- c()
  for (tiss in tissues) {
    if (is.null(conv_ogs[[tiss]])) next
    og_up   <- conv_ogs[[tiss]]$dry_up
    og_down <- conv_ogs[[tiss]]$dry_down
    all_conv_ogs <- c(all_conv_ogs, og_up, og_down)
  }
  all_conv_ogs <- unique(na.omit(all_conv_ogs))
  
  if (length(all_conv_ogs) < 5) {
    warning("Too few convergent OGs overall (", length(all_conv_ogs),
            ") to run a global GO analysis.")
    return(NULL)
  }
  
  message("Global GO: combining ", length(all_conv_ogs),
          " convergent OGs across all tissues/directions.")
  
  # 2) Load repGeneName tables and collect gene_name for these OGs
  all_rep <- list()
  for (tiss in tissues) {
    fn_rep <- paste0("OGs_DE_in_multiple_pairs_", tiss, "_repGeneName.csv")
    if (!file.exists(fn_rep)) {
      warning("Rep gene file not found for tissue ", tiss, ": ", fn_rep)
      next
    }
    df_rep <- read.csv(fn_rep, stringsAsFactors = FALSE)
    if (!all(c("OG","gene_name") %in% names(df_rep))) {
      warning("No OG/gene_name columns in ", fn_rep)
      next
    }
    all_rep[[tiss]] <- df_rep
  }
  if (length(all_rep) == 0) {
    warning("No repGeneName files could be read; aborting global GO.")
    return(NULL)
  }
  
  rep_all <- dplyr::bind_rows(all_rep)
  rep_all <- rep_all %>%
    dplyr::filter(OG %in% all_conv_ogs)
  
  gene_list <- unique(na.omit(rep_all$gene_name))
  message("Global GO: ", length(gene_list), " unique gene_name entries.")
  
  if (length(gene_list) < 5) {
    warning("Too few gene symbols (", length(gene_list),
            ") after mapping OGs; aborting global GO.")
    return(NULL)
  }
  
  # 3) SYMBOL -> ENTREZ
  conv <- tryCatch(
    bitr(gene_list,
         fromType = "SYMBOL",
         toType   = "ENTREZID",
         OrgDb    = OrgDb),
    error = function(e) {
      warning("Global GO: bitr conversion failed: ", e$message)
      return(NULL)
    }
  )
  if (is.null(conv) || nrow(conv) == 0) {
    warning("Global GO: no symbols converted to Entrez.")
    return(NULL)
  }
  
  # 4) Run enrichGO
  ego <- tryCatch(
    enrichGO(
      gene          = conv$ENTREZID,
      OrgDb         = OrgDb,
      keyType       = "ENTREZID",
      ont           = ont,
      pAdjustMethod = "BH",
      pvalueCutoff  = pcut,
      qvalueCutoff  = pcut,
      readable      = TRUE
    ),
    error = function(e) {
      warning("Global GO: enrichGO failed: ", e$message)
      return(NULL)
    }
  )
  
  if (is.null(ego) || nrow(as.data.frame(ego)) == 0) {
    warning("Global GO: no enriched GO terms.")
    return(NULL)
  }
  
  ego_df <- as.data.frame(ego)
  
  out_tab <- "GO_BP_Hs_ALL_tissues_DRY_global_enrichment.csv"
  write.csv(ego_df, out_tab, row.names = FALSE)
  message("Global GO: wrote table ", out_tab)
  
  # 5) Dotplot with smaller fonts
  p_dot <- dotplot(ego, showCategory = 20) +
    ggtitle("GO BP – ALL tissues (human, convergent dry)") +
    theme_bw() +
    theme(
      plot.title   = element_text(size = 11),
      axis.text.x  = element_text(size = 8),
      axis.text.y  = element_text(size = 8),
      axis.title.x = element_text(size = 9),
      axis.title.y = element_text(size = 9)
    )
  
  print(p_dot)
  
  out_dot <- "GO_BP_Hs_ALL_tissues_DRY_global_dotplot.png"
  ggsave(out_dot, p_dot, width = 8, height = 6, dpi = 300)
  message("Global GO: wrote dotplot ", out_dot)
  
  invisible(ego)
}

## Run the global GO once after all conv_ogs and repGeneName files are ready:
ego_global <- run_global_go_all_dry()





# # test label swap
# ## Fix in coldata
# fix_map <- c(
#   "drySqu03_Liver" = "Brain",
#   "drySqu03_Heart" = "Liver",
#   "drySqu03_Muscle" = "Heart",
#   "drySqu03_Kidney" = "Muscle",
#   "drySqu03_Brain" = "Kidney"
# )
# coldata$Tissue <- ifelse(
#   coldata$Sample %in% names(fix_map),
#   fix_map[coldata$Sample],
#   coldata$Tissue
# )
# 
# ## Recompute PCA data frame and plots
# ntop <- 500
# rv   <- rowVars(vsd_mat)
# select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
# pca_res <- prcomp(t(vsd_mat[select, ]), center = TRUE, scale. = FALSE)
# pca_var <- 100 * (pca_res$sdev^2) / sum(pca_res$sdev^2)
# pca_df  <- as.data.frame(pca_res$x)
# pca_df$Sample <- rownames(pca_df)
# pca_df <- merge(pca_df, coldata, by = "Sample")  # now uses fixed Tissue
# 
# ## Re-draw global PC1xPC2, PC2xPC3, PC1xPC3 as in your script
# # (you can rerun just the three ggplot sections you already have)
# 
# ## Global PC1–PC2
# p12 <- ggplot(pca_df, aes(PC1, PC2, color = Species, shape = Tissue)) +
#   geom_point(size = 3, alpha = 0.8) +
#   xlab(paste0("PC1: ", round(pca_var[1],1), "% variance")) +
#   ylab(paste0("PC2: ", round(pca_var[2],1), "% variance")) +
#   coord_fixed() +
#   theme_bw() +
#   scale_color_manual(values = pal_warmcool, name = "Species")
# p12
# 
# 
# ## Global PC2–PC3
# p23 <- ggplot(pca_df, aes(PC2, PC3, color = Species, shape = Tissue)) +
#   geom_point(size = 3, alpha = 0.8) +
#   xlab(paste0("PC2: ", round(pca_var[2],1), "% variance")) +
#   ylab(paste0("PC3: ", round(pca_var[3],1), "% variance")) +
#   coord_fixed() +
#   theme_bw() +
#   scale_color_manual(values = pal_warmcool, name = "Species")
# p23
# 
# 
# ## Global PC1–PC3
# p13 <- ggplot(pca_df, aes(PC1, PC3, color = Species, shape = Tissue)) +
#   geom_point(size = 3, alpha = 0.8) +
#   xlab(paste0("PC1: ", round(pca_var[1],1), "% variance")) +
#   ylab(paste0("PC3: ", round(pca_var[3],1), "% variance")) +
#   coord_fixed() +
#   theme_bw() +
#   scale_color_manual(values = pal_warmcool, name = "Species")
# p13
# 
# 
# 
# 
# p12$data
# 





# separate up regulated and down regulated


## ========================================================
## 10b) Global GO: convergent dry-up vs dry-down separately
## ========================================================

run_global_go_by_direction <- function(dir_tag = c("UPdry", "DOWNdry"),
                                       OrgDb   = org.Hs.eg.db,
                                       ont     = "BP",
                                       pcut    = 0.1) {
  dir_tag <- match.arg(dir_tag)   # "UPdry" or "DOWNdry"
  tissues <- c("Brain", "Heart", "Kidney", "Liver", "Muscle")
  
  # 1) Collect convergent OGs across tissues for the chosen direction
  all_conv_ogs <- c()
  for (tiss in tissues) {
    if (is.null(conv_ogs[[tiss]])) next
    if (dir_tag == "UPdry") {
      og_dir <- conv_ogs[[tiss]]$dry_up
    } else {
      og_dir <- conv_ogs[[tiss]]$dry_down
    }
    all_conv_ogs <- c(all_conv_ogs, og_dir)
  }
  all_conv_ogs <- unique(na.omit(all_conv_ogs))
  
  if (length(all_conv_ogs) < 5) {
    warning("Too few convergent OGs overall (", length(all_conv_ogs),
            ") to run a global GO analysis for ", dir_tag, ".")
    return(NULL)
  }
  
  message("Global GO (", dir_tag, "): combining ",
          length(all_conv_ogs),
          " convergent OGs across all tissues.")
  
  # 2) Load repGeneName tables and collect gene_name for these OGs
  all_rep <- list()
  for (tiss in tissues) {
    fn_rep <- paste0("OGs_DE_in_multiple_pairs_", tiss, "_repGeneName.csv")
    if (!file.exists(fn_rep)) {
      warning("Rep gene file not found for tissue ", tiss, ": ", fn_rep)
      next
    }
    df_rep <- read.csv(fn_rep, stringsAsFactors = FALSE)
    if (!all(c("OG","gene_name") %in% names(df_rep))) {
      warning("No OG/gene_name columns in ", fn_rep)
      next
    }
    all_rep[[tiss]] <- df_rep
  }
  if (length(all_rep) == 0) {
    warning("No repGeneName files could be read; aborting global GO.")
    return(NULL)
  }
  
  rep_all <- dplyr::bind_rows(all_rep) %>%
    dplyr::filter(OG %in% all_conv_ogs)
  
  gene_list <- unique(na.omit(rep_all$gene_name))
  message("Global GO (", dir_tag, "): ",
          length(gene_list), " unique gene_name entries.")
  
  if (length(gene_list) < 5) {
    warning("Too few gene symbols (", length(gene_list),
            ") after mapping OGs; aborting global GO for ", dir_tag, ".")
    return(NULL)
  }
  
  # 3) SYMBOL -> ENTREZ
  conv <- tryCatch(
    bitr(gene_list,
         fromType = "SYMBOL",
         toType   = "ENTREZID",
         OrgDb    = OrgDb),
    error = function(e) {
      warning("Global GO (", dir_tag, "): bitr conversion failed: ", e$message)
      return(NULL)
    }
  )
  if (is.null(conv) || nrow(conv) == 0) {
    warning("Global GO (", dir_tag, "): no symbols converted to Entrez.")
    return(NULL)
  }
  
  # 4) Run enrichGO
  ego <- tryCatch(
    enrichGO(
      gene          = conv$ENTREZID,
      OrgDb         = OrgDb,
      keyType       = "ENTREZID",
      ont           = ont,
      pAdjustMethod = "BH",
      pvalueCutoff  = pcut,
      qvalueCutoff  = pcut,
      readable      = TRUE
    ),
    error = function(e) {
      warning("Global GO (", dir_tag, "): enrichGO failed: ", e$message)
      return(NULL)
    }
  )
  
  if (is.null(ego) || nrow(as.data.frame(ego)) == 0) {
    warning("Global GO (", dir_tag, "): no enriched GO terms.")
    return(NULL)
  }
  
  ego_df <- as.data.frame(ego)
  tag <- paste0("DRY_", dir_tag)  # DRY_UPdry or DRY_DOWNdry
  
  out_tab <- paste0("GO_BP_Hs_ALL_tissues_", tag, "_global_enrichment.csv")
  write.csv(ego_df, out_tab, row.names = FALSE)
  message("Global GO (", dir_tag, "): wrote table ", out_tab)
  
  # 5) Dotplot
  p_dot <- dotplot(ego, showCategory = 20) +
    ggtitle(paste("GO BP – ALL tissues (human, ", tag, ")", sep = "")) +
    theme_bw() +
    theme(
      plot.title   = element_text(size = 11),
      axis.text.x  = element_text(size = 8),
      axis.text.y  = element_text(size = 8),
      axis.title.x = element_text(size = 9),
      axis.title.y = element_text(size = 9)
    )
  
  print(p_dot)
  
  out_dot <- paste0("GO_BP_Hs_ALL_tissues_", tag, "_global_dotplot.png")
  ggsave(out_dot, p_dot, width = 8, height = 6, dpi = 300)
  message("Global GO (", dir_tag, "): wrote dotplot ", out_dot)
  
  invisible(ego)
}


ego_global_up   <- run_global_go_by_direction("UPdry")
ego_global_down <- run_global_go_by_direction("DOWNdry")


p12
p23






#
# I need to build tissue specific background gene universes 
#


## ---------- Global OG -> representative gene_name mapping (all OGs) ----------

# og_map, gal_map, tae_map, clean_transcript_root(), choose_rep_gene()
# are already defined above in your script.

# Annotate ALL OGs (not just DE in multiple pairs)
og_annot_ref_all <- og_map %>%
  mutate(transcript_root = clean_transcript_root(gene_id)) %>%
  left_join(gal_map, by = c("transcript_root" = "transcriptID")) %>%
  left_join(tae_map, by = c("transcript_root" = "transcriptID")) %>%
  arrange(OG, species)

# Choose one representative gene_name per OG (no Tissue this time)
rep_names_all <- og_annot_ref_all %>%
  dplyr::group_by(OG) %>%
  dplyr::do(choose_rep_gene(.)) %>%
  dplyr::ungroup()

# This table will be used to map ANY OG to a gene_name
og_rep_all <- rep_names_all %>%
  dplyr::select(OG, gene_name, name_source)

# Optional: save for inspection
write.csv(og_rep_all,
          "OGs_ALL_repGeneName.csv",
          row.names = FALSE)


## Helper: GO for one tissue / OG subset (UPdry or DOWNdry, human only)
run_go_for_tissue_dry <- function(tiss,
                                  og_subset,
                                  OrgDb = org.Hs.eg.db,
                                  ont = "BP",
                                  pcut = 0.05,
                                  dir_tag = c("UPdry","DOWNdry")) {
  dir_tag <- match.arg(dir_tag)
  
  ## ---------- 1) Foreground: convergent OGs in this tissue ----------
  fn_rep <- paste0("OGs_DE_in_multiple_pairs_", tiss, "_repGeneName.csv")
  if (!file.exists(fn_rep)) {
    warning("Rep gene file not found for tissue ", tiss, ": ", fn_rep)
    return(NULL)
  }
  df_rep <- read.csv(fn_rep, stringsAsFactors = FALSE)
  if (!all(c("OG","gene_name") %in% names(df_rep))) {
    warning("No OG/gene_name columns in ", fn_rep)
    return(NULL)
  }
  
  if (!is.null(og_subset)) {
    df_rep <- df_rep %>% dplyr::filter(OG %in% og_subset)
  }
  
  fg_symbols <- unique(na.omit(df_rep$gene_name))
  message("Tissue ", tiss, " – GO (", dir_tag, ", human) with ",
          length(fg_symbols), " foreground gene_name entries.")
  
  if (length(fg_symbols) < 5) {
    warning("Too few foreground genes (", length(fg_symbols),
            ") in tissue ", tiss, " for GO (", dir_tag, "); skipping.")
    return(NULL)
  }
  
  ## ---------- 2) Background: all OGs expressed in this tissue ----------
  # dds: DESeqDataSet with OG rows
  if (!exists("dds")) {
    stop("dds object not found in environment; needed for tissue-specific background.")
  }
  
  # Columns (samples) belonging to this tissue
  cols_tiss <- which(dds$Tissue == tiss)
  if (length(cols_tiss) == 0) {
    warning("No samples for tissue ", tiss, " in dds; cannot build background.")
    return(NULL)
  }
  
  # OGs with any non-zero counts in this tissue
  cts_tiss <- counts(dds)[, cols_tiss, drop = FALSE]
  bg_ogs   <- rownames(cts_tiss)[rowSums(cts_tiss > 0) > 0]
  
  if (length(bg_ogs) < 5) {
    warning("Too few background OGs (", length(bg_ogs),
            ") in tissue ", tiss, "; skipping GO.")
    return(NULL)
  }
  
  # Map those OGs to gene_name via og_rep_all
  if (!exists("og_rep_all")) {
    stop("og_rep_all (global OG->gene_name table) not found; build it before GO.")
  }
  
  bg_map <- og_rep_all %>%
    dplyr::filter(OG %in% bg_ogs, !is.na(gene_name), gene_name != "")
  
  bg_symbols <- unique(bg_map$gene_name)
  message("Tissue ", tiss, " – background has ",
          length(bg_symbols), " gene_name entries (expressed in ", tiss, ").")
  
  if (length(bg_symbols) < 5) {
    warning("Too few background symbols (", length(bg_symbols),
            ") in tissue ", tiss, "; skipping GO.")
    return(NULL)
  }
  
  ## ---------- 3) Map both foreground and background to ENTREZ ----------
  fg_conv <- tryCatch(
    bitr(fg_symbols,
         fromType = "SYMBOL",
         toType   = "ENTREZID",
         OrgDb    = OrgDb),
    error = function(e) {
      warning("bitr conversion failed for FOREGROUND in ", tiss, " (", dir_tag, "): ", e$message)
      return(NULL)
    }
  )
  if (is.null(fg_conv) || nrow(fg_conv) == 0) {
    warning("No foreground symbols converted to Entrez for human in ", tiss, " (", dir_tag, ")")
    return(NULL)
  }
  
  bg_conv <- tryCatch(
    bitr(bg_symbols,
         fromType = "SYMBOL",
         toType   = "ENTREZID",
         OrgDb    = OrgDb),
    error = function(e) {
      warning("bitr conversion failed for BACKGROUND in ", tiss, " (", dir_tag, "): ", e$message)
      return(NULL)
    }
  )
  if (is.null(bg_conv) || nrow(bg_conv) == 0) {
    warning("No background symbols converted to Entrez for human in ", tiss, " (", dir_tag, ")")
    return(NULL)
  }
  
  ## ---------- 4) Run enrichGO with custom universe ----------
  ego <- tryCatch(
    enrichGO(
      gene          = fg_conv$ENTREZID,
      universe      = bg_conv$ENTREZID,   # <-- key change
      OrgDb         = OrgDb,
      keyType       = "ENTREZID",
      ont           = ont,
      pAdjustMethod = "BH",
      pvalueCutoff  = pcut,
      qvalueCutoff  = pcut,
      readable      = TRUE
    ),
    error = function(e) {
      warning("enrichGO failed for human in ", tiss, " (", dir_tag, "): ", e$message)
      return(NULL)
    }
  )
  
  if (is.null(ego) || nrow(as.data.frame(ego)) == 0) {
    warning("No enriched GO terms for human in ", tiss, " (", dir_tag, ")")
    return(NULL)
  }
  
  ego_df <- as.data.frame(ego)
  tag <- paste0("DRY_", dir_tag)  # e.g. DRY_UPdry
  
  out_tab <- paste0("GO_BP_Hs_", tiss, "_", tag, "_enrichment.csv")
  write.csv(ego_df, out_tab, row.names = FALSE)
  message("  Wrote GO table: ", out_tab)
  
  p_dot <- dotplot(ego, showCategory = 20) +
    ggtitle(paste("GO BP –", tiss, "(human,", tag, ")")) +
    theme_bw() +
    theme(
      plot.title   = element_text(size = 11),
      axis.text.x  = element_text(size = 8),
      axis.text.y  = element_text(size = 8),
      axis.title.x = element_text(size = 9),
      axis.title.y = element_text(size = 9)
    )
  out_dot <- paste0("GO_BP_Hs_", tiss, "_", tag, "_dotplot.png")
  ggsave(out_dot, p_dot, width = 8, height = 6, dpi = 300)
  message("  Wrote dotplot: ", out_dot)
  
  ego
}


## Run GO for convergent dry-up / dry-down in each tissue (human only)
for (tiss in c("Brain", "Heart", "Kidney", "Liver", "Muscle")) {
  og_up   <- conv_ogs[[tiss]]$dry_up
  og_down <- conv_ogs[[tiss]]$dry_down
  
  if (!is.null(og_up) && length(og_up) >= 5) {
    message("\n=== ", tiss, " – convergently UP in dry (human GO) ===")
    run_go_for_tissue_dry(tiss, og_up, dir_tag = "UPdry")
  }
  
  if (!is.null(og_down) && length(og_down) >= 5) {
    message("\n=== ", tiss, " – convergently DOWN in dry (human GO) ===")
    run_go_for_tissue_dry(tiss, og_down, dir_tag = "DOWNdry")
  }
}


view_go_tissue <- function(tiss,
                           dir_tag = c("UPdry","DOWNdry"),
                           ont = "BP",
                           pcut = 0.05) {
  dir_tag <- match.arg(dir_tag)
  
  og_subset <- if (dir_tag == "UPdry") {
    conv_ogs[[tiss]]$dry_up
  } else {
    conv_ogs[[tiss]]$dry_down
  }
  
  if (is.null(og_subset) || length(og_subset) < 5) {
    warning("Not enough convergent OGs for ", tiss, " (", dir_tag, ")")
    return(invisible(NULL))
  }
  
  message("Re-running GO for ", tiss, " (", dir_tag, ", human only)")
  ego <- run_go_for_tissue_dry(
    tiss      = tiss,
    og_subset = og_subset,
    OrgDb     = org.Hs.eg.db,
    ont       = ont,
    pcut      = pcut,
    dir_tag   = dir_tag
  )
  # run_go_for_tissue_dry already prints a dotplot and saves PNG;
  # this function just returns the enrichResult object for inspection.
  invisible(ego)
}

# Examples:
view_go_tissue("Kidney", "UPdry")
view_go_tissue("Kidney", "DOWNdry")
view_go_tissue("Muscle", "UPdry")
view_go_tissue("Muscle", "DOWNdry")
view_go_tissue("Brain", "UPdry")
view_go_tissue("Brain", "DOWNdry")
view_go_tissue("Liver", "UPdry")
view_go_tissue("Liver", "DOWNdry")
view_go_tissue("Heart", "UPdry")
view_go_tissue("Heart", "DOWNdry")









# make a heat map to see convergence in GO terms 




############################################################
## Per‑pair, per‑tissue GO with tissue‑specific background
## and heatmaps of GO term overlap
############################################################

library(dplyr)
library(tidyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(pheatmap)

## 0) Global OG -> representative gene_name mapping (for ALL OGs)
##    (independent of DE / tissue)

og_annot_ref_all <- og_map %>%
  mutate(transcript_root = clean_transcript_root(gene_id)) %>%
  left_join(gal_map, by = c("transcript_root" = "transcriptID")) %>%
  left_join(tae_map, by = c("transcript_root" = "transcriptID")) %>%
  arrange(OG, species)

rep_names_all <- og_annot_ref_all %>%
  dplyr::group_by(OG) %>%
  dplyr::do(choose_rep_gene(.)) %>%
  dplyr::ungroup()

og_rep_all <- rep_names_all %>%
  dplyr::select(OG, gene_name, name_source)

## Optional: inspect
# write.csv(og_rep_all, "OGs_ALL_repGeneName.csv", row.names = FALSE)


############################################################
## 1) Tissue-specific background (universe)
############################################################

get_tissue_background_entrez <- function(tiss,
                                         dds,
                                         og_rep_all,
                                         OrgDb) {
  # columns (samples) belonging to this tissue
  cols_tiss <- which(dds$Tissue == tiss)
  if (length(cols_tiss) == 0) {
    warning("No samples for tissue ", tiss, " in dds; cannot build background.")
    return(NULL)
  }
  
  # OGs with any non-zero counts in this tissue
  cts_tiss <- counts(dds)[, cols_tiss, drop = FALSE]
  bg_ogs   <- rownames(cts_tiss)[rowSums(cts_tiss > 0) > 0]
  if (length(bg_ogs) < 5) {
    warning("Too few background OGs (", length(bg_ogs),
            ") in tissue ", tiss, "; skipping background.")
    return(NULL)
  }
  
  # Map those OGs to gene_name
  bg_map <- og_rep_all %>%
    dplyr::filter(OG %in% bg_ogs, !is.na(gene_name), gene_name != "")
  bg_symbols <- unique(bg_map$gene_name)
  
  if (length(bg_symbols) < 5) {
    warning("Too few background symbols (", length(bg_symbols),
            ") in tissue ", tiss, "; skipping background.")
    return(NULL)
  }
  
  # SYMBOL -> ENTREZ
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
## 2) Dry-relative DE OGs for one tissue & pair
############################################################

get_dry_direction_for_pair <- function(tiss,
                                       pair_row,   # one row from pairs_arid
                                       de_sig_tiss # de_sig[[tiss]]
) {
  pid <- pair_row$pair_id
  spA <- pair_row$spA
  spB <- pair_row$spB
  dry <- pair_row$dry_species
  
  if (!pid %in% names(de_sig_tiss)) return(NULL)
  
  sig_ogs <- de_sig_tiss[[pid]]
  if (is.null(sig_ogs) || length(sig_ogs) == 0) return(NULL)
  
  fname <- paste0("DE_OG_selected_", spA, "_vs_", spB, "_", tiss, "_kallisto.csv")
  if (!file.exists(fname)) return(NULL)
  
  de_tab <- read.csv(fname, row.names = 1, check.names = FALSE)
  de_tab$OG <- rownames(de_tab)
  
  de_sub <- de_tab[de_tab$OG %in% sig_ogs, c("OG", "log2FoldChange")]
  if (nrow(de_sub) == 0) return(NULL)
  
  # direction relative to dry species
  dir_rel <- if (dry == spA) {
    ifelse(de_sub$log2FoldChange > 0, "dry_up",
           ifelse(de_sub$log2FoldChange < 0, "dry_down", "zero"))
  } else if (dry == spB) {
    ifelse(de_sub$log2FoldChange > 0, "dry_down",
           ifelse(de_sub$log2FoldChange < 0, "dry_up", "zero"))
  } else {
    stop("dry_species not equal to spA or spB for pair ", pid)
  }
  
  data.frame(
    Tissue       = tiss,
    pair_id      = pid,
    OG           = de_sub$OG,
    directionDry = dir_rel,
    stringsAsFactors = FALSE
  )
}


############################################################
## 3) Run GO per tissue × pair × direction
############################################################

run_go_pairwise_all <- function(
    tissues    = c("Brain", "Heart", "Kidney", "Liver", "Muscle"),
    pcut       = 0.05,
    OrgDb      = org.Hs.eg.db,
    de_sig_in  = de_sig,
    dds_in     = dds,
    og_rep_all_in = og_rep_all,
    pairs_arid_in = pairs_arid
) {
  results_go <- list()
  
  for (tiss in tissues) {
    message("\n=== Tissue: ", tiss, " ===")
    
    # 3A) tissue-specific background (universe)
    bg_universe <- get_tissue_background_entrez(
      tiss       = tiss,
      dds        = dds_in,
      og_rep_all = og_rep_all_in,
      OrgDb      = OrgDb
    )
    if (is.null(bg_universe)) {
      warning("Skipping tissue ", tiss, " (no background).")
      next
    }
    
    results_go[[tiss]] <- list(UPdry = list(), DOWNdry = list())
    
    de_sig_tiss_all <- de_sig_in[[tiss]]
    if (is.null(de_sig_tiss_all)) {
      message("  No DE results stored for tissue ", tiss, "; skipping.")
      next
    }
    
    # one row per pair_id in pairs_arid
    for (i in seq_len(nrow(pairs_arid_in))) {
      pair_row <- pairs_arid_in[i, ]
      pid      <- pair_row$pair_id
      
      message("  Pair: ", pid)
      if (!(pid %in% names(de_sig_tiss_all))) {
        message("    No DE sig OGs; skipping.")
        next
      }
      
      # 3B) get dry-relative direction for this pair
      dir_df <- get_dry_direction_for_pair(tiss, pair_row, de_sig_tiss_all)
      if (is.null(dir_df) || nrow(dir_df) == 0) {
        message("    No OGs with dry direction; skipping.")
        next
      }
      
      for (dir_tag in c("UPdry", "DOWNdry")) {
        og_dir <- dir_df$OG[dir_df$directionDry ==
                              ifelse(dir_tag == "UPdry", "dry_up", "dry_down")]
        og_dir <- unique(og_dir)
        message("    Direction ", dir_tag, ": ", length(og_dir), " OGs.")
        if (length(og_dir) < 5) {
          message("      Too few OGs for GO; skipping.")
          next
        }
        
        ## Map OGs to gene_name
        fg_map <- og_rep_all_in %>%
          dplyr::filter(OG %in% og_dir, !is.na(gene_name), gene_name != "")
        fg_symbols <- unique(fg_map$gene_name)
        if (length(fg_symbols) < 5) {
          message("      Too few mapped gene symbols; skipping.")
          next
        }
        
        # SYMBOL -> ENTREZ
        fg_conv <- tryCatch(
          bitr(fg_symbols,
               fromType = "SYMBOL",
               toType   = "ENTREZID",
               OrgDb    = OrgDb),
          error = function(e) {
            warning("Foreground bitr failed for ", tiss, " ", pid, " ", dir_tag,
                    ": ", e$message)
            return(NULL)
          }
        )
        if (is.null(fg_conv) || nrow(fg_conv) == 0) {
          message("      No foreground symbols converted; skipping.")
          next
        }
        
        # enrichGO with tissue-specific universe
        ego <- tryCatch(
          enrichGO(
            gene          = fg_conv$ENTREZID,
            universe      = bg_universe,
            OrgDb         = OrgDb,
            keyType       = "ENTREZID",
            ont           = "BP",
            pAdjustMethod = "BH",
            pvalueCutoff  = 1,    # keep all, threshold later
            qvalueCutoff  = 1,
            readable      = TRUE
          ),
          error = function(e) {
            warning("enrichGO failed for ", tiss, " ", pid, " ", dir_tag,
                    ": ", e$message)
            return(NULL)
          }
        )
        if (is.null(ego) || nrow(as.data.frame(ego)) == 0) {
          message("      No GO terms returned.")
          next
        }
        
        results_go[[tiss]][[dir_tag]][[pid]] <- ego
        
        # Optional: write out tables
        out_tab <- paste0("GO_BP_Hs_", tiss, "_", pid, "_", dir_tag,
                          "_pair_enrichment.csv")
        write.csv(as.data.frame(ego), out_tab, row.names = FALSE)
      }
    }
  }
  
  results_go
}

## Now run it:
go_pairwise_results <- run_go_pairwise_all()


############################################################
## 4) Build GO overlap heatmap for one tissue & direction
############################################################

## For a given tissue and direction (UPdry / DOWNdry), this:
## - collects GO terms across all pairs in that tissue,
## - builds a matrix of –log10(p.adjust) per GO term × pair,
##   with 0 = term present but padj > alpha (non-significant),
##   NA = term not present in that pair’s result.

build_go_overlap_matrix <- function(tiss,
                                    dir_tag = c("UPdry", "DOWNdry"),
                                    go_results = go_pairwise_results,
                                    alpha = 0.05) {
  dir_tag <- match.arg(dir_tag)
  
  if (is.null(go_results[[tiss]]) ||
      is.null(go_results[[tiss]][[dir_tag]])) {
    stop("No GO results for tissue ", tiss, " and direction ", dir_tag)
  }
  
  res_list <- go_results[[tiss]][[dir_tag]]
  if (length(res_list) == 0) {
    stop("No GO results for ", tiss, " ", dir_tag)
  }
  
  # collect all GO IDs across pairs
  all_terms <- unique(unlist(
    lapply(res_list, function(ego) as.data.frame(ego)$ID)
  ))
  
  if (length(all_terms) == 0) {
    stop("No GO terms found for ", tiss, " ", dir_tag)
  }
  
  pairs <- names(res_list)
  mat   <- matrix(NA_real_,
                  nrow = length(all_terms),
                  ncol = length(pairs),
                  dimnames = list(all_terms, pairs))
  
  # fill matrix
  for (pid in pairs) {
    ego <- res_list[[pid]]
    if (is.null(ego)) next
    df <- as.data.frame(ego)
    rownames(df) <- df$ID
    
    # for each term, get padj (if present)
    common_ids <- intersect(all_terms, rownames(df))
    if (length(common_ids) == 0) next
    
    padj <- df[common_ids, "p.adjust"]
    # encode:
    # padj <= alpha -> -log10(padj) (significant)
    # padj > alpha  -> 0 (non-significant but present)
    # NA -> NA
    vals <- rep(NA_real_, length(common_ids))
    names(vals) <- common_ids
    sig_idx     <- which(!is.na(padj) & padj <= alpha)
    nonsig_idx  <- which(!is.na(padj) & padj > alpha)
    if (length(sig_idx) > 0) {
      vals[sig_idx] <- -log10(padj[sig_idx])
    }
    if (length(nonsig_idx) > 0) {
      vals[nonsig_idx] <- 0  # present but non-significant
    }
    
    mat[names(vals), pid] <- vals
  }
  
  mat
}

## Plot heatmap: 
##  - white: term not present for that pair (NA)
##  - grey: term present but not significant (0)
##  - colored: significant, intensity ~ -log10(padj)




plot_go_overlap_heatmap <- function(tiss,
                                    dir_tag = c("UPdry", "DOWNdry"),
                                    go_results = go_pairwise_results,
                                    alpha = 0.05,
                                    top_n_terms = 50) {
  dir_tag <- match.arg(dir_tag)
  mat <- build_go_overlap_matrix(tiss, dir_tag, go_results, alpha)
  
  # optionally select top N terms by how often they’re significant
  sig_counts <- rowSums(mat > 0, na.rm = TRUE)
  keep_ids   <- names(sort(sig_counts, decreasing = TRUE))[1:min(top_n_terms,
                                                                 length(sig_counts))]
  mat_filt   <- mat[keep_ids, , drop = FALSE]
  
  ## ---------------------------------------------
  ## Map GO IDs (rownames) -> Descriptions
  ## ---------------------------------------------
  res_list <- go_results[[tiss]][[dir_tag]]
  # collect ID/Description from all pairs
  id_desc_df <- do.call(
    rbind,
    lapply(res_list, function(ego) {
      if (is.null(ego)) return(NULL)
      df <- as.data.frame(ego)
      df[, c("ID", "Description")]
    })
  )
  id_desc_df <- unique(id_desc_df)
  
  # vector: ID -> Description
  desc_map <- id_desc_df$Description
  names(desc_map) <- id_desc_df$ID
  
  # for each row in mat_filt, pick its Description (or keep ID if missing)
  ids <- rownames(mat_filt)
  row_labels <- ifelse(ids %in% names(desc_map), desc_map[ids], ids)
  
  # make descriptions unique for plotting (pheatmap requires unique rownames)
  # append ID where duplicates occur
  dup <- duplicated(row_labels)
  if (any(dup)) {
    row_labels[dup] <- paste0(row_labels[dup], " (", ids[dup], ")")
  }
  rownames(mat_filt) <- row_labels
  
  ## ---------------------------------------------
  ## Plot: handle "no significant" vs "some significant"
  ## ---------------------------------------------
  
  # check if there are any significant values (>0)
  has_pos <- any(mat_filt > 0, na.rm = TRUE)
  
  if (!has_pos) {
    message("No significant GO terms for ", tiss, " ", dir_tag,
            " (all cells are 0 or NA). Plotting presence/absence only.")
    
    # simple breaks: NA = white, 0 = grey
    breaks <- c(-0.5, 0.5)   # everything 0 will fall into the single bin
    colors <- c("grey80")
    
    pheatmap(mat_filt,
             color          = colors,
             breaks         = breaks,
             na_col         = "white",
             cluster_rows   = TRUE,
             cluster_cols   = TRUE,
             main           = paste("GO BP overlap –", tiss, dir_tag,
                                    "\nno significant terms (0 = present, >0 none)"),
             fontsize_row   = 6,
             fontsize_col   = 8,
             border_color   = NA)
    
  } else {
    # continuous palette: NA=white, 0=grey, >0 = yellow->red
    max_val <- max(mat_filt, na.rm = TRUE)
    if (!is.finite(max_val) || max_val <= 0) max_val <- 1
    
    # ensure strictly increasing, unique breaks
    seq_part <- seq(0.00001, max_val, length.out = 100)
    breaks <- unique(c(0, 0.00001, seq_part))
    
    colors <- c("grey80", colorRampPalette(c("yellow", "red"))(length(breaks) - 2))
    
    pheatmap(mat_filt,
             color          = colors,
             breaks         = breaks,
             na_col         = "white",
             cluster_rows   = TRUE,
             cluster_cols   = TRUE,
             main           = paste("GO BP overlap –", tiss, dir_tag,
                                    "\n0 = overlap but padj > ", alpha,
                                    "; color = -log10(padj)"),
             fontsize_row   = 6,
             fontsize_col   = 8,
             border_color   = NA)
  }
}


## EXAMPLES:
## For each tissue, you can do (UP in dry):
plot_go_overlap_heatmap("Kidney", "UPdry")
plot_go_overlap_heatmap("Kidney", "DOWNdry")
plot_go_overlap_heatmap("Liver",  "UPdry")
plot_go_overlap_heatmap("Liver",  "DOWNdry")
plot_go_overlap_heatmap("Brain",  "UPdry")
plot_go_overlap_heatmap("Brain",  "DOWNdry")
plot_go_overlap_heatmap("Heart",  "UPdry")
plot_go_overlap_heatmap("Heart",  "DOWNdry")
plot_go_overlap_heatmap("Muscle",  "UPdry")
plot_go_overlap_heatmap("Muscle",  "DOWNdry")

                    












############################################################
## GO per tissue × pair (all DE OGs, no up/down split)
############################################################

run_go_pairwise_all_nodir <- function(
    tissues       = c("Brain", "Heart", "Kidney", "Liver", "Muscle"),
    OrgDb         = org.Hs.eg.db,
    de_sig_in     = de_sig,
    dds_in        = dds,
    og_rep_all_in = og_rep_all,
    pairs_in      = pairs_arid
) {
  results_go <- list()
  
  for (tiss in tissues) {
    message("\n=== Tissue: ", tiss, " ===")
    
    # 1) tissue-specific background (universe)
    bg_universe <- get_tissue_background_entrez(
      tiss       = tiss,
      dds        = dds_in,
      og_rep_all = og_rep_all_in,
      OrgDb      = OrgDb
    )
    if (is.null(bg_universe)) {
      warning("Skipping tissue ", tiss, " (no background).")
      next
    }
    
    results_go[[tiss]] <- list()
    
    de_sig_tiss <- de_sig_in[[tiss]]
    if (is.null(de_sig_tiss)) {
      message("  No DE results stored for tissue ", tiss, "; skipping.")
      next
    }
    
    # 2) loop over pairs
    for (i in seq_len(nrow(pairs_in))) {
      pair_row <- pairs_in[i, ]
      pid      <- pair_row$pair_id
      
      message("  Pair: ", pid)
      if (!(pid %in% names(de_sig_tiss))) {
        message("    No DE sig OGs; skipping.")
        next
      }
      
      og_dir <- unique(de_sig_tiss[[pid]])  # all DE OGs, any direction
      message("    DE OGs (both directions): ", length(og_dir))
      if (length(og_dir) < 5) {
        message("      Too few OGs for GO; skipping.")
        next
      }
      
      ## Map OGs to gene_name
      fg_map <- og_rep_all_in %>%
        dplyr::filter(OG %in% og_dir,
                      !is.na(gene_name),
                      gene_name != "")
      fg_symbols <- unique(fg_map$gene_name)
      if (length(fg_symbols) < 5) {
        message("      Too few mapped gene symbols; skipping.")
        next
      }
      
      # SYMBOL -> ENTREZ
      fg_conv <- tryCatch(
        bitr(fg_symbols,
             fromType = "SYMBOL",
             toType   = "ENTREZID",
             OrgDb    = OrgDb),
        error = function(e) {
          warning("Foreground bitr failed for ", tiss, " ", pid, ": ", e$message)
          return(NULL)
        }
      )
      if (is.null(fg_conv) || nrow(fg_conv) == 0) {
        message("      No foreground symbols converted; skipping.")
        next
      }
      
      # enrichGO with tissue-specific universe
      ego <- tryCatch(
        enrichGO(
          gene          = fg_conv$ENTREZID,
          universe      = bg_universe,
          OrgDb         = OrgDb,
          keyType       = "ENTREZID",
          ont           = "BP",
          pAdjustMethod = "BH",
          pvalueCutoff  = 1,    # keep all; threshold later
          qvalueCutoff  = 1,
          readable      = TRUE
        ),
        error = function(e) {
          warning("enrichGO failed for ", tiss, " ", pid, ": ", e$message)
          return(NULL)
        }
      )
      if (is.null(ego) || nrow(as.data.frame(ego)) == 0) {
        message("      No GO terms returned.")
        next
      }
      
      results_go[[tiss]][[pid]] <- ego
      
      # Optional: write out tables
      out_tab <- paste0("GO_BP_Hs_", tiss, "_", pid, "_ALLdir_pair_enrichment.csv")
      write.csv(as.data.frame(ego), out_tab, row.names = FALSE)
    }
  }
  
  results_go
}

## Run this once:
go_pairwise_results_all <- run_go_pairwise_all_nodir()



############################################################
## GO overlap matrix & heatmap (all DE, no up/down split)
############################################################

build_go_overlap_matrix_all <- function(tiss,
                                        go_results = go_pairwise_results_all,
                                        alpha = 0.05) {
  if (is.null(go_results[[tiss]])) {
    stop("No GO results for tissue ", tiss)
  }
  
  res_list <- go_results[[tiss]]
  if (length(res_list) == 0) {
    stop("No GO results for tissue ", tiss)
  }
  
  # collect all GO IDs across pairs
  all_terms <- unique(unlist(
    lapply(res_list, function(ego) as.data.frame(ego)$ID)
  ))
  
  if (length(all_terms) == 0) {
    stop("No GO terms found for tissue ", tiss)
  }
  
  pairs <- names(res_list)
  mat   <- matrix(NA_real_,
                  nrow = length(all_terms),
                  ncol = length(pairs),
                  dimnames = list(all_terms, pairs))
  
  # fill matrix
  for (pid in pairs) {
    ego <- res_list[[pid]]
    if (is.null(ego)) next
    df <- as.data.frame(ego)
    rownames(df) <- df$ID
    
    common_ids <- intersect(all_terms, rownames(df))
    if (length(common_ids) == 0) next
    
    padj <- df[common_ids, "p.adjust"]
    vals <- rep(NA_real_, length(common_ids))
    names(vals) <- common_ids
    
    sig_idx    <- which(!is.na(padj) & padj <= alpha)
    nonsig_idx <- which(!is.na(padj) & padj > alpha)
    
    if (length(sig_idx) > 0) {
      vals[sig_idx] <- -log10(padj[sig_idx])     # significant
    }
    if (length(nonsig_idx) > 0) {
      vals[nonsig_idx] <- 0                      # present but non-significant
    }
    
    mat[names(vals), pid] <- vals
  }
  
  mat
}


plot_go_overlap_heatmap_all <- function(tiss,
                                        go_results = go_pairwise_results_all,
                                        alpha = 0.05,
                                        top_n_terms = 50) {
  mat <- build_go_overlap_matrix_all(tiss, go_results, alpha)
  
  # optionally select top N terms by how often they’re significant
  sig_counts <- rowSums(mat > 0, na.rm = TRUE)
  keep_ids   <- names(sort(sig_counts, decreasing = TRUE))[1:min(top_n_terms,
                                                                 length(sig_counts))]
  mat_filt   <- mat[keep_ids, , drop = FALSE]
  
  ## Map GO IDs -> Descriptions
  res_list <- go_results[[tiss]]
  id_desc_df <- do.call(
    rbind,
    lapply(res_list, function(ego) {
      if (is.null(ego)) return(NULL)
      df <- as.data.frame(ego)
      df[, c("ID", "Description")]
    })
  )
  id_desc_df <- unique(id_desc_df)
  
  desc_map <- id_desc_df$Description
  names(desc_map) <- id_desc_df$ID
  
  ids <- rownames(mat_filt)
  row_labels <- ifelse(ids %in% names(desc_map), desc_map[ids], ids)
  
  dup <- duplicated(row_labels)
  if (any(dup)) {
    row_labels[dup] <- paste0(row_labels[dup], " (", ids[dup], ")")
  }
  rownames(mat_filt) <- row_labels
  
  # check if any significant values
  has_pos <- any(mat_filt > 0, na.rm = TRUE)
  
  if (!has_pos) {
    message("No significant GO terms for ", tiss,
            " (all cells are 0 or NA). Plotting presence/absence only.")
    
    breaks <- c(-0.5, 0.5)
    colors <- c("grey80")
    
    pheatmap(mat_filt,
             color          = colors,
             breaks         = breaks,
             na_col         = "white",
             cluster_rows   = TRUE,
             cluster_cols   = TRUE,
             main           = paste("GO BP overlap –", tiss,
                                    "\nno significant terms (0 = present, >0 none)"),
             fontsize_row   = 6,
             fontsize_col   = 8,
             border_color   = NA)
    
  } else {
    max_val <- max(mat_filt, na.rm = TRUE)
    if (!is.finite(max_val) || max_val <= 0) max_val <- 1
    
    seq_part <- seq(0.00001, max_val, length.out = 100)
    breaks <- unique(c(0, 0.00001, seq_part))
    
    colors <- c("grey80", colorRampPalette(c("yellow", "red"))(length(breaks) - 2))
    
    pheatmap(mat_filt,
             color          = colors,
             breaks         = breaks,
             na_col         = "white",
             cluster_rows   = TRUE,
             cluster_cols   = TRUE,
             main           = paste("GO BP overlap –", tiss,
                                    "\n0 = overlap but padj > ", alpha,
                                    "; color = -log10(padj)"),
             fontsize_row   = 6,
             fontsize_col   = 8,
             border_color   = NA)
  }
}



# Kidney
plot_go_overlap_heatmap_all("Kidney")
# Liver
plot_go_overlap_heatmap_all("Liver")
# Brain
plot_go_overlap_heatmap_all("Brain")
# Heart
plot_go_overlap_heatmap_all("Heart")
# Muscle
plot_go_overlap_heatmap_all("Muscle")










############################################################
## Plot: # DE genes vs divergence time, colored by tissue
## Requires:
##  - subosc.nex at:
##      /n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/nf_pipeline/subosc.nex
##  - DE CSVs:
##      DE_OG_selected_<spA>_vs_<spB>_<Tissue>_kallisto.csv
############################################################

############################################################
## DE genes vs divergence time (patristic), by tissue
## - one point per pair × tissue
## - y: # DE genes (padj < 0.05, |log2FC| >= 1)
## - x: divergence (cophenetic distance)
## - color: tissue
## - label: species pair (pair_id) with ggrepel
############################################################

setwd("/n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/nf_pipeline/thaDol_mapping_test")

library(ape)
library(phytools)
library(dplyr)
library(ggplot2)
library(readr)
library(ggrepel)

##----------------------------------------------------------
## 1) Read tree and compute divergence for your 5 pairs
##----------------------------------------------------------

tree_path <- "/n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/nf_pipeline/subosc.nex"
tree <- read.nexus(tree_path)

# Map short species codes to tip labels in the tree
spcode_to_tip <- c(
  dysSti = "Dysith_stirax_CPRB145",
  sakCri = "Sakeides_critus_CPJB132",
  sakCan = "Sakerus_canad_K5831",
  sakLuc = "Sakerus_luctu_MZU92813",
  thaAmb = "Thamphil_ambus_MZU86040",
  thaPel = "Thamphil_pelzel_MZU94631",
  thaAtr = "Thamphil_atrcha_L46549",
  thaBer = "Thamphil_berna_LSU45997",
  thaRuf = "Thamphil_rufcap_LSU58333",
  thaTor = "Thamphil_tortus_LSU13900"
)

# Species pairs (must match DE filenames DE_OG_selected_<spA>_vs_<spB>_...)
pairs_of_interest <- tribble(
  ~pair_id,           ~spA,     ~spB,
  "dysSti_vs_sakCri", "dysSti", "sakCri",
  "sakLuc_vs_sakCan", "sakLuc", "sakCan",
  "thaAmb_vs_thaPel", "thaAmb", "thaPel",
  "thaBer_vs_thaAtr", "thaBer", "thaAtr",
  "thaRuf_vs_thaTor", "thaRuf", "thaTor"
)

# subset tree to only required tips
needed_tips <- unique(spcode_to_tip[c(pairs_of_interest$spA,
                                      pairs_of_interest$spB)])
tree_sub <- keep.tip(tree, needed_tips)

# pairwise patristic distances
pat <- cophenetic.phylo(tree_sub)

# divergence table
divergence_df <- pairs_of_interest %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    tipA      = spcode_to_tip[spA],
    tipB      = spcode_to_tip[spB],
    divergence = pat[tipA, tipB]
  ) %>%
  dplyr::ungroup() %>%
  dplyr::select(pair_id, spA, spB, divergence)

print(divergence_df)

##----------------------------------------------------------
## 2) Count DE genes per pair × tissue
##----------------------------------------------------------

tissues_to_run <- c("Brain", "Heart", "Liver", "Muscle", "Kidney")

# helper: count DE genes for one pair/tissue
# uses same thresholds as your DE script: padj < 0.05, |log2FC| >= 1
count_de_for_pair_tissue <- function(tiss, spA, spB,
                                     prefix = "DE_OG_selected",
                                     suffix = "_kallisto.csv") {
  fname <- paste0(prefix, "_", spA, "_vs_", spB, "_", tiss, suffix)
  if (!file.exists(fname)) {
    warning("DE file not found: ", fname)
    return(NA_integer_)
  }
  tab <- suppressMessages(read_csv(fname, show_col_types = FALSE))
  if (!all(c("padj", "log2FoldChange") %in% colnames(tab))) {
    warning("Missing padj/log2FoldChange in ", fname)
    return(NA_integer_)
  }
  
  tab_filt <- tab %>%
    dplyr::filter(!is.na(padj),
                  padj < 0.05,
                  abs(log2FoldChange) >= 1)
  
  nrow(tab_filt)
}

de_counts_list <- list()

for (tiss in tissues_to_run) {
  for (i in seq_len(nrow(pairs_of_interest))) {
    pid <- pairs_of_interest$pair_id[i]
    spA <- pairs_of_interest$spA[i]
    spB <- pairs_of_interest$spB[i]
    
    n_de <- count_de_for_pair_tissue(tiss, spA, spB)
    
    de_counts_list[[length(de_counts_list) + 1]] <- data.frame(
      Tissue   = tiss,
      pair_id  = pid,
      spA      = spA,
      spB      = spB,
      n_DE     = n_de,
      stringsAsFactors = FALSE
    )
  }
}

de_counts_df <- bind_rows(de_counts_list) %>%
  dplyr::filter(!is.na(n_DE))

print(de_counts_df)

##----------------------------------------------------------
## 3) Merge DE counts with divergence
##----------------------------------------------------------

plot_df <- de_counts_df %>%
  dplyr::left_join(divergence_df, by = c("pair_id", "spA", "spB"))

print(plot_df)

##----------------------------------------------------------
## 4) Plot: #DE genes vs divergence, colored by tissue,
##          labeled by pair_id with ggrepel
##----------------------------------------------------------

p <- ggplot(plot_df, aes(x = divergence, y = n_DE, color = Tissue)) +
  geom_point(size = 3, alpha = 0.8) +
  # one regression line per tissue; use linewidth instead of size
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed", linewidth = 0.5) +
  # label each point with pair_id
  ggrepel::geom_text_repel(
    aes(label = pair_id),
    size = 3,
    show.legend = FALSE,
    max.overlaps = Inf
  ) +
  theme_bw() +
  xlab("Divergence time / patristic distance") +
  ylab("Number of DE genes (padj < 0.05, |log2FC| ≥ 1)") +
  ggtitle("DE genes vs divergence time across species pairs") +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "right"
  )

# On a cluster without graphics, don't print(p) to screen; just save:
# print(p)  # works in RStudio / local machine

ggsave("DE_genes_vs_divergence_by_tissue_labeled.png",
       p, width = 7, height = 5, dpi = 300)
