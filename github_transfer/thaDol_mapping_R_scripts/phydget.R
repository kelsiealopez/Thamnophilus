############################################################
# prepare data for phydget – using *.kallisto.merged.OG_all9_counts.tsv
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
  
  gene_counts_list[[run]] <- df
}

stopifnot(length(gene_counts_list) > 0)

## Merge on gene_id
merged <- Reduce(function(x, y) merge(x, y, by = "gene_id", all = FALSE),
                 gene_counts_list)

cat("Merged matrix has", nrow(merged), "genes and", ncol(merged) - 1, "samples\n")

rownames(merged) <- merged$gene_id
merged$gene_id <- NULL
cts <- merged   # transcripts × samples

## Collapse multiple transcripts per OG to a single row if needed
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
## EXCLUDE SAMPLES HERE (same as DE analysis)
## =========================================
samples_to_drop <- c("thaPal01_Liver")

keep_cols <- setdiff(colnames(cts), samples_to_drop)
cts <- cts[, keep_cols, drop = FALSE]
## =========================================

## Build coldata from sample names
samples    <- colnames(cts)
tissue     <- sub(".*_", "", samples)           # after last "_"
id_species <- sub("_[^_]*$", "", samples)       # before last "_"
species    <- sub("[0-9]+.*$", "", id_species)  # leading letters

coldata <- data.frame(
  Sample  = samples,
  Species = species,
  Tissue  = tissue,
  stringsAsFactors = FALSE
)

## Fix thaCap -> thaDol in Species only
coldata$Species[coldata$Species == "thaCap"] <- "thaDol"

stopifnot(all(colnames(cts) == coldata$Sample))

############################################################
## 2) Write per‑tissue PhyDGET input tables
############################################################
## PhyDGET expects:
##  - first column: OG (gene ID)
##  - remaining columns: samples (counts)
##  - tab‑separated text, with header

write_phydget_counts_for_tissue <- function(tiss, cts_mat, coldata_df) {
  idx <- which(coldata_df$Tissue == tiss)
  if (length(idx) == 0) {
    warning("No samples for tissue ", tiss)
    return(NULL)
  }
  cts_tiss <- cts_mat[, idx, drop = FALSE]
  
  # Add OG column as first column
  out <- data.frame(
    OG = rownames(cts_tiss),
    cts_tiss,
    check.names = FALSE
  )
  
  fn <- paste0("Phydget_counts_", tiss, "_OG_all9.tsv")
  write.table(out, fn, sep = "\t", quote = FALSE, row.names = FALSE)
  message("Wrote ", fn, " with ",
          nrow(out), " OGs and ", ncol(out) - 1, " samples.")
  invisible(fn)
}

## Tissues you want to run in PhyDGET
tissues_to_run <- c("Brain", "Heart", "Kidney", "Liver", "Muscle")

written_files <- lapply(tissues_to_run, write_phydget_counts_for_tissue,
                        cts_mat = cts, coldata_df = coldata)

written_files




# prepare ultrametric trees!

############################################################
# Make Thamnophilus species trees ultrametric for PhyDGET
############################################################
library(ape)

setwd("/n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/fastoma")

## List of input trees (Newick) and desired output names
tree_files <- c(
  KIDNEY = "species_tree_KIDNEY_13_sp.nwk",
  MUSCLE = "species_tree_MUSCLE_13_sp.nwk",
  ALL14  = "species_tree_14_sp.nwk"
)

## Directory for ultrametric output
out_dir <- "ultrametric_trees"
if (!dir.exists(out_dir)) dir.create(out_dir)

make_ultrametric_and_write <- function(infile, out_prefix) {
  cat("Reading tree:", infile, "\n")
  tr <- read.tree(infile)
  
  cat("  is.ultrametric before:", is.ultrametric(tr), "\n")
  
  ## Use chronos to enforce ultrametricity
  cat("  Running chronos() to enforce ultrametricity...\n")
  tr_ultra <- chronos(tr, lambda = 1)  # adjust lambda if desired
  
  cat("  is.ultrametric after:", is.ultrametric(tr_ultra), "\n")
  
  ## Write both Newick and Nexus versions
  out_nwk <- file.path(out_dir, paste0(out_prefix, "_ultrametric.nwk"))
  out_nex <- file.path(out_dir, paste0(out_prefix, "_ultrametric.nex"))
  
  write.tree(tr_ultra, file = out_nwk)
  write.nexus(tr_ultra, file = out_nex)
  
  cat("  Wrote ultrametric Newick to:", out_nwk, "\n")
  cat("  Wrote ultrametric Nexus  to:", out_nex, "\n\n")
  
  invisible(tr_ultra)
}

## Run for all three trees
ultra_trees <- mapply(
  FUN        = make_ultrametric_and_write,
  infile     = tree_files,
  out_prefix = names(tree_files),
  SIMPLIFY   = FALSE
)

ultra_trees







# list out thhe samples for phydget

############################################################
# Build PhyDGET --samples lines for Kidney automatically
############################################################

setwd("/n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/nf_pipeline/thaDol_mapping_test")

# 1) Read just the header of the Kidney PhyDGET counts file
kidney_header <- read.delim("Phydget_counts_Kidney_OG_all9.tsv",
                            nrows = 0, check.names = FALSE)
kidney_samples <- colnames(kidney_header)[-1]  # drop "OG"

# 2) Parse species and tissue from sample names
#    Format is e.g. sakCri01_Kidney, thaAtr04_Kidney, etc.
kidney_tissue <- sub(".*_", "", kidney_samples)                 # after last "_"
id_species    <- sub("_[^_]*$", "", kidney_samples)             # before last "_"
kidney_species <- sub("[0-9]+.*$", "", id_species)              # leading letters only

kidney_coldata <- data.frame(
  Sample  = kidney_samples,
  Species = kidney_species,
  Tissue  = kidney_tissue,
  stringsAsFactors = FALSE
)

# Fix thaCap -> thaDol (as in your other scripts)
kidney_coldata$Species[kidney_coldata$Species == "thaCap"] <- "thaDol"

# Sanity check: all tissue labels should be "Kidney"
print(table(kidney_coldata$Tissue))

# 3) For each species, collect its Kidney samples and build --samples line
species_order <- sort(unique(kidney_coldata$Species))

sample_lines <- sapply(species_order, function(sp) {
  sp_samples <- kidney_coldata$Sample[kidney_coldata$Species == sp]
  # Ensure there is at least one sample
  if (length(sp_samples) == 0) return(NA_character_)
  paste0("--samples ", sp, ":", paste(sp_samples, collapse = ","))
})

# Drop any NA (shouldn’t be any, but just in case)
sample_lines <- sample_lines[!is.na(sample_lines)]

# 4) Print to console so you can copy/paste into the PhyDGET job file
cat("# Auto-generated --samples lines for Kidney\n")
cat(paste(sample_lines, collapse = "\n"), "\n")

# 5) Optionally, write to a file for convenience
writeLines(c("# Auto-generated --samples lines for Kidney",
             sample_lines),
           con = "Phydget_Kidney_samples_block.txt")

message("Wrote Kidney sample block to Phydget_Kidney_samples_block.txt")

