setwd("/n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/nf_pipeline/thaDol_mapping_test")

library(dplyr)
library(VennDiagram)
library(grid)   # for grid.draw

flags_all <- read.csv("OG_flags_ALL_tissues_adaptive_and_convergentDry.csv",
                      stringsAsFactors = FALSE)

tissues <- c("Brain","Heart","Liver","Muscle","Kidney")
thresholds <- c(2, 3, 4, 5)

for (thr in thresholds) {
  up_col   <- paste0("dry_up_ge", thr)
  down_col <- paste0("dry_down_ge", thr)
  
  for (tiss in tissues) {
    df_tiss <- flags_all %>% filter(Tissue == tiss)
    
    # Safety: skip if those columns don't exist for this file
    if (!all(c(up_col, down_col) %in% colnames(df_tiss))) {
      warning("Missing ", up_col, " or ", down_col, " for ", tiss,
              " at threshold ", thr, "; skipping.")
      next
    }
    
    adaptive <- df_tiss$OG[df_tiss$is_adaptive == 1]
    conv     <- df_tiss$OG[df_tiss[[up_col]] == 1 | df_tiss[[down_col]] == 1]
    
    sets <- list(
      Adaptive      = adaptive,
      ConvergentDry = conv
    )
    
    venn.grob <- venn.diagram(
      x = sets,
      filename = NULL,
      fill = c("#63285C", "#DC9471"),
      alpha = c(0.6, 0.6),
      cex = 1.5,
      cat.cex = 1.5,
      cat.pos = c(-20, 20),
      cat.dist = 0.05,
      main = paste0(tiss, " (≥", thr, " dry-direction pairs)"),
      main.cex = 1.4
    )
    
    pdf(paste0("VennVD_Adaptive_vs_ConvergentDry_",
               tiss, "_ge", thr, "pairs.pdf"),
        width = 5, height = 5)
    grid.draw(venn.grob)
    dev.off()
  }
}






setwd("/n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/nf_pipeline/thaDol_mapping_test")

library(dplyr)

evog_files <- list.files(pattern = "^EvoGeneX_.*_results\\.csv$")

evog_all <- bind_rows(
  lapply(evog_files, function(f) read.csv(f, stringsAsFactors = FALSE))
)

adaptive_ogs_all <- evog_all %>%
  filter(class == "adaptive") %>%
  distinct(OG)

n_unique_adaptive_OGs <- nrow(adaptive_ogs_all)
n_unique_adaptive_OGs



adaptive_by_tissue <- evog_all %>%
  filter(class == "adaptive") %>%
  group_by(Tissue) %>%
  summarise(n_adaptive_unique = n_distinct(OG), .groups = "drop")

adaptive_by_tissue
n_unique_adaptive_OGs  # global unique across all tissues




# If not already in memory, reload:
combined_ogs_rep <- read.csv("Adaptive_and_ConvergentDry_gene_names.txt", 
                             stringsAsFactors = FALSE)  # this is genes, not OGs

# Better: use combined_ogs (if you saved it), or rebuild from *_GE* byOG files.
# Assuming combined_ogs is still in your session, you can do:
unique_overlap_OGs <- combined_ogs %>%
  distinct(OG)

n_unique_overlap_OGs <- nrow(unique_overlap_OGs)
n_unique_overlap_OGs









setwd("/n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/nf_pipeline/thaDol_mapping_test")

library(dplyr)

flags_all <- read.csv("OG_flags_ALL_tissues_adaptive_and_convergentDry.csv",
                      stringsAsFactors = FALSE)

## choose threshold in # of dry-direction pairs
threshold <- 3  # change to 2, 4, or 5 as needed

## mark OGs that are convergently DE at this threshold
conv_flags <- flags_all %>%
  mutate(
    is_convergent = (n_dry_up   >= threshold) |
      (n_dry_down >= threshold)
  ) %>%
  filter(is_convergent)

## 1) Unique convergently DE OGs across all tissues
n_unique_conv_OGs <- n_distinct(conv_flags$OG)
n_unique_conv_OGs

## 2) Unique convergently DE gene names across all tissues
conv_genes <- conv_flags %>%
  filter(!is.na(gene_name), gene_name != "") %>%
  distinct(gene_name)

n_unique_conv_genes <- nrow(conv_genes)
n_unique_conv_genes

## Optional: see per-tissue counts as well
conv_by_tissue <- conv_flags %>%
  group_by(Tissue) %>%
  summarise(
    n_conv_OGs    = n_distinct(OG),
    n_conv_genes  = n_distinct(gene_name[!is.na(gene_name) & gene_name != ""]),
    .groups = "drop"
  )

conv_by_tissue

