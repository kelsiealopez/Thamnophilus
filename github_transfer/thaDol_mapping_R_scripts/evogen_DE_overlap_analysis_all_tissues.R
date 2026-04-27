############################################################
## Overlap: EvoGeneX adaptive genes vs convergently DE
## in same direction in dry species, varying #pairs
############################################################

setwd("/n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/nf_pipeline/thaDol_mapping_test")

library(dplyr)
library(tidyr)
library(tibble)

## 1) Pairs and dry species --------------------------------

pairs_arid <- tribble(
  ~pair_id,             ~spA,     ~spB,     ~dry_species, ~mesic_species,
  "dysSti_vs_sakCri",   "dysSti", "sakCri", "sakCri",     "dysSti",
  "sakLuc_vs_sakCan",   "sakLuc", "sakCan", "sakCan",     "sakLuc",
  "thaAmb_vs_thaPel",   "thaAmb", "thaPel", "thaPel",     "thaAmb",
  "thaBer_vs_thaAtr",   "thaBer", "thaAtr", "thaBer",     "thaAtr",
  "thaRuf_vs_thaTor",   "thaRuf", "thaTor", "thaTor",     "thaRuf"
)

tissues <- c("Brain","Heart","Liver","Muscle","Kidney")
thresholds <- c(2L, 3L, 4L, 5L)  # ≥2, ≥3, ≥4, ≥5 pairs

## 2) Function: per-tissue direction wrt dry ---------------

get_dry_direction <- function(tissue_name) {
  dir_list <- list()
  
  for (i in seq_len(nrow(pairs_arid))) {
    pid <- pairs_arid$pair_id[i]
    spA <- pairs_arid$spA[i]
    spB <- pairs_arid$spB[i]
    dry <- pairs_arid$dry_species[i]
    
    fname <- paste0("DE_OG_selected_", spA, "_vs_", spB, "_", tissue_name, "_kallisto.csv")
    if (!file.exists(fname)) {
      warning("DE table not found: ", fname)
      next
    }
    
    de_tab <- read.csv(fname, row.names = 1, check.names = FALSE)
    de_tab$OG <- rownames(de_tab)
    
    de_sub <- de_tab %>%
      dplyr::filter(
        !is.na(padj),
        padj < 0.05,
        !is.na(log2FoldChange),
        abs(log2FoldChange) >= 1
      ) %>%     # same sig criteria
      dplyr::select(OG, log2FoldChange)
    
    if (nrow(de_sub) == 0) next
    
    directionDry <- if (dry == spA) {
      ifelse(de_sub$log2FoldChange > 0, "dry_up",
             ifelse(de_sub$log2FoldChange < 0, "dry_down", "zero"))
    } else if (dry == spB) {
      ifelse(de_sub$log2FoldChange > 0, "dry_down",
             ifelse(de_sub$log2FoldChange < 0, "dry_up", "zero"))
    } else {
      stop("dry_species not equal to spA or spB for pair ", pid)
    }
    
    dir_list[[pid]] <- data.frame(
      Tissue       = tissue_name,
      OG           = de_sub$OG,
      pair_id      = pid,
      directionDry = directionDry,
      stringsAsFactors = FALSE
    )
  }
  
  if (length(dir_list) == 0) return(NULL)
  bind_rows(dir_list)
}

## 3) Loop over tissues and thresholds ---------------------

for (tiss in tissues) {
  message("=== Tissue: ", tiss, " ===")
  
  # EvoGeneX adaptive OGs for this tissue
  evog_file <- paste0("EvoGeneX_", tiss, "_results.csv")
  if (!file.exists(evog_file)) {
    warning("No EvoGeneX results for ", tiss)
    next
  }
  evog_res <- read.csv(evog_file, stringsAsFactors = FALSE)
  
  adaptive <- evog_res %>%
    dplyr::filter(class == "adaptive", Tissue == tiss) %>%
    dplyr::select(
      OG, Tissue, loglik_BM, loglik_OU, loglik_OU2,
      p_OU_vs_BM, p_OU2_vs_BM, p_OU2_vs_OU
    )
  
  message("  Adaptive OGs: ", nrow(adaptive))
  
  # Direction wrt dry species from per-pair DE
  df_dir <- get_dry_direction(tiss)
  if (is.null(df_dir)) {
    warning("No directional DE info for ", tiss)
    next
  }
  
  counts <- df_dir %>%
    dplyr::filter(directionDry %in% c("dry_up","dry_down")) %>%
    group_by(OG) %>%
    summarise(
      n_dry_up   = sum(directionDry == "dry_up"),
      n_dry_down = sum(directionDry == "dry_down"),
      n_pairs    = n(),
      .groups    = "drop"
    )
  
  for (thr in thresholds) {
    message("  Threshold: ≥", thr, " pairs, same direction")
    
    conv_up   <- counts %>% dplyr::filter(n_dry_up   >= thr)
    conv_down <- counts %>% dplyr::filter(n_dry_down >= thr)
    
    # UP in dry and adaptive
    up_adapt <- adaptive %>%
      inner_join(conv_up, by = "OG")
    
    # DOWN in dry and adaptive
    down_adapt <- adaptive %>%
      inner_join(conv_down, by = "OG")
    
    message("    UP in dry & adaptive:   ", nrow(up_adapt))
    message("    DOWN in dry & adaptive: ", nrow(down_adapt))
    
    ## Optional: annotate with convergent DE summary file if available
    de_conv_file <- paste0("OGs_DE_in_multiple_pairs_", tiss, "_with_refGeneIDs.csv")
    if (file.exists(de_conv_file)) {
      de_conv <- read.csv(de_conv_file, stringsAsFactors = FALSE)
      
      up_adapt_annot <- up_adapt %>%
        inner_join(de_conv, by = c("OG","Tissue"))
      
      down_adapt_annot <- down_adapt %>%
        inner_join(de_conv, by = c("OG","Tissue"))
      
      # Collapse to one row per OG (UP)
      up_OG <- up_adapt_annot %>%
        group_by(Tissue, OG) %>%
        summarise(
          n_pairs_total   = dplyr::first(n_pairs.x),  # from counts
          n_dry_up        = dplyr::first(n_dry_up),
          Pairs           = dplyr::first(Pairs),
          example_gene_id = dplyr::first(gene_id),
          example_galGal  = dplyr::first(galGal_geneID),
          example_taeGut  = dplyr::first(taeGut_geneID),
          .groups = "drop"
        )
      
      # Collapse to one row per OG (DOWN)
      down_OG <- down_adapt_annot %>%
        group_by(Tissue, OG) %>%
        summarise(
          n_pairs_total   = dplyr::first(n_pairs.x),
          n_dry_down      = dplyr::first(n_dry_down),
          Pairs           = dplyr::first(Pairs),
          example_gene_id = dplyr::first(gene_id),
          example_galGal  = dplyr::first(galGal_geneID),
          example_taeGut  = dplyr::first(taeGut_geneID),
          .groups = "drop"
        )
      
      # Write out
      up_file_full  <- paste0(tiss, "_Adaptive_and_ConvUP_in_dry_GE", thr, "pairs_full.csv")
      up_file_og    <- paste0(tiss, "_Adaptive_and_ConvUP_in_dry_GE", thr, "pairs_byOG.csv")
      down_file_full<- paste0(tiss, "_Adaptive_and_ConvDOWN_in_dry_GE", thr, "pairs_full.csv")
      down_file_og  <- paste0(tiss, "_Adaptive_and_ConvDOWN_in_dry_GE", thr, "pairs_byOG.csv")
      
      write.csv(up_adapt_annot,   up_file_full,   row.names = FALSE)
      write.csv(up_OG,            up_file_og,     row.names = FALSE)
      write.csv(down_adapt_annot, down_file_full, row.names = FALSE)
      write.csv(down_OG,          down_file_og,   row.names = FALSE)
    }
  }
}



## ============================================
## GO for adaptive + convergently DE in dry
##   using representative gene_name from og_rep
##   - Brain, Heart, Liver, Muscle: ≥3 pairs
##   - Kidney: ≥2 pairs
##   Human GO (org.Hs.eg.db), BP
## ============================================

setwd("/n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/nf_pipeline/thaDol_mapping_test")

library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)

# 0) If og_rep is not in memory, read it
#    (this file was written earlier by your big script)
if (!exists("og_rep")) {
  og_rep <- read.csv("OGs_DE_in_multiple_pairs_ALL_tissues_repGeneName.csv",
                     stringsAsFactors = FALSE)
}

# 1) Tissue-specific thresholds for # of pairs
tissues <- c("Brain","Heart","Liver","Muscle","Kidney")
thr_by_tissue <- c(Brain = 3L,
                   Heart = 3L,
                   Liver = 3L,
                   Muscle = 3L,
                   Kidney = 2L)

# 2) Read in per-tissue by-OG files (UP + DOWN) at chosen thresholds
all_ogs <- list()

for (tiss in tissues) {
  thr <- thr_by_tissue[[tiss]]
  
  up_file_og   <- paste0(tiss, "_Adaptive_and_ConvUP_in_dry_GE", thr, "pairs_byOG.csv")
  down_file_og <- paste0(tiss, "_Adaptive_and_ConvDOWN_in_dry_GE", thr, "pairs_byOG.csv")
  
  if (!file.exists(up_file_og) && !file.exists(down_file_og)) {
    message("No by-OG files found for ", tiss, " at threshold ≥", thr, " pairs; skipping.")
    next
  }
  
  up_df   <- if (file.exists(up_file_og))   read.csv(up_file_og,   stringsAsFactors = FALSE) else NULL
  down_df <- if (file.exists(down_file_og)) read.csv(down_file_og, stringsAsFactors = FALSE) else NULL
  
  if (!is.null(up_df))   up_df$direction   <- "dry_up"
  if (!is.null(down_df)) down_df$direction <- "dry_down"
  
  comb_tiss <- dplyr::bind_rows(up_df, down_df)
  
  if (nrow(comb_tiss) > 0) {
    all_ogs[[tiss]] <- comb_tiss
  }
}

if (length(all_ogs) == 0) {
  stop("No per-tissue by-OG files were found. Check file names and thresholds.")
}

combined_ogs <- dplyr::bind_rows(all_ogs)

# combined_ogs columns now include:
# "Tissue","OG","n_pairs_total","n_dry_up","Pairs",
# "example_gene_id","example_galGal","example_taeGut","direction","n_dry_down"

# 3) Attach representative gene_name from og_rep
#    og_rep has columns: Tissue, OG, n_pairs, Pairs, gene_name, name_source
combined_ogs_rep <- combined_ogs %>%
  dplyr::left_join(
    og_rep %>% dplyr::select(Tissue, OG, gene_name, name_source),
    by = c("Tissue","OG")
  )

# Quick sanity check
head(combined_ogs_rep[, c("Tissue","OG","gene_name","name_source")])

# 4) Get unique gene_name values
gene_symbols <- combined_ogs_rep %>%
  dplyr::filter(!is.na(gene_name) & gene_name != "") %>%
  dplyr::pull(gene_name) %>%
  unique()

length(gene_symbols)
head(gene_symbols, 20)

# 5) Map gene_name -> human ENTREZID (assuming gene_name ~ SYMBOL)
gene_map <- clusterProfiler::bitr(
  gene_symbols,
  fromType = "SYMBOL",
  toType   = "ENTREZID",
  OrgDb    = org.Hs.eg.db
)

nrow(gene_map)

if (nrow(gene_map) == 0) {
  stop("No gene_name entries could be mapped to human ENTREZ IDs. 
       Check what gene_name actually contains.")
}

# 6) Run GO enrichment (BP) on all mapped genes
ego_all <- enrichGO(
  gene          = gene_map$ENTREZID,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.5,
  qvalueCutoff  = 0.2,
  readable      = TRUE
)

# 7) Plot in RStudio (Plots pane)
if (is.null(ego_all) || nrow(as.data.frame(ego_all)) == 0) {
  message("No significant GO terms found at the chosen cutoffs.")
} else {
  p <- dotplot(ego_all, showCategory = 20) +
    ggtitle("GO BP – adaptive & convergently DE (all tissues, human)")
  print(p)
}

# ego_all now contains the GO results for further inspection:
# head(as.data.frame(ego_all))


gene_map$SYMBOL







#
#
#
#
# Build a per‑tissue / per‑threshold overlap summary
#
#
#
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)

setwd("/n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/nf_pipeline/thaDol_mapping_test")

pairs_arid <- tribble(
  ~pair_id,             ~spA,     ~spB,     ~dry_species, ~mesic_species,
  "dysSti_vs_sakCri",   "dysSti", "sakCri", "sakCri",     "dysSti",
  "sakLuc_vs_sakCan",   "sakLuc", "sakCan", "sakCan",     "sakLuc",
  "thaAmb_vs_thaPel",   "thaAmb", "thaPel", "thaPel",     "thaAmb",
  "thaBer_vs_thaAtr",   "thaBer", "thaAtr", "thaBer",     "thaAtr",
  "thaRuf_vs_thaTor",   "thaRuf", "thaTor", "thaTor",     "thaRuf"
)

tissues    <- c("Brain","Heart","Liver","Muscle","Kidney")
thresholds <- c(2L, 3L, 4L, 5L)

## 1) Direction wrt dry species: reuse your get_dry_direction() ---------

get_dry_direction <- function(tissue_name) {
  dir_list <- list()
  
  for (i in seq_len(nrow(pairs_arid))) {
    pid <- pairs_arid$pair_id[i]
    spA <- pairs_arid$spA[i]
    spB <- pairs_arid$spB[i]
    dry <- pairs_arid$dry_species[i]
    
    fname <- paste0("DE_OG_selected_", spA, "_vs_", spB, "_", tissue_name, "_kallisto.csv")
    if (!file.exists(fname)) {
      warning("DE table not found: ", fname)
      next
    }
    de_tab <- read.csv(fname, row.names = 1, check.names = FALSE)
    de_tab$OG <- rownames(de_tab)
    
    de_sub <- de_tab %>%
      dplyr::filter(
        !is.na(padj),
        padj < 0.05,
        !is.na(log2FoldChange),
        abs(log2FoldChange) >= 1
      ) %>%
      dplyr::select(OG, log2FoldChange)
    
    if (nrow(de_sub) == 0) next
    
    directionDry <- if (dry == spA) {
      ifelse(de_sub$log2FoldChange > 0, "dry_up",
             ifelse(de_sub$log2FoldChange < 0, "dry_down", "zero"))
    } else if (dry == spB) {
      ifelse(de_sub$log2FoldChange > 0, "dry_down",
             ifelse(de_sub$log2FoldChange < 0, "dry_up", "zero"))
    } else {
      stop("dry_species not equal to spA or spB for pair ", pid)
    }
    
    dir_list[[pid]] <- data.frame(
      Tissue       = tissue_name,
      OG           = de_sub$OG,
      pair_id      = pid,
      directionDry = directionDry,
      stringsAsFactors = FALSE
    )
  }
  
  if (length(dir_list) == 0) return(NULL)
  bind_rows(dir_list)
}

## 2) Build summary of overlaps for each tissue & threshold --------------

overlap_summary <- list()

for (tiss in tissues) {
  message("=== Tissue: ", tiss, " ===")
  
  evog_file <- paste0("EvoGeneX_", tiss, "_results.csv")
  if (!file.exists(evog_file)) {
    warning("No EvoGeneX results for ", tiss)
    next
  }
  evog_res <- read.csv(evog_file, stringsAsFactors = FALSE)
  
  adaptive <- evog_res %>%
    dplyr::filter(class == "adaptive", Tissue == tiss) %>%
    dplyr::select(OG, Tissue)
  n_adapt <- nrow(adaptive)
  
  df_dir <- get_dry_direction(tiss)
  if (is.null(df_dir)) {
    warning("No directional DE info for ", tiss)
    next
  }
  
  counts <- df_dir %>%
    dplyr::filter(directionDry %in% c("dry_up","dry_down")) %>%
    group_by(OG) %>%
    summarise(
      n_dry_up   = sum(directionDry == "dry_up"),
      n_dry_down = sum(directionDry == "dry_down"),
      n_pairs    = n(),
      .groups    = "drop"
    )
  
  for (thr in thresholds) {
    conv_up   <- counts %>% dplyr::filter(n_dry_up   >= thr)
    conv_down <- counts %>% dplyr::filter(n_dry_down >= thr)
    
    n_conv_up   <- nrow(conv_up)
    n_conv_down <- nrow(conv_down)
    
    up_adapt   <- adaptive %>% inner_join(conv_up,   by = "OG")
    down_adapt <- adaptive %>% inner_join(conv_down, by = "OG")
    
    overlap_summary[[length(overlap_summary) + 1]] <- tibble(
      Tissue          = tiss,
      Threshold_pairs = thr,
      n_adaptive      = n_adapt,
      n_conv_up       = n_conv_up,
      n_conv_down     = n_conv_down,
      n_adapt_and_up  = nrow(up_adapt),
      n_adapt_and_down= nrow(down_adapt)
    )
  }
}

overlap_summary_df <- bind_rows(overlap_summary)

# Write summary table
write.csv(overlap_summary_df,
          "Adaptive_vs_ConvergentDE_overlap_summary_byTissue.csv",
          row.names = FALSE)

# Long-format for plotting
overlap_long <- overlap_summary_df %>%
  tidyr::pivot_longer(
    cols = c(n_adapt_and_up, n_adapt_and_down),
    names_to = "Direction",
    values_to = "Count"
  ) %>%
  mutate(
    Direction = dplyr::recode(Direction,
                              n_adapt_and_up   = "Adaptive & dry-up",
                              n_adapt_and_down = "Adaptive & dry-down")
  )

ggplot(overlap_long,
       aes(x = Threshold_pairs, y = Count, color = Direction)) +
  geom_line() +
  geom_point() +
  facet_wrap(~ Tissue, scales = "free_y") +
  xlab("Minimum number of dry-direction pairs (threshold)") +
  ylab("# OGs") +
  theme_bw() +
  theme(legend.position = "bottom")




library(dplyr)
library(tidyr)
library(tibble)

setwd("/n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/nf_pipeline/thaDol_mapping_test")

# If og_rep not in memory, load it
if (!exists("og_rep")) {
  og_rep <- read.csv("OGs_DE_in_multiple_pairs_ALL_tissues_repGeneName.csv",
                     stringsAsFactors = FALSE)
}

tissues    <- c("Brain","Heart","Liver","Muscle","Kidney")
thresholds <- c(2L, 3L, 4L, 5L)

all_flag_tables <- list()

for (tiss in tissues) {
  message("=== Building flag table for tissue: ", tiss, " ===")
  
  # 1) OU-adaptive OGs from EvoGeneX
  evog_file <- paste0("EvoGeneX_", tiss, "_results.csv")
  if (!file.exists(evog_file)) {
    warning("  No EvoGeneX results for ", tiss)
    next
  }
  evog_res <- read.csv(evog_file, stringsAsFactors = FALSE)
  
  adaptive_ogs <- evog_res %>%
    dplyr::filter(class == "adaptive", Tissue == tiss) %>%
    dplyr::pull(OG) %>%
    unique()
  
  # 2) Convergent DE direction wrt dry species
  df_dir <- get_dry_direction(tiss)
  if (is.null(df_dir)) {
    warning("  No directional DE info for ", tiss)
    next
  }
  
  counts <- df_dir %>%
    dplyr::filter(directionDry %in% c("dry_up","dry_down")) %>%
    dplyr::group_by(OG) %>%
    dplyr::summarise(
      n_dry_up   = sum(directionDry == "dry_up"),
      n_dry_down = sum(directionDry == "dry_down"),
      .groups    = "drop"
    )
  
  # 3) Start flag table with all OGs that appear in either adaptive or counts
  og_all <- union(adaptive_ogs, counts$OG)
  if (length(og_all) == 0) {
    warning("  No OGs to annotate for ", tiss)
    next
  }
  
  flag_df <- tibble(
    Tissue      = tiss,
    OG          = og_all,
    is_adaptive = as.integer(OG %in% adaptive_ogs)
  )
  
  # 4) Add counts (n_dry_up / n_dry_down), fill NAs with 0
  flag_df <- flag_df %>%
    dplyr::left_join(counts, by = "OG") %>%
    dplyr::mutate(
      n_dry_up   = ifelse(is.na(n_dry_up),   0L, n_dry_up),
      n_dry_down = ifelse(is.na(n_dry_down), 0L, n_dry_down)
    )
  
  # 5) Threshold-based flags
  for (thr in thresholds) {
    up_col   <- paste0("dry_up_ge",   thr)
    down_col <- paste0("dry_down_ge", thr)
    flag_df[[up_col]]   <- as.integer(flag_df$n_dry_up   >= thr)
    flag_df[[down_col]] <- as.integer(flag_df$n_dry_down >= thr)
  }
  
  # 6) Attach representative gene_name from og_rep
  flag_df <- flag_df %>%
    dplyr::left_join(
      og_rep %>% dplyr::select(Tissue, OG, gene_name, name_source),
      by = c("Tissue","OG")
    )
  
  # Reorder columns for readability
  flag_df <- flag_df %>%
    dplyr::select(
      Tissue, OG, gene_name, name_source,
      is_adaptive,
      n_dry_up, n_dry_down,
      dplyr::starts_with("dry_up_ge"),
      dplyr::starts_with("dry_down_ge")
    )
  
  # Store and write per-tissue file
  all_flag_tables[[tiss]] <- flag_df
  
  out_tiss <- paste0("OG_flags_", tiss, "_adaptive_and_convergentDry.csv")
  write.csv(flag_df, out_tiss, row.names = FALSE)
  message("  Wrote per-tissue flag table: ", out_tiss)
}

# 7) Combined table across tissues (optional)
if (length(all_flag_tables) > 0) {
  flags_all <- dplyr::bind_rows(all_flag_tables)
  write.csv(flags_all,
            "OG_flags_ALL_tissues_adaptive_and_convergentDry.csv",
            row.names = FALSE)
  message("Wrote combined flags table: OG_flags_ALL_tissues_adaptive_and_convergentDry.csv")
}








#
# also GO on only adaptive 
#
#



library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

# if og_rep not loaded:
if (!exists("og_rep")) {
  og_rep <- read.csv("OGs_DE_in_multiple_pairs_ALL_tissues_repGeneName.csv",
                     stringsAsFactors = FALSE)
}

run_go_adaptive_only <- function(tiss,
                                 OrgDb = org.Hs.eg.db,
                                 ont   = "BP",
                                 pcut  = 0.05) {
  evog_file <- paste0("EvoGeneX_", tiss, "_results.csv")
  if (!file.exists(evog_file)) {
    warning("No EvoGeneX results for ", tiss)
    return(NULL)
  }
  evog_res <- read.csv(evog_file, stringsAsFactors = FALSE)
  
  adaptive_ogs <- evog_res %>%
    dplyr::filter(class == "adaptive", Tissue == tiss) %>%
    dplyr::pull(OG) %>%
    unique()
  
  if (length(adaptive_ogs) < 5) {
    warning("Too few adaptive OGs for ", tiss, " (", length(adaptive_ogs), ")")
    return(NULL)
  }
  
  # map OG -> gene_name via og_rep
  genes_tiss <- og_rep %>%
    dplyr::filter(Tissue == tiss, OG %in% adaptive_ogs) %>%
    dplyr::filter(!is.na(gene_name) & gene_name != "") %>%
    dplyr::pull(gene_name) %>%
    unique()
  
  message("Tissue ", tiss, " – adaptive-only GO with ", length(genes_tiss), " gene symbols.")
  
  if (length(genes_tiss) < 5) {
    warning("Too few gene symbols for ", tiss, " after mapping OGs.")
    return(NULL)
  }
  
  gmap <- tryCatch(
    bitr(genes_tiss,
         fromType = "SYMBOL",
         toType   = "ENTREZID",
         OrgDb    = OrgDb),
    error = function(e) {
      warning("bitr failed for ", tiss, ": ", e$message)
      return(NULL)
    }
  )
  if (is.null(gmap) || nrow(gmap) == 0) {
    warning("No symbols converted to Entrez for ", tiss)
    return(NULL)
  }
  
  ego <- tryCatch(
    enrichGO(
      gene          = gmap$ENTREZID,
      OrgDb         = OrgDb,
      keyType       = "ENTREZID",
      ont           = ont,
      pAdjustMethod = "BH",
      pvalueCutoff  = pcut,
      qvalueCutoff  = pcut,
      readable      = TRUE
    ),
    error = function(e) {
      warning("enrichGO failed for ", tiss, ": ", e$message)
      return(NULL)
    }
  )
  
  if (is.null(ego) || nrow(as.data.frame(ego)) == 0) {
    warning("No enriched GO terms for ", tiss, " adaptive-only.")
    return(NULL)
  }
  
  ego_df <- as.data.frame(ego)
  out_tab <- paste0("GO_BP_Hs_", tiss, "_ADAPTIVE_only_enrichment.csv")
  write.csv(ego_df, out_tab, row.names = FALSE)
  message("  Wrote adaptive-only GO table: ", out_tab)
  
  p_dot <- dotplot(ego, showCategory = 20) +
    ggtitle(paste("GO BP –", tiss, " (adaptive-only, human)")) +
    theme_bw() +
    theme(
      plot.title   = element_text(size = 11),
      axis.text.x  = element_text(size = 8),
      axis.text.y  = element_text(size = 8),
      axis.title.x = element_text(size = 9),
      axis.title.y = element_text(size = 9)
    )
  
  print(p_dot)
  out_dot <- paste0("GO_BP_Hs_", tiss, "_ADAPTIVE_only_dotplot.png")
  ggsave(out_dot, p_dot, width = 8, height = 6, dpi = 300)
  message("  Wrote adaptive-only GO dotplot: ", out_dot)
  
  invisible(ego)
}

# run for all tissues
for (tiss in tissues) {
  message("\n=== Adaptive-only GO for ", tiss, " ===")
  run_go_adaptive_only(tiss)
}





############################################################
## Global GO: all OU-adaptive OGs combined (human, BP)
############################################################

library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)

setwd("/n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/nf_pipeline/thaDol_mapping_test")

# ensure og_rep is loaded
if (!exists("og_rep")) {
  og_rep <- read.csv("OGs_DE_in_multiple_pairs_ALL_tissues_repGeneName.csv",
                     stringsAsFactors = FALSE)
}

# if tissues vector not in scope, define it:
if (!exists("tissues")) {
  tissues <- c("Brain","Heart","Liver","Muscle","Kidney")
}

run_global_go_adaptive_only <- function(tissues,
                                        OrgDb = org.Hs.eg.db,
                                        ont   = "BP",
                                        pcut  = 0.05) {
  all_adaptive_ogs <- c()
  
  # 1) collect adaptive OGs from each tissue
  for (tiss in tissues) {
    evog_file <- paste0("EvoGeneX_", tiss, "_results.csv")
    if (!file.exists(evog_file)) {
      warning("No EvoGeneX results for ", tiss)
      next
    }
    evog_res <- read.csv(evog_file, stringsAsFactors = FALSE)
    
    adaptive_ogs_t <- evog_res %>%
      dplyr::filter(class == "adaptive", Tissue == tiss) %>%
      dplyr::pull(OG) %>%
      unique()
    
    all_adaptive_ogs <- c(all_adaptive_ogs, adaptive_ogs_t)
  }
  
  all_adaptive_ogs <- unique(all_adaptive_ogs)
  if (length(all_adaptive_ogs) < 5) {
    warning("Too few adaptive OGs overall (", length(all_adaptive_ogs),
            ") to run global GO.")
    return(NULL)
  }
  
  message("Global adaptive-only GO: ", length(all_adaptive_ogs), " OGs total across tissues.")
  
  # 2) Map OG -> gene_name via og_rep (any tissue)
  genes_all <- og_rep %>%
    dplyr::filter(OG %in% all_adaptive_ogs) %>%
    dplyr::filter(!is.na(gene_name) & gene_name != "") %>%
    dplyr::pull(gene_name) %>%
    unique()
  
  message("Global adaptive-only GO: ", length(genes_all), " unique gene_name entries.")
  
  if (length(genes_all) < 5) {
    warning("Too few gene symbols after mapping OGs globally.")
    return(NULL)
  }
  
  # 3) SYMBOL -> ENTREZID
  gmap <- tryCatch(
    bitr(genes_all,
         fromType = "SYMBOL",
         toType   = "ENTREZID",
         OrgDb    = OrgDb),
    error = function(e) {
      warning("Global adaptive-only GO: bitr failed: ", e$message)
      return(NULL)
    }
  )
  if (is.null(gmap) || nrow(gmap) == 0) {
    warning("Global adaptive-only GO: no symbols converted to Entrez.")
    return(NULL)
  }
  
  # 4) enrichGO
  ego <- tryCatch(
    enrichGO(
      gene          = gmap$ENTREZID,
      OrgDb         = OrgDb,
      keyType       = "ENTREZID",
      ont           = ont,
      pAdjustMethod = "BH",
      pvalueCutoff  = pcut,
      qvalueCutoff  = pcut,
      readable      = TRUE
    ),
    error = function(e) {
      warning("Global adaptive-only GO: enrichGO failed: ", e$message)
      return(NULL)
    }
  )
  
  if (is.null(ego) || nrow(as.data.frame(ego)) == 0) {
    warning("Global adaptive-only GO: no enriched GO terms.")
    return(NULL)
  }
  
  ego_df <- as.data.frame(ego)
  out_tab <- "GO_BP_Hs_ALL_tissues_ADAPTIVE_only_enrichment.csv"
  write.csv(ego_df, out_tab, row.names = FALSE)
  message("Global adaptive-only GO: wrote table ", out_tab)
  
  # 5) Dotplot
  p_dot <- dotplot(ego, showCategory = 20) +
    ggtitle("GO BP – ALL tissues (adaptive-only, human)") +
    theme_bw() +
    theme(
      plot.title   = element_text(size = 11),
      axis.text.x  = element_text(size = 8),
      axis.text.y  = element_text(size = 8),
      axis.title.x = element_text(size = 9),
      axis.title.y = element_text(size = 9)
    )
  
  print(p_dot)
  
  out_dot <- "GO_BP_Hs_ALL_tissues_ADAPTIVE_only_dotplot.png"
  ggsave(out_dot, p_dot, width = 8, height = 6, dpi = 300)
  message("Global adaptive-only GO: wrote dotplot ", out_dot)
  
  invisible(ego)
}

# Run the global adaptive-only GO
ego_global_adaptive <- run_global_go_adaptive_only(tissues)





#### venn diagrams
library(dplyr)
library(eulerr)
library(RColorBrewer)

setwd("/n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/nf_pipeline/thaDol_mapping_test")

flags_all <- read.csv("OG_flags_ALL_tissues_adaptive_and_convergentDry.csv",
                      stringsAsFactors = FALSE)

tissues    <- c("Brain","Heart","Liver","Muscle","Kidney")
thresholds <- c(2L, 3L, 4L, 5L)

# Helper: draws one Euler diagram directly to the current device
draw_euler_panel <- function(df_tiss, tiss, thr) {
  thr_col_up   <- paste0("dry_up_ge",   thr)
  thr_col_down <- paste0("dry_down_ge", thr)
  
  if (!all(c(thr_col_up, thr_col_down) %in% colnames(df_tiss))) {
    warning("Columns ", thr_col_up, " or ", thr_col_down, " missing for ", tiss)
    return(invisible(NULL))
  }
  
  adaptive_ogs <- df_tiss %>%
    dplyr::filter(is_adaptive == 1) %>%
    dplyr::pull(OG) %>%
    unique()
  
  up_ogs <- df_tiss %>%
    dplyr::filter(.data[[thr_col_up]] == 1) %>%
    dplyr::pull(OG) %>%
    unique()
  
  down_ogs <- df_tiss %>%
    dplyr::filter(.data[[thr_col_down]] == 1) %>%
    dplyr::pull(OG) %>%
    unique()
  
  sets <- list(
    Adaptive   = adaptive_ogs,
    "Dry-up"   = up_ogs,
    "Dry-down" = down_ogs
  )
  
  # Drop empty sets
  sets <- sets[vapply(sets, function(x) length(x) > 0, logical(1))]
  if (length(sets) < 2) {
    plot.new()
    title(main = paste0(tiss, " (≥", thr, "): <2 non-empty sets"))
    return(invisible(NULL))
  }
  
  fit  <- eulerr::euler(sets)
  cols <- RColorBrewer::brewer.pal(max(3, length(sets)), "Set3")
  
  plot(
    fit,
    fills = list(fill = cols[seq_along(sets)], alpha = 0.6),
    labels = list(font = 2),
    edges = list(lwd = 1),
    quantities = list(type = "counts", font = 1),
    main = paste0("≥", thr, " pairs")
  )
}

# For each tissue, make a 2x2 grid (thr = 2,3,4,5) in one PDF
for (tiss in tissues) {
  df_tiss <- flags_all %>%
    dplyr::filter(Tissue == tiss)
  
  if (nrow(df_tiss) == 0) {
    warning("No flags for tissue: ", tiss)
    next
  }
  
  pdf_file <- paste0("Euler_Adaptive_vs_ConvDE_", tiss, "_grid.pdf")
  pdf(pdf_file, width = 8, height = 8)
  par(mfrow = c(2, 2), mar = c(4, 4, 4, 2))  # 2x2 grid
  
  for (thr in thresholds) {
    draw_euler_panel(df_tiss, tiss, thr)
  }
  
  dev.off()
  message("Wrote grid PDF: ", pdf_file)
}



# ilters by the tissue‑specific thresholds
# Brain/Heart/Liver/Muscle: ≥3 pairs
# Kidney: ≥2 pairs

# unique gene symbols that are adaptive + convergently DE
genes_adaptive_conv <- combined_ogs_rep %>%
  dplyr::filter(!is.na(gene_name), gene_name != "") %>%
  dplyr::pull(gene_name) %>%
  unique()

length(genes_adaptive_conv)
head(genes_adaptive_conv)

# if you want to save them:
write.table(genes_adaptive_conv,
            file = "Adaptive_and_ConvergentDry_gene_names.txt",
            quote = FALSE, row.names = FALSE, col.names = FALSE)
