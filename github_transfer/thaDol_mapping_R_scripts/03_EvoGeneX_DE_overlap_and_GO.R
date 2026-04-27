############################################################
## 03_EvoGeneX_DE_overlap_and_GO.R
##
## EvoGeneX post-processing:
## Adaptive vs convergent DE + GO with tissue-specific background
##
## Prerequisites:
##   - 01_DE_and_GO.R has been run, producing:
##       obj_dds_DESeq2_OG_all_tissues.rds
##       obj_og_rep_ALL_OGs_repGeneName.rds
##       obj_og_rep_multiple_pairs_byTissue.rds
##   - 02_EvoGeneX_all_tissues.R has been run, producing:
##       EvoGeneX_<Tissue>_results.csv
##       obj_EvoGeneX_all_tissues_results.rds (optional)
############################################################

setwd("/n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/nf_pipeline/thaDol_mapping_test")

library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(eulerr)
library(RColorBrewer)

############################################################
## 0) Load core objects from previous scripts
############################################################

## dds: DESeqDataSet for OGs across all tissues
if (file.exists("obj_dds_DESeq2_OG_all_tissues.rds")) {
  dds <- readRDS("obj_dds_DESeq2_OG_all_tissues.rds")
} else {
  stop("obj_dds_DESeq2_OG_all_tissues.rds not found. Run 01_DE_and_GO.R first.")
}

## Global OG -> representative gene_name mapping (ALL OGs)
if (file.exists("obj_og_rep_ALL_OGs_repGeneName.rds")) {
  og_rep_all <- readRDS("obj_og_rep_ALL_OGs_repGeneName.rds")
} else if (file.exists("OGs_ALL_repGeneName.csv")) {
  og_rep_all <- read.csv("OGs_ALL_repGeneName.csv", stringsAsFactors = FALSE)
} else {
  stop("Neither obj_og_rep_ALL_OGs_repGeneName.rds nor OGs_ALL_repGeneName.csv found.")
}

## Per-tissue OG -> repGeneName for OGs DE in multiple pairs
if (file.exists("obj_og_rep_multiple_pairs_byTissue.rds")) {
  og_rep <- readRDS("obj_og_rep_multiple_pairs_byTissue.rds")
} else if (file.exists("OGs_DE_in_multiple_pairs_ALL_tissues_repGeneName.csv")) {
  og_rep <- read.csv("OGs_DE_in_multiple_pairs_ALL_tissues_repGeneName.csv",
                     stringsAsFactors = FALSE)
} else {
  stop("Neither obj_og_rep_multiple_pairs_byTissue.rds nor OGs_DE_in_multiple_pairs_ALL_tissues_repGeneName.csv found.")
}

## EvoGeneX tissues vector
tissues <- c("Brain","Heart","Liver","Muscle","Kidney")

## Pairs: dry vs mesic (same as in DE/EvoGeneX scripts)
pairs_arid <- tribble(
  ~pair_id,             ~spA,     ~spB,     ~dry_species, ~mesic_species,
  "dysSti_vs_sakCri",   "dysSti", "sakCri", "sakCri",     "dysSti",
  "sakLuc_vs_sakCan",   "sakLuc", "sakCan", "sakCan",     "sakLuc",
  "thaAmb_vs_thaPel",   "thaAmb", "thaPel", "thaPel",     "thaAmb",
  "thaBer_vs_thaAtr",   "thaBer", "thaAtr", "thaBer",     "thaAtr",
  "thaRuf_vs_thaTor",   "thaRuf", "thaTor", "thaTor",     "thaRuf"
)

############################################################
## Helper functions: background universes
############################################################

## Tissue-specific background based on expressed OGs in dds
get_tissue_background_entrez <- function(tiss,
                                         dds,
                                         og_rep_all,
                                         OrgDb) {
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

## Global background = all expressed OGs in dds
get_global_background_entrez <- function(dds,
                                         og_rep_all,
                                         OrgDb) {
  cts_all <- counts(dds)
  bg_ogs  <- rownames(cts_all)[rowSums(cts_all > 0) > 0]
  if (length(bg_ogs) < 5) {
    warning("Too few background OGs globally (", length(bg_ogs), "); skipping.")
    return(NULL)
  }
  bg_map <- og_rep_all %>%
    dplyr::filter(OG %in% bg_ogs, !is.na(gene_name), gene_name != "")
  bg_symbols <- unique(bg_map$gene_name)
  if (length(bg_symbols) < 5) {
    warning("Too few background symbols globally (", length(bg_symbols), "); skipping.")
    return(NULL)
  }
  bg_conv <- tryCatch(
    bitr(bg_symbols,
         fromType = "SYMBOL",
         toType   = "ENTREZID",
         OrgDb    = OrgDb),
    error = function(e) {
      warning("Global background bitr failed: ", e$message)
      return(NULL)
    }
  )
  if (is.null(bg_conv) || nrow(bg_conv) == 0) {
    warning("No background symbols converted globally.")
    return(NULL)
  }
  unique(bg_conv$ENTREZID)
}

############################################################
## 1) Overlap: EvoGeneX adaptive genes vs convergently DE
##    in same direction in dry species, varying #pairs
############################################################

thresholds <- c(2L, 3L, 4L, 5L)  # ≥2, ≥3, ≥4, ≥5 pairs

## Function: per-tissue DE direction wrt dry species
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

## 1A) Main overlap loop: adaptive vs convergent, write by-OG files

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
    dplyr::select(
      OG, Tissue, loglik_BM, loglik_OU, loglik_OU2,
      p_OU_vs_BM, p_OU2_vs_BM, p_OU2_vs_OU
    )
  
  message("  Adaptive OGs: ", nrow(adaptive))
  
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
    
    up_adapt <- adaptive %>% inner_join(conv_up,   by = "OG")
    down_adapt <- adaptive %>% inner_join(conv_down, by = "OG")
    
    message("    UP in dry & adaptive:   ", nrow(up_adapt))
    message("    DOWN in dry & adaptive: ", nrow(down_adapt))
    
    if (nrow(up_adapt) > 0) {
      up_OG <- up_adapt %>%
        group_by(Tissue, OG) %>%
        summarise(
          n_pairs_total = first(n_pairs),
          n_dry_up      = first(n_dry_up),
          .groups       = "drop"
        )
      
      up_file_full <- paste0(tiss, "_Adaptive_and_ConvUP_in_dry_GE", thr, "pairs_full.csv")
      up_file_og   <- paste0(tiss, "_Adaptive_and_ConvUP_in_dry_GE", thr, "pairs_byOG.csv")
      
      write.csv(up_adapt, up_file_full, row.names = FALSE)
      write.csv(up_OG,     up_file_og,   row.names = FALSE)
    }
    
    if (nrow(down_adapt) > 0) {
      down_OG <- down_adapt %>%
        group_by(Tissue, OG) %>%
        summarise(
          n_pairs_total = first(n_pairs),
          n_dry_down    = first(n_dry_down),
          .groups       = "drop"
        )
      
      down_file_full <- paste0(tiss, "_Adaptive_and_ConvDOWN_in_dry_GE", thr, "pairs_full.csv")
      down_file_og   <- paste0(tiss, "_Adaptive_and_ConvDOWN_in_dry_GE", thr, "pairs_byOG.csv")
      
      write.csv(down_adapt, down_file_full, row.names = FALSE)
      write.csv(down_OG,     down_file_og,   row.names = FALSE)
    }
  }
}

############################################################
## 2) GO for adaptive + convergently DE in dry (all tissues)
##    with expressed-gene background (global universe)
############################################################

if (!exists("combined_ogs_rep")) {
  tissues_go <- c("Brain","Heart","Liver","Muscle","Kidney")
  thr_by_tissue <- c(Brain = 3L,
                     Heart = 3L,
                     Liver = 3L,
                     Muscle = 3L,
                     Kidney = 2L)
  
  all_ogs <- list()
  for (tiss in tissues_go) {
    thr <- thr_by_tissue[[tiss]]
    up_file_og   <- paste0(tiss, "_Adaptive_and_ConvUP_in_dry_GE", thr, "pairs_byOG.csv")
    down_file_og <- paste0(tiss, "_Adaptive_and_ConvDOWN_in_dry_GE", thr, "pairs_byOG.csv")
    
    if (!file.exists(up_file_og) && !file.exists(down_file_og)) next
    
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
    stop("No per-tissue by-OG files found for adaptive+convergent; cannot build combined_ogs_rep.")
  }
  combined_ogs <- dplyr::bind_rows(all_ogs)
  combined_ogs_rep <- combined_ogs %>%
    dplyr::left_join(
      og_rep %>% dplyr::select(Tissue, OG, gene_name, name_source),
      by = c("Tissue","OG")
    )
}

gene_symbols_ac <- combined_ogs_rep %>%
  dplyr::filter(!is.na(gene_name) & gene_name != "") %>%
  dplyr::pull(gene_name) %>%
  unique()

bg_universe_global <- get_global_background_entrez(dds, og_rep_all, org.Hs.eg.db)

if (!is.null(bg_universe_global)) {
  gene_map_ac <- clusterProfiler::bitr(
    gene_symbols_ac,
    fromType = "SYMBOL",
    toType   = "ENTREZID",
    OrgDb    = org.Hs.eg.db
  )
  
  if (nrow(gene_map_ac) > 0) {
    ego_all_ac <- enrichGO(
      gene          = gene_map_ac$ENTREZID,
      universe      = bg_universe_global,
      OrgDb         = org.Hs.eg.db,
      keyType       = "ENTREZID",
      ont           = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.5,
      qvalueCutoff  = 0.2,
      readable      = TRUE
    )
    
    if (!is.null(ego_all_ac) && nrow(as.data.frame(ego_all_ac)) > 0) {
      p_ac <- dotplot(ego_all_ac, showCategory = 20) +
        ggtitle("GO BP – adaptive & convergently DE (all tissues, human)\n(expresssed-gene global background)")
      print(p_ac)
      write.csv(as.data.frame(ego_all_ac),
                "GO_BP_Hs_ALL_tissues_ADAPTIVEplusCONV_DRY_enrichment.csv",
                row.names = FALSE)
      ggsave("GO_BP_Hs_ALL_tissues_ADAPTIVEplusCONV_DRY_dotplot.png",
             p_ac, width = 8, height = 6, dpi = 300)
    } else {
      message("No significant GO terms for adaptive+convergent (global).")
    }
  } else {
    message("No adaptive+convergent symbols could be mapped to ENTREZ.")
  }
}

############################################################
## 3) Per-tissue overlap summary (adaptive vs convergent)
############################################################

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
write.csv(overlap_summary_df,
          "Adaptive_vs_ConvergentDE_overlap_summary_byTissue.csv",
          row.names = FALSE)

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

############################################################
## 4) Per-tissue flag tables (adaptive & convergent)
############################################################

all_flag_tables <- list()

for (tiss in tissues) {
  message("=== Building flag table for tissue: ", tiss, " ===")
  
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
  
  flag_df <- flag_df %>%
    dplyr::left_join(counts, by = "OG") %>%
    dplyr::mutate(
      n_dry_up   = ifelse(is.na(n_dry_up),   0L, n_dry_up),
      n_dry_down = ifelse(is.na(n_dry_down), 0L, n_dry_down)
    )
  
  for (thr in thresholds) {
    up_col   <- paste0("dry_up_ge",   thr)
    down_col <- paste0("dry_down_ge", thr)
    flag_df[[up_col]]   <- as.integer(flag_df$n_dry_up   >= thr)
    flag_df[[down_col]] <- as.integer(flag_df$n_dry_down >= thr)
  }
  
  flag_df <- flag_df %>%
    dplyr::left_join(
      og_rep %>% dplyr::select(Tissue, OG, gene_name, name_source),
      by = c("Tissue","OG")
    ) %>%
    dplyr::select(
      Tissue, OG, gene_name, name_source,
      is_adaptive,
      n_dry_up, n_dry_down,
      dplyr::starts_with("dry_up_ge"),
      dplyr::starts_with("dry_down_ge")
    )
  
  all_flag_tables[[tiss]] <- flag_df
  
  out_tiss <- paste0("OG_flags_", tiss, "_adaptive_and_convergentDry.csv")
  write.csv(flag_df, out_tiss, row.names = FALSE)
  message("  Wrote per-tissue flag table: ", out_tiss)
}

if (length(all_flag_tables) > 0) {
  flags_all <- dplyr::bind_rows(all_flag_tables)
  write.csv(flags_all,
            "OG_flags_ALL_tissues_adaptive_and_convergentDry.csv",
            row.names = FALSE)
  message("Wrote combined flags table: OG_flags_ALL_tissues_adaptive_and_convergentDry.csv")
}

############################################################
## 5) GO on adaptive-only (per tissue) with expr-gene background
############################################################

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
  
  bg_universe <- get_tissue_background_entrez(tiss, dds, og_rep_all, OrgDb)
  if (is.null(bg_universe)) {
    warning("No background universe for ", tiss, " adaptive-only GO.")
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
      universe      = bg_universe,
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
  out_tab <- paste0("GO_BP_Hs_", tiss, "_ADAPTIVE_only_enrichment_exprBG.csv")
  write.csv(ego_df, out_tab, row.names = FALSE)
  message("  Wrote adaptive-only GO table (expr BG): ", out_tab)
  
  p_dot <- dotplot(ego, showCategory = 20) +
    ggtitle(paste("GO BP –", tiss, " (adaptive-only, human)\n(expr. genes background)")) +
    theme_bw() +
    theme(
      plot.title   = element_text(size = 11),
      axis.text.x  = element_text(size = 8),
      axis.text.y  = element_text(size = 8),
      axis.title.x = element_text(size = 9),
      axis.title.y = element_text(size = 9)
    )
  
  print(p_dot)
  out_dot <- paste0("GO_BP_Hs_", tiss, "_ADAPTIVE_only_exprBG_dotplot.png")
  ggsave(out_dot, p_dot, width = 8, height = 6, dpi = 300)
  message("  Wrote adaptive-only GO dotplot (expr BG): ", out_dot)
  
  invisible(ego)
}

for (tiss in tissues) {
  message("\n=== Adaptive-only GO (expr BG) for ", tiss, " ===")
  run_go_adaptive_only(tiss)
}

############################################################
## 6) Global GO: all OU-adaptive OGs combined (expr BG)
############################################################

run_global_go_adaptive_only <- function(tissues,
                                        OrgDb = org.Hs.eg.db,
                                        ont   = "BP",
                                        pcut  = 0.05) {
  all_adaptive_ogs <- c()
  
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
  
  bg_universe_global <- get_global_background_entrez(dds, og_rep_all, OrgDb)
  if (is.null(bg_universe_global)) {
    warning("No global background universe; aborting global GO.")
    return(NULL)
  }
  
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
  
  ego <- tryCatch(
    enrichGO(
      gene          = gmap$ENTREZID,
      universe      = bg_universe_global,
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
  out_tab <- "GO_BP_Hs_ALL_tissues_ADAPTIVE_only_enrichment_exprBG.csv"
  write.csv(ego_df, out_tab, row.names = FALSE)
  message("Global adaptive-only GO (expr BG): wrote table ", out_tab)
  
  p_dot <- dotplot(ego, showCategory = 20) +
    ggtitle("GO BP – ALL tissues (adaptive-only, human)\n(expresssed-gene background)") +
    theme_bw() +
    theme(
      plot.title   = element_text(size = 11),
      axis.text.x  = element_text(size = 8),
      axis.text.y  = element_text(size = 8),
      axis.title.x = element_text(size = 9),
      axis.title.y = element_text(size = 9)
    )
  
  print(p_dot)
  
  out_dot <- "GO_BP_Hs_ALL_tissues_ADAPTIVE_only_exprBG_dotplot.png"
  ggsave(out_dot, p_dot, width = 8, height = 6, dpi = 300)
  message("Global adaptive-only GO (expr BG): wrote dotplot ", out_dot)
  
  invisible(ego)
}

ego_global_adaptive_exprBG <- run_global_go_adaptive_only(tissues)

############################################################
## 7) Euler diagrams of Adaptive vs ConvDry per tissue/threshold
############################################################

flags_all <- read.csv("OG_flags_ALL_tissues_adaptive_and_convergentDry.csv",
                      stringsAsFactors = FALSE)

tissues    <- c("Brain","Heart","Liver","Muscle","Kidney")
thresholds <- c(2L, 3L, 4L, 5L)

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

for (tiss in tissues) {
  df_tiss <- flags_all %>%
    dplyr::filter(Tissue == tiss)
  
  if (nrow(df_tiss) == 0) {
    warning("No flags for tissue: ", tiss)
    next
  }
  
  pdf_file <- paste0("Euler_Adaptive_vs_ConvDE_", tiss, "_grid.pdf")
  pdf(pdf_file, width = 8, height = 8)
  par(mfrow = c(2, 2), mar = c(4, 4, 4, 2))
  
  for (thr in thresholds) {
    draw_euler_panel(df_tiss, tiss, thr)
  }
  
  dev.off()
  message("Wrote grid PDF: ", pdf_file)
}

############################################################
## 8) Save list of genes adaptive + convergentDry (by thresholds)
############################################################

genes_adaptive_conv <- combined_ogs_rep %>%
  dplyr::filter(!is.na(gene_name), gene_name != "") %>%
  dplyr::pull(gene_name) %>%
  unique()

write.table(genes_adaptive_conv,
            file = "Adaptive_and_ConvergentDry_gene_names.txt",
            quote = FALSE, row.names = FALSE, col.names = FALSE)

############################################################
## 9) GO on fixed list: Adaptive_and_ConvergentDry_gene_names.txt
############################################################

genes_list <- scan("Adaptive_and_ConvergentDry_gene_names.txt",
                   what = character(), quiet = TRUE)
genes_list <- unique(genes_list[genes_list != ""])

if (!exists("bg_universe_global") || is.null(bg_universe_global)) {
  bg_universe_global <- get_global_background_entrez(dds, og_rep_all, org.Hs.eg.db)
}

if (!is.null(bg_universe_global) && length(genes_list) >= 5) {
  gene_map_list <- clusterProfiler::bitr(
    genes_list,
    fromType = "SYMBOL",
    toType   = "ENTREZID",
    OrgDb    = org.Hs.eg.db
  )
  
  if (!is.null(gene_map_list) && nrow(gene_map_list) > 0) {
    ego_list <- enrichGO(
      gene          = gene_map_list$ENTREZID,
      universe      = bg_universe_global,
      OrgDb         = org.Hs.eg.db,
      keyType       = "ENTREZID",
      ont           = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.5,
      qvalueCutoff  = 0.2,
      readable      = TRUE
    )
    
    if (!is.null(ego_list) && nrow(as.data.frame(ego_list)) > 0) {
      p_list <- dotplot(ego_list, showCategory = 20) +
        ggtitle("GO BP – Adaptive & ConvergentDry gene list\n(global expressed-gene background)")
      print(p_list)
      
      write.csv(as.data.frame(ego_list),
                "GO_BP_Hs_Adaptive_and_ConvergentDry_LIST_enrichment.csv",
                row.names = FALSE)
      ggsave("GO_BP_Hs_Adaptive_and_ConvergentDry_LIST_dotplot.png",
             p_list, width = 8, height = 6, dpi = 300)
    } else {
      message("No significant GO terms for Adaptive_and_ConvergentDry gene list.")
    }
  } else {
    message("No symbols from Adaptive_and_ConvergentDry list could be mapped to ENTREZ.")
  }
} else {
  message("Background universe missing or too few genes in Adaptive_and_ConvergentDry list.")
}

############################################################
## 10) Grouped bar plot: EvoGeneX class proportions per tissue
############################################################

library(scales)  # for percent_format

evog_files <- list.files(pattern = "^EvoGeneX_.*_results\\.csv$")

evog_all <- dplyr::bind_rows(
  lapply(evog_files, function(f) {
    read.csv(f, stringsAsFactors = FALSE)
  })
)

evog_all$class <- factor(
  evog_all$class,
  levels = c("adaptive", "constrained", "neutral", "failed_fit")
)

class_summary <- evog_all %>%
  dplyr::group_by(Tissue, class) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop_last") %>%
  dplyr::mutate(
    total = sum(n),
    prop  = n / total
  ) %>%
  dplyr::ungroup()

p_classes <- ggplot(class_summary,
                    aes(x = Tissue, y = prop, fill = class)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    x = "Tissue",
    y = "Proportion of OGs",
    fill = "EvoGeneX class",
    title = "EvoGeneX model classes per tissue"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )

print(p_classes)
ggsave("EvoGeneX_class_proportions_by_tissue.png",
       p_classes, width = 7, height = 5, dpi = 300)

############################################################
## 11–16) Scatter plots & phylogeny/boxplots
## (kept mostly as in your original script; they depend only
##  on files and objects already created above/earlier scripts)
##
## If you want, we can split these plotting sections into a
## separate "04_EvoGeneX_plots.R" later.
############################################################

## ---- The rest of your plotting functions and loops (sections 11–16) ---
##      can stay as they are, since they already work off:
##       - DE_OG_selected_* CSVs
##       - EvoGeneX_*_results.csv
##       - dds, og_rep_all, og_rep
##       - direct kallisto OG count files
##
## You can keep them below this comment block unchanged,
## or move them into a separate script for plotting only.
