############################################################
## EvoGeneX post-processing:
## Adaptive vs convergent DE + GO with tissue-specific background
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
## 0) Prereqs: og_rep_all, og_rep, dds, helper functions
############################################################

## If global OG->repGeneName mapping not in memory, load or build it
if (!exists("og_rep_all")) {
  if (file.exists("OGs_ALL_repGeneName.csv")) {
    og_rep_all <- read.csv("OGs_ALL_repGeneName.csv", stringsAsFactors = FALSE)
  } else {
    # Build from OG_gene_map_all9 + isoforms
    og_map_file <- "/n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/fastoma/out_folder_9_species_mamba/OG_gene_map_all9.tsv"
    og_map <- read.delim(og_map_file, header = TRUE, stringsAsFactors = FALSE)
    
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
    
    # clean_transcript_root() and choose_rep_gene() must already be defined
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
    
    write.csv(og_rep_all,
              "OGs_ALL_repGeneName.csv",
              row.names = FALSE)
  }
}

## Per-tissue OG->repGeneName (from DE script)
if (!exists("og_rep")) {
  og_rep <- read.csv("OGs_DE_in_multiple_pairs_ALL_tissues_repGeneName.csv",
                     stringsAsFactors = FALSE)
}

## EvoGeneX tissues vector
tissues <- c("Brain","Heart","Liver","Muscle","Kidney")

## Pairs: dry vs mesic
pairs_arid <- tribble(
  ~pair_id,             ~spA,     ~spB,     ~dry_species, ~mesic_species,
  "dysSti_vs_sakCri",   "dysSti", "sakCri", "sakCri",     "dysSti",
  "sakLuc_vs_sakCan",   "sakLuc", "sakCan", "sakCan",     "sakLuc",
  "thaAmb_vs_thaPel",   "thaAmb", "thaPel", "thaPel",     "thaAmb",
  "thaBer_vs_thaAtr",   "thaBer", "thaAtr", "thaBer",     "thaAtr",
  "thaRuf_vs_thaTor",   "thaRuf", "thaTor", "thaTor",     "thaRuf"
)

## dds is needed for tissue-specific background
if (!exists("dds")) {
  stop("dds (DESeqDataSet) is not in memory. 
       Source your DE script or load an RData containing dds before running this.")
}

## Helper: tissue-specific background based on expressed OGs in dds
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

## Helper: global background = all expressed OGs in dds
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

## Function: per-tissue direction wrt dry species (DE, not EvoGeneX)
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
    
    # UP in dry and adaptive
    up_adapt <- adaptive %>%
      inner_join(conv_up, by = "OG")
    
    # DOWN in dry and adaptive
    down_adapt <- adaptive %>%
      inner_join(conv_down, by = "OG")
    
    message("    UP in dry & adaptive:   ", nrow(up_adapt))
    message("    DOWN in dry & adaptive: ", nrow(down_adapt))
    
    ## Build simple per‑OG tables using 'counts' only
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
## 2) GO for adaptive + convergently DE in dry (all tissues combined)
##    with expressed-gene background (global universe)
############################################################

# Rebuild combined_ogs_rep if needed
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
## 3) Build per-tissue overlap summary (adaptive vs convergent)
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
## 4) Build per-tissue flag tables (adaptive & convergent)
############################################################

# If og_rep not in memory, load it again (safety)
if (!exists("og_rep")) {
  og_rep <- read.csv("OGs_DE_in_multiple_pairs_ALL_tissues_repGeneName.csv",
                     stringsAsFactors = FALSE)
}

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

length(genes_adaptive_conv)
head(genes_adaptive_conv)

write.table(genes_adaptive_conv,
            file = "Adaptive_and_ConvergentDry_gene_names.txt",
            quote = FALSE, row.names = FALSE, col.names = FALSE)




############################################################
## 9) GO on fixed list: Adaptive_and_ConvergentDry_gene_names.txt
############################################################

# Read gene symbols from file
genes_list <- scan("Adaptive_and_ConvergentDry_gene_names.txt",
                   what = character(), quiet = TRUE)
genes_list <- unique(genes_list[genes_list != ""])

# Ensure global background is available
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
      # dotplot
      p_list <- dotplot(ego_list, showCategory = 20) +
        ggtitle("GO BP – Adaptive & ConvergentDry gene list\n(global expressed-gene background)")
      print(p_list)
      
      # write table + plot
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
## 10) Grouped bar plot: proportion of EvoGeneX classes per tissue
############################################################

library(scales)  # for percent_format

# Collect all EvoGeneX result tables
evog_files <- list.files(pattern = "^EvoGeneX_.*_results\\.csv$")

evog_all <- dplyr::bind_rows(
  lapply(evog_files, function(f) {
    read.csv(f, stringsAsFactors = FALSE)
  })
)

# Ensure 'class' is a factor with desired order
evog_all$class <- factor(
  evog_all$class,
  levels = c("adaptive", "constrained", "neutral", "failed_fit")
)

# Count and compute proportions per tissue
class_summary <- evog_all %>%
  dplyr::group_by(Tissue, class) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop_last") %>%
  dplyr::mutate(
    total = sum(n),
    prop  = n / total
  ) %>%
  dplyr::ungroup()

# Grouped bar plot
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
## 11) Scatter: DE log2FC vs EvoGeneX signal for adaptive OGs
############################################################

plot_adaptive_scatter <- function(tiss,
                                  spA,
                                  spB,
                                  dry_label = NULL,
                                  y_col = "p_OU2_vs_BM") {
  # EvoGeneX results for this tissue
  evog_file <- paste0("EvoGeneX_", tiss, "_results.csv")
  if (!file.exists(evog_file)) {
    warning("EvoGeneX file not found for ", tiss, ": ", evog_file)
    return(invisible(NULL))
  }
  evog_res <- read.csv(evog_file, stringsAsFactors = FALSE)
  
  # Keep only adaptive OGs and required numeric column
  if (!y_col %in% colnames(evog_res)) {
    warning("Column ", y_col, " not found in EvoGeneX results for ", tiss)
    return(invisible(NULL))
  }
  
  evog_adapt <- evog_res %>%
    dplyr::filter(class == "adaptive") %>%
    dplyr::select(OG, Tissue, !!y_col)
  
  if (nrow(evog_adapt) == 0) {
    warning("No adaptive OGs for ", tiss)
    return(invisible(NULL))
  }
  
  # DE results for this pair/tissue (same naming as in your script)
  de_file <- paste0("DE_OG_selected_", spA, "_vs_", spB, "_", tiss, "_kallisto.csv")
  if (!file.exists(de_file)) {
    warning("DE file not found: ", de_file)
    return(invisible(NULL))
  }
  
  de_tab <- read.csv(de_file, row.names = 1, check.names = FALSE)
  de_tab$OG <- rownames(de_tab)
  
  # Standard DE filter (same as earlier in your script)
  de_sig <- de_tab %>%
    dplyr::filter(
      !is.na(padj),
      padj < 0.05,
      !is.na(log2FoldChange),
      abs(log2FoldChange) >= 1
    ) %>%
    dplyr::select(OG, log2FoldChange)
  
  if (nrow(de_sig) == 0) {
    warning("No significant DE OGs for ", spA, " vs ", spB, " in ", tiss)
    return(invisible(NULL))
  }
  
  # Join DE with adaptive EvoGeneX results
  merged <- de_sig %>%
    dplyr::inner_join(evog_adapt, by = "OG") %>%
    dplyr::filter(!is.na(.data[[y_col]]), .data[[y_col]] > 0)
  
  if (nrow(merged) == 0) {
    message("No overlap of DE & adaptive OGs for ", spA, "_vs_", spB, " in ", tiss)
    return(invisible(NULL))
  }
  
  # Define y as -log10 of the chosen p-value column
  merged <- merged %>%
    dplyr::mutate(y_value = -log10(.data[[y_col]]))
  
  pair_label <- if (!is.null(dry_label)) {
    paste0(spA, "_vs_", spB, " (dry = ", dry_label, ")")
  } else {
    paste0(spA, "_vs_", spB)
  }
  
  p_scatter <- ggplot(merged,
                      aes(x = log2FoldChange, y = y_value)) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
    geom_vline(xintercept = c(-1, 1), linetype = "dotted", color = "grey70") +
    geom_point(alpha = 0.7) +
    labs(
      x = "log2 fold change (DESeq2)",
      y = paste0("-log10(", y_col, ")"),
      title = paste0("Adaptive OGs – DE vs EvoGeneX signal\n",
                     "Tissue: ", tiss, ", Pair: ", pair_label)
    ) +
    theme_bw()
  
  print(p_scatter)
  
  out_name <- paste0(
    "Scatter_DE_vs_EvoGeneX_", tiss, "_", spA, "_vs_", spB, "_", y_col, ".png"
  )
  ggsave(out_name, p_scatter, width = 6, height = 5, dpi = 300)
  message("Wrote scatter plot: ", out_name)
  
  invisible(merged)
}

## Example: loop over all tissues and all pairs in 'pairs_arid'
for (tiss in tissues) {
  for (i in seq_len(nrow(pairs_arid))) {
    spA <- pairs_arid$spA[i]
    spB <- pairs_arid$spB[i]
    dry <- pairs_arid$dry_species[i]
    
    message("Scatter for ", tiss, " | ", spA, "_vs_", spB)
    plot_adaptive_scatter(
      tiss      = tiss,
      spA       = spA,
      spB       = spB,
      dry_label = dry,
      y_col     = "p_OU2_vs_BM"   # change here if you want a different EvoGeneX column
    )
  }
}







#
#
# show convergent, like one plot per tissue 




############################################################
## 12) Per-tissue scatter: DE log2FC vs EvoGeneX loglik_OU2
##     Colored by species pair, highlighting convergent genes
############################################################

plot_scatter_per_tissue <- function(tiss) {
  message("=== Scatter plot for tissue: ", tiss, " ===")
  
  # 1) EvoGeneX results for this tissue
  evog_file <- paste0("EvoGeneX_", tiss, "_results.csv")
  if (!file.exists(evog_file)) {
    warning("  EvoGeneX file not found: ", evog_file)
    return(invisible(NULL))
  }
  evog_res <- read.csv(evog_file, stringsAsFactors = FALSE)
  
  if (!"loglik_OU2" %in% colnames(evog_res)) {
    warning("  loglik_OU2 column not found in EvoGeneX results for ", tiss)
    return(invisible(NULL))
  }
  
  evog_sub <- evog_res %>%
    dplyr::select(OG, Tissue, loglik_OU2, class)
  
  # 2) Build DE table across all species pairs for this tissue
  de_list <- list()
  for (i in seq_len(nrow(pairs_arid))) {
    pid <- pairs_arid$pair_id[i]
    spA <- pairs_arid$spA[i]
    spB <- pairs_arid$spB[i]
    
    de_file <- paste0("DE_OG_selected_", spA, "_vs_", spB, "_", tiss, "_kallisto.csv")
    if (!file.exists(de_file)) {
      warning("  DE file not found: ", de_file)
      next
    }
    
    de_tab <- read.csv(de_file, row.names = 1, check.names = FALSE)
    de_tab$OG <- rownames(de_tab)
    
    de_sig <- de_tab %>%
      dplyr::filter(
        !is.na(padj),
        padj < 0.05,
        !is.na(log2FoldChange),
        abs(log2FoldChange) >= 1
      ) %>%
      dplyr::select(OG, log2FoldChange)
    
    if (nrow(de_sig) == 0) next
    
    de_sig$pair_id <- pid
    de_list[[pid]] <- de_sig
  }
  
  if (length(de_list) == 0) {
    warning("  No DE data for any pair in tissue ", tiss)
    return(invisible(NULL))
  }
  
  de_all <- dplyr::bind_rows(de_list)
  
  # 3) Join DE with EvoGeneX (loglik_OU2)
  merged_all <- de_all %>%
    dplyr::inner_join(evog_sub, by = "OG")
  
  if (nrow(merged_all) == 0) {
    warning("  No overlap of DE and EvoGeneX OGs for ", tiss)
    return(invisible(NULL))
  }
  
  # 4) Get convergent genes for this tissue (dry-direction DE across pairs)
  df_dir <- get_dry_direction(tiss)
  if (is.null(df_dir)) {
    warning("  No directional DE info for ", tiss, " to define convergence.")
    # We still can plot without highlighting
    merged_all$conv_flag <- "Not_convergent"
  } else {
    counts <- df_dir %>%
      dplyr::filter(directionDry %in% c("dry_up","dry_down")) %>%
      dplyr::group_by(OG) %>%
      dplyr::summarise(
        n_dry_up   = sum(directionDry == "dry_up"),
        n_dry_down = sum(directionDry == "dry_down"),
        .groups    = "drop"
      )
    
    # Threshold: ≥3 pairs normally, ≥2 pairs for Kidney
    thr_conv <- if (tiss == "Kidney") 2L else 3L
    
    conv_ogs_thr <- counts %>%
      dplyr::filter(n_dry_up >= thr_conv | n_dry_down >= thr_conv) %>%
      dplyr::pull(OG) %>%
      unique()
    
    merged_all <- merged_all %>%
      dplyr::mutate(
        conv_flag = ifelse(OG %in% conv_ogs_thr,
                           paste0("Convergent_≥", thr_conv, "pairs"),
                           "Not_convergent")
      )
  }
  
  merged_all$conv_flag <- factor(
    merged_all$conv_flag,
    levels = c("Not_convergent",
               paste0("Convergent_≥", if (tiss == "Kidney") 2L else 3L, "pairs"))
  )
  
  # 5) Plot:
  #    - all points colored by pair_id
  #    - convergent points overplotted as bigger, outlined points
  # 5) Plot:
  #    - all points colored by pair_id
  #    - convergent points overplotted with same color + black outline
  p_scatter <- ggplot(merged_all,
                      aes(x = log2FoldChange, y = loglik_OU2, color = pair_id)) +
    # base layer: all DE+EvoGeneX points
    geom_point(alpha = 0.5, size = 1.5) +
    # highlight convergent genes: same fill color as pair, black outline
    geom_point(
      data = merged_all %>% dplyr::filter(conv_flag != "Not_convergent"),
      aes(x = log2FoldChange, y = loglik_OU2, fill = pair_id),
      inherit.aes = FALSE,
      size  = 1.5,
      shape = 21,      # filled circle with outline
      color = "black", # outline color
      stroke = 0.6
    ) +
    scale_fill_discrete(guide = "none") +  # don't duplicate legend
    labs(
      x = "log2 fold change (DESeq2)",
      y = "EvoGeneX loglik_OU2",
      color = "Species pair",
      title = paste0(
        "Tissue: ", tiss,
        " – DE vs EvoGeneX (loglik_OU2)\n",
        "Highlight: convergent DE in dry (≥",
        if (tiss == "Kidney") 2L else 3L, " pairs)"
      )
    ) +
    theme_bw() +
    theme(
      legend.position = "right",
      plot.title = element_text(size = 11)
    )
  
  print(p_scatter)
  
  out_file <- paste0(
    "Scatter_log2FC_vs_loglikOU2_", tiss, "_allPairs_convergentHighl.png"
  )
  ggsave(out_file, p_scatter, width = 7, height = 5, dpi = 300)
  message("  Wrote scatter plot: ", out_file)
  
  invisible(merged_all)
}

## Run for all tissues
for (tiss in tissues) {
  plot_scatter_per_tissue(tiss)
}







# want to get top percentiles 
#
#
#
#


############################################################
## 12) Per-tissue scatter: DE log2FC vs EvoGeneX loglik_OU2
##     Colored by species pair, highlighting convergent genes
##     + top 15% (most negative loglik) thresholds and gene list
############################################################

plot_scatter_per_tissue <- function(tiss) {
  message("=== Scatter plot for tissue: ", tiss, " ===")
  
  # 1) EvoGeneX results for this tissue
  evog_file <- paste0("EvoGeneX_", tiss, "_results.csv")
  if (!file.exists(evog_file)) {
    warning("  EvoGeneX file not found: ", evog_file)
    return(invisible(NULL))
  }
  evog_res <- read.csv(evog_file, stringsAsFactors = FALSE)
  
  if (!"loglik_OU2" %in% colnames(evog_res)) {
    warning("  loglik_OU2 column not found in EvoGeneX results for ", tiss)
    return(invisible(NULL))
  }
  
  evog_sub <- evog_res %>%
    dplyr::select(OG, Tissue, loglik_OU2, class)
  
  # 2) Build DE table across all species pairs for this tissue
  de_list <- list()
  for (i in seq_len(nrow(pairs_arid))) {
    pid <- pairs_arid$pair_id[i]
    spA <- pairs_arid$spA[i]
    spB <- pairs_arid$spB[i]
    
    de_file <- paste0("DE_OG_selected_", spA, "_vs_", spB, "_", tiss, "_kallisto.csv")
    if (!file.exists(de_file)) {
      warning("  DE file not found: ", de_file)
      next
    }
    
    de_tab <- read.csv(de_file, row.names = 1, check.names = FALSE)
    de_tab$OG <- rownames(de_tab)
    
    de_sig <- de_tab %>%
      dplyr::filter(
        !is.na(padj),
        padj < 0.05,
        !is.na(log2FoldChange),
        abs(log2FoldChange) >= 1
      ) %>%
      dplyr::select(OG, log2FoldChange)
    
    if (nrow(de_sig) == 0) next
    
    de_sig$pair_id <- pid
    de_list[[pid]] <- de_sig
  }
  
  if (length(de_list) == 0) {
    warning("  No DE data for any pair in tissue ", tiss)
    return(invisible(NULL))
  }
  
  de_all <- dplyr::bind_rows(de_list)
  
  # 3) Join DE with EvoGeneX (loglik_OU2)
  merged_all <- de_all %>%
    dplyr::inner_join(evog_sub, by = "OG") %>%
    dplyr::filter(!is.na(loglik_OU2))
  
  if (nrow(merged_all) == 0) {
    warning("  No overlap of DE and EvoGeneX OGs for ", tiss)
    return(invisible(NULL))
  }
  
  # 4) Get convergent genes for this tissue (dry-direction DE across pairs)
  df_dir <- get_dry_direction(tiss)
  if (is.null(df_dir)) {
    warning("  No directional DE info for ", tiss, " to define convergence.")
    merged_all$conv_flag <- "Not_convergent"
  } else {
    counts <- df_dir %>%
      dplyr::filter(directionDry %in% c("dry_up","dry_down")) %>%
      dplyr::group_by(OG) %>%
      dplyr::summarise(
        n_dry_up   = sum(directionDry == "dry_up"),
        n_dry_down = sum(directionDry == "dry_down"),
        .groups    = "drop"
      )
    
    # Threshold: ≥3 pairs normally, ≥2 pairs for Kidney
    thr_conv <- if (tiss == "Kidney") 2L else 3L
    
    conv_ogs_thr <- counts %>%
      dplyr::filter(n_dry_up >= thr_conv | n_dry_down >= thr_conv) %>%
      dplyr::pull(OG) %>%
      unique()
    
    merged_all <- merged_all %>%
      dplyr::mutate(
        conv_flag = ifelse(OG %in% conv_ogs_thr,
                           paste0("Convergent_≥", thr_conv, "pairs"),
                           "Not_convergent")
      )
  }
  
  merged_all$conv_flag <- factor(
    merged_all$conv_flag,
    levels = c("Not_convergent",
               paste0("Convergent_≥", if (tiss == "Kidney") 2L else 3L, "pairs"))
  )
  
  # 5) Compute thresholds:
  #    - y_thr: 15th percentile (MOST NEGATIVE 15% have loglik_OU2 <= y_thr)
  #    - fc_thr: 85th percentile of |log2FC| (top 15% by magnitude)
  y_thr  <- stats::quantile(merged_all$loglik_OU2, 0.15, na.rm = TRUE)  # more negative = "better"
  fc_thr <- stats::quantile(abs(merged_all$log2FoldChange), 0.85, na.rm = TRUE)
  
  merged_all <- merged_all %>%
    dplyr::mutate(
      y_high  = loglik_OU2 <= y_thr,                     # most negative
      fc_high = abs(log2FoldChange) >= fc_thr
    )
  
  # 6) Plot with dashed lines
  p_scatter <- ggplot(merged_all,
                      aes(x = log2FoldChange, y = loglik_OU2, color = pair_id)) +
    geom_point(alpha = 0.5, size = 1.5) +
    geom_point(
      data = merged_all %>% dplyr::filter(conv_flag != "Not_convergent"),
      aes(x = log2FoldChange, y = loglik_OU2, fill = pair_id),
      inherit.aes = FALSE,
      size  = 1.5,
      shape = 21,
      color = "black",
      stroke = 0.6
    ) +
    geom_hline(yintercept = y_thr, linetype = "dashed", color = "grey40") +
    geom_vline(xintercept = c(-fc_thr, fc_thr), linetype = "dashed", color = "grey40") +
    scale_fill_discrete(guide = "none") +
    labs(
      x = "log2 fold change (DESeq2)",
      y = "EvoGeneX loglik_OU2",
      color = "Species pair",
      title = paste0(
        "Tissue: ", tiss,
        " – DE vs EvoGeneX (loglik_OU2)\n",
        "Top 15%: most negative loglik_OU2 & largest |log2FC|"
      )
    ) +
    theme_bw() +
    theme(
      legend.position = "right",
      plot.title = element_text(size = 11)
    )
  
  print(p_scatter)
  
  out_file <- paste0(
    "Scatter_log2FC_vs_loglikOU2_", tiss, "_allPairs_convergentHighl_top15pct.png"
  )
  ggsave(out_file, p_scatter, width = 7, height = 5, dpi = 300)
  message("  Wrote scatter plot: ", out_file)
  
  # 7) Classify genes by percentile category
  high_genes <- merged_all %>%
    dplyr::mutate(
      high_category = dplyr::case_when(
        y_high  & fc_high ~ "both_high",
        fc_high & !y_high ~ "fc_only_high",
        y_high  & !fc_high ~ "ll_only_high",
        TRUE ~ "none"
      )
    ) %>%
    dplyr::filter(high_category != "none") %>%
    dplyr::group_by(OG, Tissue, conv_flag, high_category) %>%
    dplyr::summarise(
      log2FC_values     = paste(round(log2FoldChange, 3), collapse = ";"),
      loglik_OU2_values = paste(round(loglik_OU2,       3), collapse = ";"),
      pair_ids          = paste(unique(pair_id),        collapse = ";"),
      .groups           = "drop"
    )
  
  if (nrow(high_genes) > 0) {
    if (exists("og_rep")) {
      high_genes <- high_genes %>%
        dplyr::left_join(
          og_rep %>%
            dplyr::select(Tissue, OG, gene_name, name_source),
          by = c("Tissue", "OG")
        )
    } else {
      high_genes$gene_name   <- NA_character_
      high_genes$name_source <- NA_character_
    }
    
    out_tab <- paste0("High_loglikOU2_log2FC_genes_", tiss, "_top15pct_withCategories.csv")
    write.csv(high_genes, out_tab, row.names = FALSE)
    message("  Wrote high-OF-interest gene table with categories: ", out_tab)
  } else {
    message("  No genes in top 15% of either |log2FC| or loglik_OU2 for ", tiss)
  }
  
  invisible(list(
    merged = merged_all,
    high_genes = high_genes
  ))
}

## Run for all tissues (15% version)
for (tiss in tissues) {
  plot_scatter_per_tissue(tiss)
}






############################################################
## 12b) Per-tissue scatter: DE log2FC vs EvoGeneX loglik_OU2
##      + top 5% (most negative loglik) thresholds and gene list
############################################################

plot_scatter_per_tissue <- function(tiss) {
  message("=== Scatter plot for tissue (top 5%): ", tiss, " ===")
  
  evog_file <- paste0("EvoGeneX_", tiss, "_results.csv")
  if (!file.exists(evog_file)) {
    warning("  EvoGeneX file not found: ", evog_file)
    return(invisible(NULL))
  }
  evog_res <- read.csv(evog_file, stringsAsFactors = FALSE)
  
  if (!"loglik_OU2" %in% colnames(evog_res)) {
    warning("  loglik_OU2 column not found in EvoGeneX results for ", tiss)
    return(invisible(NULL))
  }
  
  evog_sub <- evog_res %>%
    dplyr::select(OG, Tissue, loglik_OU2, class)
  
  de_list <- list()
  for (i in seq_len(nrow(pairs_arid))) {
    pid <- pairs_arid$pair_id[i]
    spA <- pairs_arid$spA[i]
    spB <- pairs_arid$spB[i]
    
    de_file <- paste0("DE_OG_selected_", spA, "_vs_", spB, "_", tiss, "_kallisto.csv")
    if (!file.exists(de_file)) {
      warning("  DE file not found: ", de_file)
      next
    }
    
    de_tab <- read.csv(de_file, row.names = 1, check.names = FALSE)
    de_tab$OG <- rownames(de_tab)
    
    de_sig <- de_tab %>%
      dplyr::filter(
        !is.na(padj),
        padj < 0.05,
        !is.na(log2FoldChange),
        abs(log2FoldChange) >= 1
      ) %>%
      dplyr::select(OG, log2FoldChange)
    
    if (nrow(de_sig) == 0) next
    
    de_sig$pair_id <- pid
    de_list[[pid]] <- de_sig
  }
  
  if (length(de_list) == 0) {
    warning("  No DE data for any pair in tissue ", tiss)
    return(invisible(NULL))
  }
  
  de_all <- dplyr::bind_rows(de_list)
  
  merged_all <- de_all %>%
    dplyr::inner_join(evog_sub, by = "OG") %>%
    dplyr::filter(!is.na(loglik_OU2))
  
  if (nrow(merged_all) == 0) {
    warning("  No overlap of DE and EvoGeneX OGs for ", tiss)
    return(invisible(NULL))
  }
  
  df_dir <- get_dry_direction(tiss)
  if (is.null(df_dir)) {
    warning("  No directional DE info for ", tiss, " to define convergence.")
    merged_all$conv_flag <- "Not_convergent"
  } else {
    counts <- df_dir %>%
      dplyr::filter(directionDry %in% c("dry_up","dry_down")) %>%
      dplyr::group_by(OG) %>%
      dplyr::summarise(
        n_dry_up   = sum(directionDry == "dry_up"),
        n_dry_down = sum(directionDry == "dry_down"),
        .groups    = "drop"
      )
    
    thr_conv <- if (tiss == "Kidney") 2L else 3L
    
    conv_ogs_thr <- counts %>%
      dplyr::filter(n_dry_up >= thr_conv | n_dry_down >= thr_conv) %>%
      dplyr::pull(OG) %>%
      unique()
    
    merged_all <- merged_all %>%
      dplyr::mutate(
        conv_flag = ifelse(OG %in% conv_ogs_thr,
                           paste0("Convergent_≥", thr_conv, "pairs"),
                           "Not_convergent")
      )
  }
  
  merged_all$conv_flag <- factor(
    merged_all$conv_flag,
    levels = c("Not_convergent",
               paste0("Convergent_≥", if (tiss == "Kidney") 2L else 3L, "pairs"))
  )
  
  # Top 5% = most negative 5% for loglik_OU2
  y_thr  <- stats::quantile(merged_all$loglik_OU2, 0.05, na.rm = TRUE)
  fc_thr <- stats::quantile(abs(merged_all$log2FoldChange), 0.95, na.rm = TRUE)
  
  merged_all <- merged_all %>%
    dplyr::mutate(
      y_high  = loglik_OU2 <= y_thr,
      fc_high = abs(log2FoldChange) >= fc_thr
    )
  
  p_scatter <- ggplot(merged_all,
                      aes(x = log2FoldChange, y = loglik_OU2, color = pair_id)) +
    geom_point(alpha = 0.5, size = 1.5) +
    geom_point(
      data = merged_all %>% dplyr::filter(conv_flag != "Not_convergent"),
      aes(x = log2FoldChange, y = loglik_OU2, fill = pair_id),
      inherit.aes = FALSE,
      size  = 1.5,
      shape = 21,
      color = "black",
      stroke = 0.6
    ) +
    geom_hline(yintercept = y_thr, linetype = "dashed", color = "grey40") +
    geom_vline(xintercept = c(-fc_thr, fc_thr), linetype = "dashed", color = "grey40") +
    scale_fill_discrete(guide = "none") +
    labs(
      x = "log2 fold change (DESeq2)",
      y = "EvoGeneX loglik_OU2",
      color = "Species pair",
      title = paste0(
        "Tissue: ", tiss,
        " – DE vs EvoGeneX (loglik_OU2)\n",
        "Top 5%: most negative loglik_OU2 & largest |log2FC|"
      )
    ) +
    theme_bw() +
    theme(
      legend.position = "right",
      plot.title = element_text(size = 11)
    )
  
  print(p_scatter)
  
  out_file <- paste0(
    "Scatter_log2FC_vs_loglikOU2_", tiss, "_allPairs_convergentHighl_top5pct.png"
  )
  ggsave(out_file, p_scatter, width = 7, height = 5, dpi = 300)
  message("  Wrote scatter plot: ", out_file)
  
  high_genes <- merged_all %>%
    dplyr::mutate(
      high_category = dplyr::case_when(
        y_high  & fc_high ~ "both_high",
        fc_high & !y_high ~ "fc_only_high",
        y_high  & !fc_high ~ "ll_only_high",
        TRUE ~ "none"
      )
    ) %>%
    dplyr::filter(high_category != "none") %>%
    dplyr::group_by(OG, Tissue, conv_flag, high_category) %>%
    dplyr::summarise(
      log2FC_values     = paste(round(log2FoldChange, 3), collapse = ";"),
      loglik_OU2_values = paste(round(loglik_OU2,       3), collapse = ";"),
      pair_ids          = paste(unique(pair_id),        collapse = ";"),
      .groups           = "drop"
    )
  
  if (nrow(high_genes) > 0) {
    if (exists("og_rep")) {
      high_genes <- high_genes %>%
        dplyr::left_join(
          og_rep %>%
            dplyr::select(Tissue, OG, gene_name, name_source),
          by = c("Tissue", "OG")
        )
    } else {
      high_genes$gene_name   <- NA_character_
      high_genes$name_source <- NA_character_
    }
    
    out_tab <- paste0("High_loglikOU2_log2FC_genes_", tiss, "_top5pct_withCategories.csv")
    write.csv(high_genes, out_tab, row.names = FALSE)
    message("  Wrote high-OF-interest gene table with categories: ", out_tab)
  } else {
    message("  No genes in top 5% of either |log2FC| or loglik_OU2 for ", tiss)
  }
  
  invisible(list(
    merged = merged_all,
    high_genes = high_genes
  ))
}

## Run for all tissues (5% version)
for (tiss in tissues) {
  plot_scatter_per_tissue(tiss)
}
############################################################
## Fill missing gene_name/name_source in *_top5pct_withCategories.csv
## using global OG->gene mapping (og_rep_all)
############################################################

# Safety: load og_rep_all if needed
if (!exists("og_rep_all")) {
  if (file.exists("OGs_ALL_repGeneName.csv")) {
    og_rep_all <- read.csv("OGs_ALL_repGeneName.csv", stringsAsFactors = FALSE)
  } else {
    stop("OGs_ALL_repGeneName.csv not found and og_rep_all not in memory.")
  }
}

# Find all the high-percentile files you just created
high_files <- list.files(
  pattern = "^High_loglikOU2_log2FC_genes_.*_top5pct_withCategories\\.csv$"
)

for (f in high_files) {
  message("Updating gene names in: ", f)
  
  df <- read.csv(f, stringsAsFactors = FALSE)
  
  # Join global mapping by OG
  df2 <- df %>%
    dplyr::left_join(
      og_rep_all %>%
        dplyr::select(OG, gene_name_global = gene_name,
                      name_source_global = name_source),
      by = "OG"
    ) %>%
    dplyr::mutate(
      gene_name   = ifelse(is.na(gene_name)   | gene_name   == "",
                           gene_name_global,   gene_name),
      name_source = ifelse(is.na(name_source) | name_source == "",
                           name_source_global, name_source)
    ) %>%
    dplyr::select(-gene_name_global, -name_source_global)
  
  # Overwrite or write to new file; here we overwrite
  write.csv(df2, f, row.names = FALSE)
}







############################################################
## 13) Phylogeny + per-species CPM boxplots
##     for OGs that are BOTH-high (top 5% in |log2FC| & loglik_OU2)
##     One PNG per OG, title includes OG ID + gene name
############################################################

library(ggtree)
library(ape)
library(phytools)

## Base directory where kallisto OG-level count files live
nf_base <- "/n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/nf_pipeline/thaDol_mapping_test"

############################################################
## 13.1 Read OG-level kallisto counts and compute CPM
############################################################

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

og_counts_list <- list()
for (run in names(files)) {
  path <- file.path(nf_base, files[[run]])
  if (!file.exists(path)) {
    warning("Missing OG count file for run: ", run, " at ", path)
    next
  }
  cat("Reading:", path, "for run", run, "\n")
  df <- read.delim(path, header = TRUE, check.names = FALSE)
  stopifnot(colnames(df)[1] == "OG")
  og_counts_list[[run]] <- df
}
stopifnot(length(og_counts_list) > 0)

merged <- Reduce(function(x, y) merge(x, y, by = "OG", all = FALSE), og_counts_list)
cat("Merged matrix has", nrow(merged), "OGs and", ncol(merged) - 1, "samples\n")

rownames(merged) <- merged$OG
merged$OG <- NULL
cts <- as.matrix(merged)  # OG x sample raw counts

## Build coldata from sample names
samples    <- colnames(cts)
tissue_smp <- sub(".*_", "", samples)           # after last "_"
id_species <- sub("_[^_]*$", "", samples)       # before last "_"
species    <- sub("[0-9]+.*$", "", id_species)  # leading letters

coldata <- data.frame(
  Sample  = samples,
  Species = species,
  Tissue  = tissue_smp,
  stringsAsFactors = FALSE
)

## Fix thaCap -> thaDol to match EvoGeneX / DE scripts
coldata$Species[coldata$Species == "thaCap"] <- "thaDol"

stopifnot(all(coldata$Sample == colnames(cts)))

## Compute CPM
lib_sizes <- colSums(cts)
cpm_mat   <- t( t(cts) / lib_sizes * 1e6 )  # OG x sample

############################################################
## 13.2 Tree path helper (same as in phylogeny script)
############################################################

get_tree_path_for_tissue <- function(tiss) {
  if (tiss == "Kidney") {
    "/n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/fastoma/species_tree_KIDNEY_13_sp.nwk"
  } else if (tiss == "Muscle") {
    "/n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/fastoma/species_tree_MUSCLE_13_sp.nwk"
  } else {
    # Brain, Heart, Liver default 14-species tree
    "/n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/fastoma/species_tree_14_sp.nwk"
  }
}

############################################################
## 13.3 Helper: tidy expression for a gene & tissue (CPM)
############################################################

############################################################
## SAFE get_expr_long: handle missing OG rows / samples
############################################################

get_expr_long <- function(og_id, tiss, expr_mat, coldata) {
  # 1) Check OG present in expression matrix
  if (!(og_id %in% rownames(expr_mat))) {
    warning("OG ", og_id, " not found in expression matrix; skipping.")
    return(NULL)
  }
  
  # 2) Samples for this tissue
  idx <- which(coldata$Tissue == tiss)
  if (length(idx) == 0) {
    warning("No samples for tissue ", tiss, " in coldata; skipping OG ", og_id)
    return(NULL)
  }
  
  samp_all <- coldata$Sample[idx]
  # Only keep samples that are actually in expr_mat
  samp <- intersect(samp_all, colnames(expr_mat))
  if (length(samp) == 0) {
    warning("No matching samples in expression matrix for tissue ", tiss,
            " (OG ", og_id, "); skipping.")
    return(NULL)
  }
  
  sp <- coldata$Species[match(samp, coldata$Sample)]
  exprv <- expr_mat[og_id, samp, drop = TRUE]  # CPM
  
  df <- data.frame(
    Sample  = samp,
    Species = sp,
    CPM     = as.numeric(exprv),
    stringsAsFactors = FALSE
  )
  
  df
}

############################################################
## 13.4 Plot function: phylogeny + per-species CPM boxplots
##      Title includes OG ID and gene_name (if available)
############################################################

plot_gene_on_tree_box <- function(tiss, og_id,
                                  expr_mat, coldata,
                                  gene_name_label = NULL,
                                  dry_species_global = c("sakCri","sakCan","thaBer","thaDol","thaPel","thaTor"),
                                  to_pdf = FALSE,
                                  outdir = "plots_adaptive_phylo_box_CPM_top5pct_both") {
  if (to_pdf && !dir.exists(outdir)) {
    dir.create(outdir, showWarnings = FALSE)
  }
  
  # Read tissue-specific tree
  tree_path <- get_tree_path_for_tissue(tiss)
  if (!file.exists(tree_path)) {
    stop("Tree file does not exist for tissue ", tiss, ": ", tree_path)
  }
  tree <- read.tree(tree_path)
  
  # Expression for this OG & tissue (CPM by sample)
  df <- get_expr_long(og_id, tiss, expr_mat, coldata)
  if (is.null(df) || nrow(df) == 0) {
    warning("No expression for OG ", og_id, " in tissue ", tiss)
    return(invisible(NULL))
  }
  
  # Restrict to species present in tree
  df <- df %>% dplyr::filter(Species %in% tree$tip.label)
  if (nrow(df) == 0) {
    warning("No overlapping species between tree and expression for OG ",
            og_id, " in tissue ", tiss)
    return(invisible(NULL))
  }
  
  # Prune tree to those species
  keep_tips <- intersect(tree$tip.label, unique(df$Species))
  tree <- keep.tip(tree, keep_tips)
  
  # Order factor for species to match tree tip order
  sp_order <- tree$tip.label
  df$Species <- factor(df$Species, levels = sp_order)
  
  # Define dry vs mesic for this tree
  dry_sp_in_tree   <- intersect(dry_species_global, sp_order)
  mesic_sp_in_tree <- setdiff(sp_order, dry_sp_in_tree)
  
  # Dry
  col_dry_main   <- "#E74504"
  col_dry_border <- col_dry_main
  col_dry_point  <- col_dry_main
  col_dry_box    <- adjustcolor(col_dry_main, alpha.f = 0.45)
  
  # Mesic
  col_mes_main   <- "#585757"
  col_mes_border <- col_mes_main
  col_mes_point  <- col_mes_main
  col_mes_box    <- adjustcolor(col_mes_main, alpha.f = 0.45)
  
  # CPM range for x-axis
  max_cpm <- max(df$CPM, na.rm = TRUE)
  x_max   <- max(1, max_cpm * 1.05)
  
  # Open PDF if requested
  if (to_pdf) {
    pdf_file <- file.path(outdir,
                          paste0("Phylo_box_CPM_", tiss, "_", og_id, ".pdf"))
    pdf(pdf_file, width = 8, height = 6)
  }
  
  # One outer title for both panels
  par(oma = c(1, 0, 1, 0))
  par(xpd = NA)
  
  ############################
  ## Panel 1: phylogeny     ##
  ############################
  par(fig = c(0.0, 0.48, 0.2475, 0.79), new = FALSE)
  par(mar = c(2, 1, 2, 0))
  
  plotTree(tree,
           direction = "rightwards",
           ftype = "i",
           fsize = 0.7)
  
  ############################
  ## Panel 2: CPM boxplots  ##
  ############################
  par(fig = c(0.48, 1.0, 0.1, 0.9), new = TRUE)
  par(mar = c(3, 0, 2, 1))
  
  n_sp <- length(sp_order)
  plot(NA,
       xlim = c(0, x_max),
       ylim = c(0.5, n_sp + 0.5),
       xlab = "CPM",
       ylab = "",
       yaxt = "n",
       main = "")
  
  # Ticks only, no labels
  axis(2,
       at = seq_len(n_sp),
       labels = FALSE,
       tck = -0.015,
       las = 2,
       cex.axis = 0.7)
  
  for (i in seq_along(sp_order)) {
    sp <- sp_order[i]
    vals <- df$CPM[df$Species == sp]
    if (length(vals) == 0) next
    
    bp_stats <- boxplot.stats(vals)$stats
    
    is_dry    <- sp %in% dry_sp_in_tree
    box_col    <- if (is_dry) col_dry_box    else col_mes_box
    border_col <- if (is_dry) col_dry_border else col_mes_border
    
    rect(xleft  = bp_stats[2],
         ybottom= i - 0.25,
         xright = bp_stats[4],
         ytop   = i + 0.25,
         col    = box_col,
         border = border_col,
         lwd    = 1.5)
    
    segments(x0 = bp_stats[3], y0 = i - 0.25,
             x1 = bp_stats[3], y1 = i + 0.25,
             col = border_col, lwd = 1.5)
    
    segments(x0 = bp_stats[1], y0 = i,
             x1 = bp_stats[2], y1 = i,
             col = border_col, lwd = 1)
    segments(x0 = bp_stats[4], y0 = i,
             x1 = bp_stats[5], y1 = i,
             col = border_col, lwd = 1)
    
    ## points: no jitter, clear fill, black outline, small
    y_jit <- rep(i, length(vals))
    points(vals, y_jit,
           pch = 21,
           bg  = NA,
           col = "black",
           lwd = 0.9,
           cex = 1.1)
  }
  
  legend("topright",
         legend = c("Dry", "Mesic"),
         pch = 21,
         pt.bg = "white",
         col = c(col_dry_point, col_mes_point),
         pt.cex = 1.4,
         bty = "n",
         title = "")
  
  # Title with OG and gene name if available
  title_txt <- paste0("Tissue: ", tiss, " | OG: ", og_id)
  if (!is.null(gene_name_label) && !is.na(gene_name_label) && gene_name_label != "") {
    title_txt <- paste0(title_txt, " (", gene_name_label, ")")
  }
  mtext(title_txt,
        outer = TRUE, side = 3, line = 0.1, cex = 1)
  
  if (to_pdf) {
    dev.off()
    message("Wrote PDF for ", tiss, " | OG ", og_id, " to ", pdf_file)
  }
  
  invisible(NULL)
}

############################################################
## 13.5 Plot ONLY OGs that are BOTH-high (top 5% both axes)
##      using the *_top5pct_withCategories.csv files
############################################################

tissues_to_plot <- c("Brain", "Heart", "Liver", "Muscle", "Kidney")

# Directory for PNGs of the BOTH-high genes
outdir_png_both <- "plots_adaptive_phylo_box_CPM_top5pct_both"
if (!dir.exists(outdir_png_both)) dir.create(outdir_png_both, showWarnings = FALSE)

# Optional: combined table of both_high genes across tissues
both_high_all <- list()

for (tiss in tissues_to_plot) {
  high_file <- paste0("High_loglikOU2_log2FC_genes_", tiss, "_top5pct_withCategories.csv")
  
  if (!file.exists(high_file)) {
    message("No high-percentile file for tissue ", tiss, ": ", high_file)
    next
  }
  
  message("\n=== Tissue: ", tiss, " ===")
  high_df <- read.csv(high_file, stringsAsFactors = FALSE)
  
  # Keep only OGs that are both_high (top 5% in both metrics)
  high_both <- high_df %>%
    dplyr::filter(high_category == "both_high")
  
  if (nrow(high_both) == 0) {
    message("  No both_high genes for ", tiss)
    next
  }
  
  message("  Found ", nrow(high_both), " OG entries with high_category == 'both_high'")
  
  # For plotting: unique OG x gene_name for this tissue
  og_list <- high_both %>%
    dplyr::select(OG, gene_name) %>%
    dplyr::distinct()
  
  # Store for combined list
  both_high_all[[tiss]] <- high_both %>%
    dplyr::mutate(Tissue = tiss)
  
  # Loop over OGs and make PNGs
  for (k in seq_len(nrow(og_list))) {
    og_id <- og_list$OG[k]
    gname <- og_list$gene_name[k]
    
    png_file <- file.path(outdir_png_both,
                          paste0("Phylo_box_CPM_", tiss, "_", og_id, ".png"))
    
    png(filename = png_file,
        width   = 418,
        height  = 482,
        units   = "px",
        res     = NA)
    
    plot_gene_on_tree_box(
      tiss            = tiss,
      og_id           = og_id,
      expr_mat        = cpm_mat,
      coldata         = coldata,
      gene_name_label = gname,
      to_pdf          = FALSE
    )
    
    dev.off()
    message("  Wrote PNG: ", png_file,
            "  (gene: ", ifelse(is.na(gname) || gname == "", "NA", gname), ")")
  }
}

# Combined BOTH-high list across tissues
if (length(both_high_all) > 0) {
  both_high_all_df <- dplyr::bind_rows(both_high_all)
  write.csv(both_high_all_df,
            "High_loglikOU2_log2FC_genes_ALL_tissues_top5pct_bothHigh.csv",
            row.names = FALSE)
  message("\nWrote combined both_high table: High_loglikOU2_log2FC_genes_ALL_tissues_top5pct_bothHigh.csv")
}






############################################################
## 14) Sanity check: log10(normalized counts+1) boxplots
##     for BOTH-high OGs (top 5% in |log2FC| & loglik_OU2)
##     One PNG per OG × tissue, species on x-axis
############################################################

library(ggplot2)

## 14.1 Build log10(normalized counts+1) from dds

if (!exists("dds")) {
  stop("dds object not found; cannot build normalized-count boxplots.")
}

norm_counts <- counts(dds, normalized = TRUE)
log_norm    <- log10(norm_counts + 1)

## Make sure we have a clean coldata from dds
coldata_dds <- as.data.frame(colData(dds))
coldata_dds$Sample  <- rownames(coldata_dds)
# Expect columns 'Species' and 'Tissue' already present in dds
stopifnot(all(colnames(norm_counts) == coldata_dds$Sample))

############################################################
## SAFE helper: extract expression for OG × tissue
############################################################

get_expr_long_dds <- function(og_id, tiss, expr_mat, coldata) {
  if (!(og_id %in% rownames(expr_mat))) {
    warning("OG ", og_id, " not found in normalized expression matrix; skipping.")
    return(NULL)
  }
  
  idx <- which(coldata$Tissue == tiss)
  if (length(idx) == 0) {
    warning("No samples for tissue ", tiss, " in coldata; skipping OG ", og_id)
    return(NULL)
  }
  
  samp_all <- coldata$Sample[idx]
  samp <- intersect(samp_all, colnames(expr_mat))
  if (length(samp) == 0) {
    warning("No matching samples in normalized matrix for tissue ", tiss,
            " (OG ", og_id, "); skipping.")
    return(NULL)
  }
  
  sp <- coldata$Species[match(samp, coldata$Sample)]
  exprv <- expr_mat[og_id, samp, drop = TRUE]
  
  data.frame(
    Sample  = samp,
    Species = sp,
    Tissue  = tiss,
    logNorm = as.numeric(exprv),
    stringsAsFactors = FALSE
  )
}

############################################################
## 14.2 Plot function: per-tissue boxplot for one OG
############################################################

plot_logNorm_box <- function(tiss, og_id,
                             expr_mat = log_norm,
                             coldata  = coldata_dds,
                             gene_name_label = NULL,
                             outdir = "plots_logNorm_box_top5pct_both") {
  if (!dir.exists(outdir)) {
    dir.create(outdir, showWarnings = FALSE)
  }
  
  df <- get_expr_long_dds(og_id, tiss, expr_mat, coldata)
  if (is.null(df) || nrow(df) == 0) {
    warning("No normalized expression for OG ", og_id, " in tissue ", tiss)
    return(invisible(NULL))
  }
  
  title_txt <- paste0("Tissue: ", tiss, " | OG: ", og_id)
  if (!is.null(gene_name_label) && !is.na(gene_name_label) && gene_name_label != "") {
    title_txt <- paste0(title_txt, " (", gene_name_label, ")")
  }
  
  p <- ggplot(df, aes(x = Species, y = logNorm, fill = Species)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.6) +
    geom_jitter(width = 0.15, size = 1.5, alpha = 0.8, color = "black") +
    labs(
      x = "Species",
      y = "log10(normalized counts + 1)",
      title = title_txt
    ) +
    theme_bw() +
    theme(
      axis.text.x  = element_text(angle = 45, hjust = 1),
      legend.position = "none",
      plot.title  = element_text(size = 10)
    )
  
  print(p)
  
  png_file <- file.path(outdir,
                        paste0("logNorm_box_", tiss, "_", og_id, ".png"))
  ggsave(png_file, p, width = 5, height = 4, dpi = 300)
  message("  Wrote logNorm boxplot: ", png_file)
  
  invisible(df)
}

############################################################
## 14.3 Loop over BOTH-high OGs (top 5%) and plot
############################################################

tissues_to_plot <- c("Brain", "Heart", "Liver", "Muscle", "Kidney")
outdir_logNorm  <- "plots_logNorm_box_top5pct_both"
if (!dir.exists(outdir_logNorm)) dir.create(outdir_logNorm, showWarnings = FALSE)

both_high_all <- list()

for (tiss in tissues_to_plot) {
  high_file <- paste0("High_loglikOU2_log2FC_genes_", tiss, "_top5pct_withCategories.csv")
  
  if (!file.exists(high_file)) {
    message("No high-percentile file for tissue ", tiss, ": ", high_file)
    next
  }
  
  message("\n=== Tissue: ", tiss, " – logNorm boxplots for BOTH-high OGs ===")
  high_df <- read.csv(high_file, stringsAsFactors = FALSE)
  
  high_both <- high_df %>%
    dplyr::filter(high_category == "both_high")
  
  if (nrow(high_both) == 0) {
    message("  No both_high genes for ", tiss)
    next
  }
  
  message("  Found ", nrow(high_both), " OG entries with high_category == 'both_high'")
  
  # For plotting: unique OG x gene_name for this tissue
  og_list <- high_both %>%
    dplyr::select(OG, gene_name) %>%
    dplyr::distinct()
  
  both_high_all[[tiss]] <- high_both %>%
    dplyr::mutate(Tissue = tiss)
  
  for (k in seq_len(nrow(og_list))) {
    og_id <- og_list$OG[k]
    gname <- og_list$gene_name[k]
    plot_logNorm_box(
      tiss            = tiss,
      og_id           = og_id,
      expr_mat        = log_norm,
      coldata         = coldata_dds,
      gene_name_label = gname,
      outdir          = outdir_logNorm
    )
  }
}

# Optional: combined table of BOTH-high OGs across tissues (already exists for phylo,
# but here we just re-create for convenience)
if (length(both_high_all) > 0) {
  both_high_all_df <- dplyr::bind_rows(both_high_all)
  write.csv(both_high_all_df,
            "High_loglikOU2_log2FC_genes_ALL_tissues_top5pct_bothHigh_forLogNorm.csv",
            row.names = FALSE)
  message("\nWrote combined both_high table (for logNorm plots): ",
          "High_loglikOU2_log2FC_genes_ALL_tissues_top5pct_bothHigh_forLogNorm.csv")
}





############################################################
## 15) Phylogeny + per-species log10(normalized+1) boxplots
##     for BOTH-high OGs (top 5% |log2FC| & most-negative loglik_OU2)
##     Same layout as CPM plots, but using DESeq2-normalized counts
############################################################

# Requires:
#  - log_norm (log10(normalized counts + 1)) from section 14
#  - coldata_dds (colData(dds) with Sample, Species, Tissue)
#  - get_tree_path_for_tissue()
#  - High_loglikOU2_log2FC_genes_*_top5pct_withCategories.csv

if (!exists("log_norm") || !exists("coldata_dds")) {
  stop("log_norm or coldata_dds not found; run section 14 first.")
}

############################################################
## 15.1 SAFE helper: extract logNorm expression for OG × tissue
############################################################

get_expr_long_logNorm <- function(og_id, tiss, expr_mat, coldata) {
  if (!(og_id %in% rownames(expr_mat))) {
    warning("OG ", og_id, " not found in log_norm matrix; skipping.")
    return(NULL)
  }
  
  idx <- which(coldata$Tissue == tiss)
  if (length(idx) == 0) {
    warning("No samples for tissue ", tiss, " in coldata; skipping OG ", og_id)
    return(NULL)
  }
  
  samp_all <- coldata$Sample[idx]
  samp <- intersect(samp_all, colnames(expr_mat))
  if (length(samp) == 0) {
    warning("No matching samples in log_norm matrix for tissue ", tiss,
            " (OG ", og_id, "); skipping.")
    return(NULL)
  }
  
  sp <- coldata$Species[match(samp, coldata$Sample)]
  exprv <- expr_mat[og_id, samp, drop = TRUE]
  
  data.frame(
    Sample  = samp,
    Species = sp,
    Tissue  = tiss,
    logNorm = as.numeric(exprv),
    stringsAsFactors = FALSE
  )
}

############################################################
## 15.2 Phylogeny + logNorm boxplots (same layout as CPM)
############################################################


############################################################
## 15.2 Phylogeny + logNorm boxplots (fixed margins / device)
############################################################

plot_gene_on_tree_box_logNorm <- function(tiss, og_id,
                                          expr_mat = log_norm,
                                          coldata = coldata_dds,
                                          gene_name_label = NULL,
                                          dry_species_global = c("sakCri","sakCan","thaBer","thaDol","thaPel","thaTor"),
                                          outdir = "plots_logNorm_phylo_box_top5pct_both") {
  if (!dir.exists(outdir)) {
    dir.create(outdir, showWarnings = FALSE)
  }
  
  tree_path <- get_tree_path_for_tissue(tiss)
  if (!file.exists(tree_path)) {
    warning("Tree file does not exist for tissue ", tiss, ": ", tree_path)
    return(invisible(NULL))
  }
  tree <- read.tree(tree_path)
  
  df <- get_expr_long_logNorm(og_id, tiss, expr_mat, coldata)
  if (is.null(df) || nrow(df) == 0) {
    warning("No logNorm expression for OG ", og_id, " in tissue ", tiss)
    return(invisible(NULL))
  }
  
  df <- df %>% dplyr::filter(Species %in% tree$tip.label)
  if (nrow(df) == 0) {
    warning("No overlapping species between tree and logNorm expression for OG ",
            og_id, " in tissue ", tiss)
    return(invisible(NULL))
  }
  
  keep_tips <- intersect(tree$tip.label, unique(df$Species))
  tree <- keep.tip(tree, keep_tips)
  
  sp_order <- tree$tip.label
  df$Species <- factor(df$Species, levels = sp_order)
  
  dry_sp_in_tree   <- intersect(dry_species_global, sp_order)
  mesic_sp_in_tree <- setdiff(sp_order, dry_sp_in_tree)
  
  col_dry_main   <- "#E74504"
  col_dry_border <- col_dry_main
  col_dry_box    <- adjustcolor(col_dry_main, alpha.f = 0.45)
  
  col_mes_main   <- "#585757"
  col_mes_border <- col_mes_main
  col_mes_box    <- adjustcolor(col_mes_main, alpha.f = 0.45)
  
  max_val <- max(df$logNorm, na.rm = TRUE)
  x_max   <- max(0.1, max_val * 1.05)
  
  png_file <- file.path(outdir,
                        paste0("Phylo_logNorm_box_", tiss, "_", og_id, ".png"))
  png(filename = png_file,
      width   = 418,
      height  = 482,
      units   = "px",
      res     = NA)  # use raw pixels, avoids margin issues
  
  # reset basic par state
  par(mfrow = c(1, 1))
  par(oma = c(1, 0, 1, 0))
  par(xpd = NA)
  
  ############################
  ## Panel 1: phylogeny     ##
  ############################
  par(fig = c(0.0, 0.48, 0.2, 0.8), new = FALSE)
  par(mar = c(2, 1, 2, 0))
  
  plotTree(tree,
           direction = "rightwards",
           ftype = "i",
           fsize = 0.7)
  
  ############################
  ## Panel 2: logNorm boxplots ##
  ############################
  par(fig = c(0.48, 1.0, 0.1, 0.9), new = TRUE)
  par(mar = c(3, 0, 2, 1))
  
  n_sp <- length(sp_order)
  plot(NA,
       xlim = c(0, x_max),
       ylim = c(0.5, n_sp + 0.5),
       xlab = "log10(normalized counts + 1)",
       ylab = "",
       yaxt = "n",
       main = "")
  
  axis(2,
       at = seq_len(n_sp),
       labels = FALSE,
       tck = -0.015,
       las = 2,
       cex.axis = 0.7)
  
  for (i in seq_along(sp_order)) {
    sp <- sp_order[i]
    vals <- df$logNorm[df$Species == sp]
    if (length(vals) == 0) next
    
    bp_stats <- boxplot.stats(vals)$stats
    
    is_dry    <- sp %in% dry_sp_in_tree
    box_col    <- if (is_dry) col_dry_box    else col_mes_box
    border_col <- if (is_dry) col_dry_border else col_mes_border
    
    rect(xleft  = bp_stats[2],
         ybottom= i - 0.25,
         xright = bp_stats[4],
         ytop   = i + 0.25,
         col    = box_col,
         border = border_col,
         lwd    = 1.5)
    
    segments(x0 = bp_stats[3], y0 = i - 0.25,
             x1 = bp_stats[3], y1 = i + 0.25,
             col = border_col, lwd = 1.5)
    
    segments(x0 = bp_stats[1], y0 = i,
             x1 = bp_stats[2], y1 = i,
             col = border_col, lwd = 1)
    segments(x0 = bp_stats[4], y0 = i,
             x1 = bp_stats[5], y1 = i,
             col = border_col, lwd = 1)
    
    y_jit <- rep(i, length(vals))
    points(vals, y_jit,
           pch = 21,
           bg  = NA,
           col = "black",
           lwd = 0.9,
           cex = 1.1)
  }
  
  legend("topright",
         legend = c("Dry", "Mesic"),
         pch = 21,
         pt.bg = "white",
         col = c(col_dry_border, col_mes_border),
         pt.cex = 1.4,
         bty = "n",
         title = "")
  
  title_txt <- paste0("Tissue: ", tiss, " | OG: ", og_id)
  if (!is.null(gene_name_label) && !is.na(gene_name_label) && gene_name_label != "") {
    title_txt <- paste0(title_txt, " (", gene_name_label, ")")
  }
  mtext(title_txt,
        outer = TRUE, side = 3, line = 0.1, cex = 1)
  
  dev.off()
  message("  Wrote phylo+logNorm PNG: ", png_file)
  
  invisible(NULL)
}


############################################################
## 15.3 Loop over BOTH-high OGs (top 5%) and plot phylo+logNorm
############################################################

tissues_to_plot <- c("Brain", "Heart", "Liver", "Muscle", "Kidney")
outdir_phylo_logNorm <- "plots_logNorm_phylo_box_top5pct_both"
if (!dir.exists(outdir_phylo_logNorm)) dir.create(outdir_phylo_logNorm, showWarnings = FALSE)

for (tiss in tissues_to_plot) {
  high_file <- paste0("High_loglikOU2_log2FC_genes_", tiss, "_top5pct_withCategories.csv")
  
  if (!file.exists(high_file)) {
    message("No high-percentile file for tissue ", tiss, ": ", high_file)
    next
  }
  
  message("\n=== Tissue: ", tiss, " – phylo+logNorm for BOTH-high OGs ===")
  high_df <- read.csv(high_file, stringsAsFactors = FALSE)
  
  high_both <- high_df %>%
    dplyr::filter(high_category == "both_high")
  
  if (nrow(high_both) == 0) {
    message("  No both_high genes for ", tiss)
    next
  }
  
  message("  Found ", nrow(high_both), " OG entries with high_category == 'both_high'")
  
  og_list <- high_both %>%
    dplyr::select(OG, gene_name) %>%
    dplyr::distinct()
  
  for (k in seq_len(nrow(og_list))) {
    og_id <- og_list$OG[k]
    gname <- og_list$gene_name[k]
    plot_gene_on_tree_box_logNorm(
      tiss            = tiss,
      og_id           = og_id,
      expr_mat        = log_norm,
      coldata         = coldata_dds,
      gene_name_label = gname,
      outdir          = outdir_phylo_logNorm
    )
  }
}





############################################################
## 16) Phylogeny + boxplots for ALL EvoGeneX-adaptive OGs
##     - CPM (orange/gray, as before)
##     - log10(normalized counts + 1)
############################################################

# Directories for adaptive OG plots
outdir_cpm_adapt     <- "plots_adaptive_phylo_box_CPM_allAdaptive"
outdir_logNorm_adapt <- "plots_logNorm_phylo_box_allAdaptive"

if (!dir.exists(outdir_cpm_adapt)) dir.create(outdir_cpm_adapt, showWarnings = FALSE)
if (!dir.exists(outdir_logNorm_adapt)) dir.create(outdir_logNorm_adapt, showWarnings = FALSE)

tissues_adapt <- c("Brain", "Heart", "Liver", "Muscle", "Kidney")

for (tiss in tissues_adapt) {
  evog_file <- paste0("EvoGeneX_", tiss, "_results.csv")
  if (!file.exists(evog_file)) {
    message("No EvoGeneX results for ", tiss, ": ", evog_file)
    next
  }
  
  evog_res <- read.csv(evog_file, stringsAsFactors = FALSE)
  adaptive_ogs <- evog_res %>%
    dplyr::filter(class == "adaptive", Tissue == tiss) %>%
    dplyr::select(OG) %>%
    dplyr::distinct()
  
  if (nrow(adaptive_ogs) == 0) {
    message("No adaptive OGs for tissue ", tiss)
    next
  }
  
  message("\n=== Tissue: ", tiss, " – phylo plots for ALL adaptive OGs (",
          nrow(adaptive_ogs), " OGs) ===")
  
  # Add gene_name (if available) from og_rep for titles
  og_names <- adaptive_ogs %>%
    dplyr::left_join(
      og_rep %>% dplyr::filter(Tissue == tiss) %>%
        dplyr::select(OG, gene_name),
      by = "OG"
    )
  
  for (k in seq_len(nrow(og_names))) {
    og_id <- og_names$OG[k]
    gname <- og_names$gene_name[k]
    
    ## CPM phylogeny+boxplot
    png_cpm <- file.path(outdir_cpm_adapt,
                         paste0("Phylo_CPM_adaptive_", tiss, "_", og_id, ".png"))
    png(filename = png_cpm,
        width   = 418,
        height  = 482,
        units   = "px",
        res     = NA)
    plot_gene_on_tree_box(
      tiss            = tiss,
      og_id           = og_id,
      expr_mat        = cpm_mat,
      coldata         = coldata,
      gene_name_label = gname,
      to_pdf          = FALSE
    )
    dev.off()
    message("  Wrote CPM phylo plot: ", png_cpm,
            "  (gene: ", ifelse(is.na(gname) || gname == "", "NA", gname), ")")
    
    ## logNorm phylogeny+boxplot
    png_ln <- file.path(outdir_logNorm_adapt,
                        paste0("Phylo_logNorm_adaptive_", tiss, "_", og_id, ".png"))
    plot_gene_on_tree_box_logNorm(
      tiss            = tiss,
      og_id           = og_id,
      expr_mat        = log_norm,
      coldata         = coldata_dds,
      gene_name_label = gname,
      outdir          = outdir_logNorm_adapt
    )
    # plot_gene_on_tree_box_logNorm opens/closes its own PNG device
  }
}
