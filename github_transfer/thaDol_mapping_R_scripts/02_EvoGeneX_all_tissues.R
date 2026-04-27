############################################################
## 02_EvoGeneX_all_tissues.R
##
## EvoGeneX on Thamnophilus OG counts
## BM vs OU (global) vs OU2 (dry vs mesic), all tissues
## Uses tissue-specific trees for Kidney and Muscle
## Re-uses OG counts & coldata saved by 01_DE_and_GO.R
############################################################

setwd("/n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/nf_pipeline/thaDol_mapping_test")

library(ape)
library(phangorn)
library(EvoGeneX)
library(dplyr)

############################################################
## 0) Load OG counts and coldata from 01_DE_and_GO.R
############################################################

# These were saved at the end of the PCA section in 01_DE_and_GO.R
cts      <- readRDS("obj_cts_OG_counts.rds")       # OG x sample (collapsed)
coldata  <- readRDS("obj_coldata_samples.rds")     # Sample, Species, Tissue

# Safety check: columns vs Sample names
stopifnot(all(colnames(cts) == coldata$Sample))

## =========================================
## EXCLUDE SAMPLES HERE (same as DE analysis)
## =========================================
samples_to_drop <- c("thaPal01_Liver")

keep_cols <- setdiff(colnames(cts), samples_to_drop)
cts <- cts[, keep_cols, drop = FALSE]
coldata <- coldata[coldata$Sample %in% keep_cols, , drop = FALSE]
stopifnot(all(colnames(cts) == coldata$Sample))
## =========================================

## log2 expression
expr_log2 <- log2(cts + 1)

############################################################
## 1) Main EvoGeneX loop over tissues
############################################################

tissues_to_run <- c("Brain", "Liver", "Heart", "Muscle", "Kidney")
ogs_to_run     <- rownames(expr_log2)

all_results <- list()

for (tiss in tissues_to_run) {
  cat("=== Running EvoGeneX for tissue:", tiss, "===\n")
  
  ## 1A) Choose tissue-specific tree -----------------------
  tree_path <- switch(
    tiss,
    "Kidney" = "/n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/fastoma/species_tree_KIDNEY_13_sp.nwk",
    "Muscle" = "/n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/fastoma/species_tree_MUSCLE_13_sp.nwk",
    # default for Brain, Heart, Liver:
    "/n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/fastoma/species_tree_14_sp.nwk"
  )
  stopifnot(file.exists(tree_path))
  
  phylo_tree <- read.tree(tree_path)
  
  ## keep only species that actually have expression in this tissue
  sp_in_expr_tiss <- unique(coldata$Species[coldata$Tissue == tiss])
  phylo_tree <- keep.tip(phylo_tree,
                         base::intersect(phylo_tree$tip.label, sp_in_expr_tiss))
  
  cat("Tree tips for", tiss, "after matching expression:\n")
  print(phylo_tree$tip.label)
  
  Ntip  <- length(phylo_tree$tip.label)
  Nnode <- phylo_tree$Nnode
  inner_nodes <- (Ntip + 1):(Ntip + Nnode)
  tip_names   <- phylo_tree$tip.label
  
  ## 1B) Build regimes (global and dry/mesic) for this tissue ----
  mc <- ape::mrca(phylo_tree)  # [tips x tips] node IDs
  
  internal_pairs <- lapply(inner_nodes, function(nd) {
    idx <- which(mc == nd, arr.ind = TRUE)
    if (nrow(idx) == 0) stop("No tip pair found for node ", nd)
    i <- idx[1, 1]; j <- idx[1, 2]
    c(tip_names[i], tip_names[j])
  })
  
  internal_pairs_df <- data.frame(
    node  = sapply(internal_pairs, `[`, 1),
    node2 = sapply(internal_pairs, `[`, 2),
    stringsAsFactors = FALSE
  )
  
  stopifnot(nrow(internal_pairs_df) == Nnode)
  stopifnot(nrow(unique(internal_pairs_df)) == Nnode)
  
  # Global regime for this tissue
  leaf_df <- data.frame(
    node   = tip_names,
    node2  = "",
    regime = "global",
    stringsAsFactors = FALSE
  )
  
  internal_df <- internal_pairs_df %>% mutate(regime = "global")
  regime_global <- bind_rows(leaf_df, internal_df)
  stopifnot(nrow(regime_global) == Ntip + Nnode)
  
  reg_global_file <- paste0("regime_global_thamnophilus_", tiss, ".csv")
  write.csv(regime_global, reg_global_file, row.names = FALSE)
  
  # Dry vs mesic regime for this tissue
  dry_species <- c("sakCri","sakCan","thaBer","thaDol","thaPel","thaTor")
  dry_species <- base::intersect(dry_species, tip_names)
  mesic_species <- base::setdiff(tip_names, dry_species)
  
  cat("Dry species for", tiss, ":\n");   print(dry_species)
  cat("Mesic species for", tiss, ":\n"); print(mesic_species)
  
  two_leaf <- data.frame(
    node   = tip_names,
    node2  = "",
    regime = ifelse(tip_names %in% dry_species, "dry", "mesic"),
    stringsAsFactors = FALSE
  )
  
  assign_internal_regime <- function(tr, node, dry_vec, mesic_vec) {
    tips_idx <- phangorn::Descendants(tr, node, type = "tips")[[1]]
    tips     <- tr$tip.label[tips_idx]
    if (all(tips %in% dry_vec))      "dry"
    else if (all(tips %in% mesic_vec)) "mesic"
    else                             "mesic"
  }
  
  internal_regime <- sapply(
    inner_nodes,
    assign_internal_regime,
    tr = phylo_tree,
    dry_vec = dry_species,
    mesic_vec = mesic_species
  )
  
  two_internal <- internal_pairs_df %>%
    mutate(regime = internal_regime)
  
  regime_drymesic <- bind_rows(two_leaf, two_internal)
  stopifnot(nrow(regime_drymesic) == Ntip + Nnode)
  
  reg_drymesic_file <- paste0("regime_dry_mesic_thamnophilus_", tiss, ".csv")
  write.csv(regime_drymesic, reg_drymesic_file, row.names = FALSE)
  
  ##########################################################
  ## 2) Expression helper + EvoGeneX model objects
  ##########################################################
  
  build_expr_tall_tiss <- function(og_id) {
    idx <- which(coldata$Tissue == tiss &
                   coldata$Species %in% phylo_tree$tip.label)
    if (length(idx) == 0) return(NULL)
    
    samp  <- coldata$Sample[idx]
    sp    <- coldata$Species[idx]
    exprv <- expr_log2[og_id, samp]
    
    df <- data.frame(
      species   = sp,
      replicate = samp,
      exprval   = as.numeric(exprv),
      stringsAsFactors = FALSE
    )
    
    tab <- table(df$species)
    min_reps <- min(tab)
    if (min_reps < 2) return(NULL)
    
    df_bal <- df %>%
      dplyr::group_by(species) %>%
      dplyr::slice_head(n = min_reps) %>%
      dplyr::ungroup()
    
    df_bal
  }
  
  evog  <- EvoGeneX()
  evog$setTree(tree_path)
  evog$setRegimes(reg_global_file)
  
  brown <- Brown()
  brown$setTree(tree_path)
  
  evog2 <- EvoGeneX()
  evog2$setTree(tree_path)
  evog2$setRegimes(reg_drymesic_file)
  
  # Degrees of freedom
  ou_dof    <- 1 + 1 + 1 + 1   # theta, alpha, sigma^2, gamma
  brown_dof <- 1 + 1 + 1       # theta, sigma^2, gamma
  ou2_dof   <- 1 + 1 + 2 + 1   # theta_dry, theta_mesic, alpha, sigma^2, gamma
  
  process_one_gene_tiss <- function(og_id) {
    expr_tall <- build_expr_tall_tiss(og_id)
    if (is.null(expr_tall)) {
      return(data.frame(
        OG           = og_id,
        Tissue       = tiss,
        loglik_BM    = NA,
        loglik_OU    = NA,
        loglik_OU2   = NA,
        p_OU_vs_BM   = NA,
        p_OU2_vs_BM  = NA,
        p_OU2_vs_OU  = NA,
        class        = "insufficient_reps",
        stringsAsFactors = FALSE
      ))
    }
    
    ou_res <- evog$fit(
      expr_tall,
      format        = "tall",
      species_col   = "species",
      replicate_col = "replicate",
      exprval_col   = "exprval",
      alpha         = 0.1,
      gamma         = 0.01
    )
    
    brown_res <- brown$fit(
      expr_tall,
      format        = "tall",
      species_col   = "species",
      replicate_col = "replicate",
      exprval_col   = "exprval",
      gamma         = 0.01
    )
    
    ou2_res <- evog2$fit(
      expr_tall,
      format        = "tall",
      species_col   = "species",
      replicate_col = "replicate",
      exprval_col   = "exprval",
      alpha         = 0.1,
      gamma         = 0.01
    )
    
    p_ou_vs_bm  <- 1 - pchisq((ou_res$loglik  - brown_res$loglik) * 2,
                              df = ou_dof  - brown_dof)
    p_ou2_vs_bm <- 1 - pchisq((ou2_res$loglik - brown_res$loglik) * 2,
                              df = ou2_dof - brown_dof)
    p_ou2_vs_ou <- 1 - pchisq((ou2_res$loglik - ou_res$loglik)    * 2,
                              df = ou2_dof - ou_dof)
    
    gene_class <- dplyr::case_when(
      is.na(p_ou_vs_bm) ~ "failed_fit",
      p_ou_vs_bm >= 0.05 ~ "neutral",
      p_ou_vs_bm < 0.05 & (p_ou2_vs_bm >= 0.05 | p_ou2_vs_ou >= 0.05) ~ "constrained",
      p_ou2_vs_bm < 0.05 & p_ou2_vs_ou < 0.05 ~ "adaptive",
      TRUE ~ "uncertain"
    )
    
    data.frame(
      OG          = og_id,
      Tissue      = tiss,
      loglik_BM   = brown_res$loglik,
      loglik_OU   = ou_res$loglik,
      loglik_OU2  = ou2_res$loglik,
      p_OU_vs_BM  = p_ou_vs_bm,
      p_OU2_vs_BM = p_ou2_vs_bm,
      p_OU2_vs_OU = p_ou2_vs_ou,
      class       = gene_class,
      stringsAsFactors = FALSE
    )
  }
  
  ##########################################################
  ## 3) Run EvoGeneX for all OGs in THIS tissue
  ##########################################################
  
  res_list <- lapply(ogs_to_run, process_one_gene_tiss)
  res_df <- bind_rows(res_list)
  
  out_csv <- paste0("EvoGeneX_", tiss, "_results.csv")
  write.csv(res_df, out_csv, row.names = FALSE)
  cat("Wrote EvoGeneX results for", tiss, "to", out_csv, "\n")
  
  all_results[[tiss]] <- res_df
}

############################################################
## 2) Save combined EvoGeneX results object
############################################################

saveRDS(all_results, file = "obj_EvoGeneX_all_tissues_results.rds")

############################################################
## End of 02_EvoGeneX_all_tissues.R
############################################################
