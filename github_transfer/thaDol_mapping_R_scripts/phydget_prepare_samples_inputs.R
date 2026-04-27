############################################################
# Build PhyDGET --samples lines for any tissue
############################################################

setwd("/n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/nf_pipeline/thaDol_mapping_test")

make_phydget_sample_block <- function(tiss) {
  fn <- paste0("Phydget_counts_", tiss, "_OG_all9.tsv")
  if (!file.exists(fn)) {
    warning("Counts file not found for tissue ", tiss, ": ", fn)
    return(NULL)
  }
  
  # 1) Read just the header
  hdr <- read.delim(fn, nrows = 0, check.names = FALSE)
  samples <- colnames(hdr)[-1]  # drop "OG"
  
  # 2) Parse species and tissue from sample names
  tissue  <- sub(".*_", "", samples)
  ids     <- sub("_[^_]*$", "", samples)
  species <- sub("[0-9]+.*$", "", ids)
  
  coldata <- data.frame(
    Sample  = samples,
    Species = species,
    Tissue  = tissue,
    stringsAsFactors = FALSE
  )
  
  # Fix thaCap -> thaDol
  coldata$Species[coldata$Species == "thaCap"] <- "thaDol"
  
  cat("Tissue:", tiss, "\n")
  print(table(coldata$Tissue))
  
  # 3) Build --sample lines by species
  species_order <- sort(unique(coldata$Species))
  
  sample_lines <- sapply(species_order, function(sp) {
    sp_samples <- coldata$Sample[coldata$Species == sp]
    if (length(sp_samples) == 0) return(NA_character_)
    paste0("--sample ", sp, ":", paste(sp_samples, collapse = ","))
  })
  
  sample_lines <- sample_lines[!is.na(sample_lines)]
  
  # 4) Write to file
  out_file <- paste0("Phydget_", tiss, "_samples_block.txt")
  writeLines(
    c(paste0("# Auto-generated --sample lines for ", tiss),
      sample_lines),
    con = out_file
  )
  message("Wrote sample block to ", out_file)
  
  invisible(sample_lines)
}

## Call for each tissue
tissues_to_run <- c("Kidney", "Muscle", "Brain", "Heart", "Liver")
lapply(tissues_to_run, make_phydget_sample_block)