
## 1. Nextflow for QC, read trimming, QC, pseudoalignment, QC, and transcript quantification 

```bash

nf_core_rnaseq_cross_resume_more_time_custom_config_KALLISTO_all_to_thaDol_from_galGal_anno.sh

#!/bin/bash
#SBATCH -p shared,edwards
#SBATCH -c 16
#SBATCH -t 0-12:00
#SBATCH --mem=150000
#SBATCH --mail-type=END
#SBATCH -o nf_core_rnaseq_%x_%j.out
#SBATCH -e nf_core_rnaseq_%x_%j.err

# Arguments:
#   1 = RNA_SPECIES    (where fastq CSV + reads are; e.g. thaTor)
#   2 = RUN_LABEL      (optional run label, default = 1)

set -euo pipefail

RNA_SPECIES="$1"        # species the RNA comes from (TSV/CSV + FASTQs)
REF_SPECIES="thaDol"    # always map TO thaDol
RUN_LABEL="${2:-1}"     # optional second arg is run label if redoing an analysis

BASE="/n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus"

############################
# 1) Directories and files
############################

# this has the full paths to the raw fastq sequencing file. structured:
# sample,R1,R2,auto

stringtie_dir="${BASE}/annotation/03_stringtie/${RNA_SPECIES}"

# CSV path
csv_fullpath="${stringtie_dir}/fastq_${RNA_SPECIES}_stranded_fullpath.csv"


# Reference: thaDol genome/annotation
fasta="${BASE}/annotation/02_toga_galGal/thaDol/target/thaDol_forLASTZ.fa"
gtf="${BASE}/annotation/04_merging/thaDol/galGal_to_thaDol_query_annotation.gtf"
# If you want the relaxed GTF instead, use:
# gtf="${BASE}/annotation/04_merging/thaDol/galGal_to_thaDol_query_annotation.relaxedGenes.gtf"

# Work/output directory
workdir="${BASE}/nf_pipeline/galGal_to_thaDol_anno/${RNA_SPECIES}_to_${REF_SPECIES}"
out_dir="${workdir}/output_${RUN_LABEL}"

mkdir -p "${workdir}"

############################
# 2) Run nf-core/rnaseq
############################

cd "${workdir}"

nextflow run nf-core/rnaseq \
    --input "${csv_fullpath}" \
    --outdir "${out_dir}" \
    --gtf "${gtf}" \
    --fasta "${fasta}" \
    -profile singularity \
    -resume \
    --pseudo_aligner kallisto \
    --skip_alignment \
    -c /n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/nfcore_rnaseq_star_longtime.config

```

#### Need to filter for the longest transcripts of the chicken orthologs 

## 2. Differential expression analysis with DESeq2

## 3 OU two-regime modeling with EVE model 

## 4. Gene Ontology of DE, EVE, DE + EVE Convergently Evolving Genes

## 5. PCA 

## 6. Weighted Gene Co-expression Network Analysis

## 7. Associations with environmental data?
