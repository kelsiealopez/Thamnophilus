
## 1. Nextflow for QC, read trimming, QC, pseudoalignment, QC, and transcript quantification 


#### a. Get stats on the RNA sequencing files 

```bash

seqkit stats -T -a *.fastq.gz > fastq_stats.tsv

```

#### b. Run RNAseq NextFlow

```bash

nf_core_rnaseq_cross_resume_more_time_custom_config_KALLISTO_all_to_thaDol_redo_sakLuc.sh

#!/bin/bash
#SBATCH -p shared,edwards
#SBATCH -c 16
#SBATCH -t 0-12:00
#SBATCH --mem=150000
#SBATCH --mail-type=END
#SBATCH -o nf_core_rnaseq_%x_%j.out
#SBATCH -e nf_core_rnaseq_%x_%j.err

# Usage examples:
#   sbatch --job-name thaDol nf_core_rnaseq_cross.sh thaDol thaDol 1
#   sbatch --job-name thaTor_to_thaRuf nf_core_rnaseq_cross.sh thaTor thaRuf 1
#
# Arguments:
#   1 = RNA_SPECIES    (where fastq TSV + reads are; e.g. thaTor)
#   2 = REF_SPECIES    (which genome/annotation to map to; e.g. thaRuf)
#   3 = RUN_LABEL      (optional run label, default = 1)

set -euo pipefail

RNA_SPECIES="$1"              # species the RNA comes from (TSV + FASTQs)
REF_SPECIES="${2:-$1}"        # species providing the reference genome/annotation
RUN_LABEL="${3:-1}"           # optional, defaults to "1"

BASE="/n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus"

############################
# 1) Directories and files
############################

# RNA-side: where the StringTie TSV and sample info live
stringtie_dir="${BASE}/annotation/03_stringtie/${RNA_SPECIES}"
tsv="${stringtie_dir}/fastq_${RNA_SPECIES}.tsv"

# Reference-side: which genome/annotation to map to
fasta="${BASE}/annotation/02_toga_galGal/${REF_SPECIES}/target/${REF_SPECIES}_forLASTZ.fa"
gtf="${BASE}/annotation/04_merging/${REF_SPECIES}/${REF_SPECIES}.toga_galGal_taeGut_rnaseq.anno.gtf"

# Work/output directory: encode both RNA and reference species so it’s obvious
workdir="${BASE}/nf_pipeline/${RNA_SPECIES}_to_${REF_SPECIES}"
out_dir="${workdir}/output_${RUN_LABEL}"

csv_fullpath="${stringtie_dir}/fastq_${RNA_SPECIES}_stranded_fullpath_correct.csv"

mkdir -p "${workdir}"

############################
# 2) Build CSV from TSV
############################

#cd "${stringtie_dir}"

#if [[ ! -f "${tsv}" ]]; then
#  echo "ERROR: TSV file not found: ${tsv}" >&2
#  exit 1
#fi

# TSV format:
# sample_id \t /path/to/R1.fastq.gz \t /path/to/R2.fastq.gz
#
# Build nf-core/rnaseq CSV with strandedness=auto using Python for safety
#python3 - << EOF
#import csv
#from pathlib import Path

#tsv_path = Path("${tsv}")
#csv_path = Path("${csv_fullpath}")

#with tsv_path.open() as fin, csv_path.open("w", newline="") as fout:
#    writer = csv.writer(fout)
#    writer.writerow(["sample", "fastq_1", "fastq_2", "strandedness"])
#    for line in fin:
#        line = line.strip()
#        if not line:
#            continue
#        parts = line.split("\\t")
#        if len(parts) < 3:
#            continue
#        sample, r1, r2 = parts[:3]
#        writer.writerow([sample, r1, r2, "auto"])
#EOF

############################
# 3) Run nf-core/rnaseq
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

#### c. Need to filter for the longest transcripts of the chicken orthologs 


```bash

# filter the thamnophilus doliatus annotation to loose set of orthologs, one to one with intact partially intact and uncertain loss (TOGA) for those that are present in at least 8 species out of 11 with genomes 
python3 filter_toga_chicken_at_least_8_species_o2o_I_PI_UL.py
```

```bash
# filter the kallisto output (which was mapped to the galGal->thaDol annotation) to those orthologs ^ 
python3 filter_kallisto_by_relaxed_TOGA_longest.py
```

## 2. Differential expression analysis with DESeq2

## 3 OU two-regime modeling with EVE model 

## 4. Gene Ontology of DE, EVE, DE + EVE Convergently Evolving Genes

## 5. PCA 

## 6. Weighted Gene Co-expression Network Analysis

## 7. Associations with environmental data?
