# Thamnophilus


# 1. Assembly

# 2. Assembly QC check

# 3. Repeat Masking 

```bash

cd /n/netscratch/edwards_lab/Lab/kelsielopez/suboscines/repeats

thaDolRepMod_RepMasking.sh

#!/bin/bash
#SBATCH -t 3-00:00
#SBATCH -p shared,edwards
#SBATCH -c 32
#SBATCH --mem=100G
#SBATCH -o thaDolRepMod_RepMasking_%j.out
#SBATCH -e thaDolRepMod_RepMasking_%j.err
#SBATCH --mail-type=END

set -euxo pipefail

# -- Paths & variables
GENOME="/n/netscratch/edwards_lab/Lab/kelsielopez/suboscines/Doliatus_pacbio/pacbio_hifi_assembly/workflow/results/assembly/ThaDol_18-293.p_ctg.fa"
PREFIX="thaDol"
SRF_LIB="/n/netscratch/edwards_lab/Lab/kelsielopez/suboscines/concat_fastq/suboscine_srf.fa"
BIRD_LIB="/n/netscratch/edwards_lab/Lab/kelsielopez/repeats/rep_masker/Passerellidae.final_TElibrary.fa"
N_THREADS=32
OUTDIR="/n/netscratch/edwards_lab/Lab/kelsielopez/suboscines/repeats/thaDol"

mkdir -p "$OUTDIR"
cd "$OUTDIR"

# -- 1. Build RepeatModeler Database --
/n/home03/kelsielopez/repeat-annotation/RepeatModeler-2.0.3/BuildDatabase -name ${PREFIX}_RepeatModelerDB -engine ncbi "$GENOME" 2>&1 | tee build_db.log

# -- 2. Run RepeatModeler --
/n/home03/kelsielopez/repeat-annotation/RepeatModeler-2.0.3/RepeatModeler -pa $N_THREADS -engine ncbi -database ${PREFIX}_RepeatModelerDB 2>&1 | tee repeatmodeler.log

# -- 3. Combine RepeatModeler output with SRF and Bird libraries --
RAW_LIB="${PREFIX}_RepeatModelerDB-families.fa"
COMBINED_LIB="${PREFIX}_combined_lib.fa"
cat "$RAW_LIB" "$SRF_LIB" "$BIRD_LIB" > "$COMBINED_LIB"

# -- 4. Run RepeatMasker with the combined custom library
RepeatMasker -pa $N_THREADS -e ncbi -xsmall -gff -a -lib "$COMBINED_LIB" -dir "$OUTDIR/mask_out" "$GENOME" 2>&1 | tee repeatmasker.log

echo "RepeatModeler and RepeatMasker run complete!"

# -a flag means output alignments. so i can look at k2p divergence. 


```
```bash

# actually i am skipping srf. it is kind of ovekill because im not focusing on repeats in this analysis 
(env_nf_rna) [kelsielopez@holy8a24402 01_repMod_repMask]$ head repMod_repMask_PacBio_genomes.sh -n 500
#!/bin/bash
#SBATCH -t 3-00:00
#SBATCH -p shared,edwards
#SBATCH -c 32
#SBATCH --mem=100G
#SBATCH -o repMod_repMask_PacBio_genomes_%A_%a.out
#SBATCH -e repMod_repMask_PacBio_genomes_%A_%a.err
#SBATCH --mail-type=END
#SBATCH --array=0-2

set -euxo pipefail

# -- Variables
GENOMES_LIST="/n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/annotation/01_repMod_repMask/PacBio_genomes.txt"
BIRD_LIB="/n/netscratch/edwards_lab/Lab/kelsielopez/repeats/rep_masker/Passerellidae.final_TElibrary.fa"
N_THREADS=32

# -- Get the genome file path for this task
GENOME=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" $GENOMES_LIST)
BASENAME=$(basename "$GENOME" .p_ctg.fa)  # e.g., "S_canadensis"

OUTDIR="/n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/annotation/01_repMod_repMask/$BASENAME"
mkdir -p "$OUTDIR"
cd "$OUTDIR"

# -- 1. Build RepeatModeler Database --
/n/home03/kelsielopez/repeat-annotation/RepeatModeler-2.0.3/BuildDatabase -name ${BASENAME}_RepeatModelerDB -engine ncbi "$GENOME" 2>&1 | tee build_db.log

# -- 2. Run RepeatModeler --
/n/home03/kelsielopez/repeat-annotation/RepeatModeler-2.0.3/RepeatModeler -pa $N_THREADS -engine ncbi -database ${BASENAME}_RepeatModelerDB 2>&1 | tee repeatmodeler.log

# -- 3. Combine RepeatModeler output with Bird library ONLY
RAW_LIB="${BASENAME}_RepeatModelerDB-families.fa"
COMBINED_LIB="${BASENAME}_combined_lib.fa"
cat "$RAW_LIB" "$BIRD_LIB" > "$COMBINED_LIB"

# -- 4. Run RepeatMasker with the combined custom library
RepeatMasker -pa $N_THREADS -e ncbi -xsmall -gff -a -lib "$COMBINED_LIB" -dir "$OUTDIR/mask_out" "$GENOME" 2>&1 | tee repeatmasker.log

echo "RepeatModeler and RepeatMasker run complete for $BASENAME!"




```




# 4. RNA-based annotation - stringtie

```bash

# prepare reference for STAR


01_thaDol_prepareGenomeIndex.sh

#!/bin/bash
#SBATCH -t 0-12:00
#SBATCH -p test
#SBATCH -c 32
#SBATCH --mem=100G
#SBATCH -o 01_thaDol_prepareGenomeIndex_%j.out
#SBATCH -e 01_thaDol_prepareGenomeIndex_%j.err
#SBATCH --mail-type=END
# Set paths
GENOME=/n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/Doliatus_pacbio/01_repeat_masking/03_known_out/ThaDol_18-293.p_ctg.known_mask.masked.fasta
GENOME_INDEX=thaDol_genome_index
mkdir -p $GENOME_INDEX

# Create STAR genome index
STAR --runThreadN 32 --runMode genomeGenerate --genomeDir $GENOME_INDEX --genomeFastaFiles $GENOME

```

```bash

# prepare a file with name, forward, and reverse path

ython_env1) [kelsielopez@holylogin06 annotation]$ ls *tsv
thaDolFastqPaths.tsv
(python_env1) [kelsielopez@holylogin06 annotation]$ head thaDolFastqPaths.tsv
thaCap01_Liver	/n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/RNA_round1_re_demultiplex/fastq/P504_KLope19066_B08v1_22-328_Liver_S58_L004_R1_001.fastq.gz	/n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/RNA_round1_re_demultiplex/fastq/P504_KLope19066_B08v1_22-328_Liver_S58_L004_R2_001.fastq.gz
thaCap01_Heart	/n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/RNA_round1_re_demultiplex/fastq/P504_KLope19066_C08v1_22-328_Heart_S59_L004_R1_001.fastq.gz	/n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/RNA_round1_re_demultiplex/fastq/P504_KLope19066_C08v1_22-328_Heart_S59_L004_R2_001.fastq.gz
thaCap01_Kidney	/n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/RNA_round1_re_demultiplex/fastq/P504_KLope19066_D08v1_22-328_Kidney_S60_L004_R1_001.fastq.gz	/n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/RNA_round1_re_demultiplex/fastq/P504_KLope19066_D08v1_22-328_Kidney_S60_L004_R2_001.fastq.gz
thaCap01_Brain	/n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/RNA_round1_re_demultiplex/fastq/P504_KLope19066_E08v1_22-328_Brain_S61_L004_R1_001.fastq.gz	/n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/RNA_round1_re_demultiplex/fastq/P504_KLope19066_E08v1_22-328_Brain_S61_L004_R2_001.fastq.gz
thaCap02_Brain	/n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/RNA_round1_re_demultiplex/fastq/P504_KLope19066_A09v1_22-330_Brain_S65_L004_R1_001.fastq.gz	/n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/RNA_round1_re_demultiplex/fastq/P504_KLope19066_A09v1_22-330_Brain_S65_L004_R2_001.fastq.gz
thaCap02_Liver	/n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/RNA_round1_re_demultiplex/fastq/P504_KLope19066_F08v1_22-330_Liver_S62_L004_R1_001.fastq.gz	/n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/RNA_round1_re_demultiplex/fastq/P504_KLope19066_F08v1_22-330_Liver_S62_L004_R2_001.fastq.gz
thaCap02_Heart	/n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/RNA_round1_re_demultiplex/fastq/P504_KLope19066_G08v1_22-330_Heart_S63_L004_R1_001.fastq.gz	/n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/RNA_round1_re_demultiplex/fastq/P504_KLope19066_G08v1_22-330_Heart_S63_L004_R2_001.fastq.gz
thaCap02_Kidney	/n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/RNA_round1_re_demultiplex/fastq/P504_KLope19066_H08v1_22-330_Kidney_S64_L004_R1_001.fastq.gz	/n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/RNA_round1_re_demultiplex/fastq/P504_KLope19066_H08v1_22-330_Kidney_S64_L004_R2_001.fastq.gz
thaCap03_Liver	/n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/RNA_round1_re_demultiplex/fastq/P504_KLope19066_B09v1_22-334_Liver_S66_L004_R1_001.fastq.gz	/n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/RNA_round1_re_demultiplex/fastq/P504_KLope19066_B09v1_22-334_Liver_S66_L004_R2_001.fastq.gz
thaCap03_Heart	/n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/RNA_round1_re_demultiplex/fastq/P504_KLope19066_C09v1_22-334_Heart_S67_L004_R1_001.fastq.gz	/n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/RNA_round1_re_demultiplex/fastq/P504_KLope19066_C09v1_22-334_Heart_S67_L004_R2_001.fastq.gz
(python_env1) [kelsielopez@holylogin06 annotation]$ wc -l thaDolFastqPaths.tsv
11 thaDolFastqPaths.tsv
(python_env1) [kelsielopez@holylogin06 annotation]$ 


```

```bash
# now batch to align everything to the index 

cd /n/netscratch/edwards_lab/Lab/kelsielopez/suboscines/annotation

02_trim_map_thaDol_to_masked_fasta_STAR.sh

#!/bin/bash
#SBATCH --job-name=02_trim_map_thaDol_to_masked_fasta_STAR
#SBATCH -p edwards,shared
#SBATCH --array=1-12           # set XX=number of samples/lines in samples.tsv
#SBATCH --cpus-per-task=8
#SBATCH --mem=100G
#SBATCH --time=24:00:00
#SBATCH -o logs/array_%A_%a.out
#SBATCH -e logs/array_%A_%a.err
#SBATCH --mail-type=END,FAIL

set -euxo pipefail


# 2. Get sample info from file
SAMPLEFILE=thaDolFastqPaths.tsv

IFS=$'\t' read SAMPLE R1 R2 <<< $(sed -n "${SLURM_ARRAY_TASK_ID}p" $SAMPLEFILE)

echo "Processing $SAMPLE: $R1 $R2"

# 3. Parameters/Paths
GENOME=/n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/Doliatus_pacbio/01_repeat_masking/03_known_out/ThaDol_18-293.p_ctg.known_mask.masked.fasta
REFERENCE=$GENOME
GENOME_INDEX=/n/netscratch/edwards_lab/Lab/kelsielopez/suboscines/annotation/thaDol_genome_index
THREADS=8
TRIMDIR=trimmed
MAPDIR=mapped
METRICSDIR=metrics
OUTPUT_METRICS_FILE=${METRICSDIR}/${SAMPLE}_alignment_summary_metrics.txt


mkdir -p $TRIMDIR $MAPDIR $METRICSDIR

# 4. Trim_Galore
trim_galore --paired -j $THREADS --retain_unpaired --phred33 --output_dir $TRIMDIR \
  --length 35 -q 0 --stringency 5 $R1 $R2

TR1=${TRIMDIR}/$(basename $R1 .fastq.gz)_val_1.fq.gz
TR2=${TRIMDIR}/$(basename $R2 .fastq.gz)_val_2.fq.gz

# 5. STAR mapping
STAR --twopassMode Basic --runThreadN $THREADS --genomeDir $GENOME_INDEX \
  --readFilesCommand zcat \
  --readFilesIn $TR1 $TR2 \
  --outSAMstrandField intronMotif \
  --outSAMtype BAM Unsorted \
  --outFileNamePrefix $MAPDIR/${SAMPLE}_

# 6. Samtools sort/index
BAM="$MAPDIR/${SAMPLE}_Aligned.out.bam"
SORTBAM="$MAPDIR/${SAMPLE}_sorted.bam"
samtools sort -@ $THREADS -o $SORTBAM $BAM
samtools index $SORTBAM


# 7. Picard CollectAlignmentSummaryMetrics
picard CollectAlignmentSummaryMetrics \
  --INPUT $SORTBAM \
  --REFERENCE_SEQUENCE $REFERENCE \
  --OUTPUT $OUTPUT_METRICS_FILE

echo "Done $SAMPLE"

```

```bash

# i don't think the metrics thing worked but thats fine


# download stringtie into net 

https://github.com/gpertea/stringtie

git clone https://github.com/gpertea/stringtie
cd stringtie
make -j4 release


# add to path 

/n/netscratch/edwards_lab/Lab/kelsielopez/stringtie

export PATH=/n/netscratch/edwards_lab/Lab/kelsielopez/stringtie:$PATH



# next run stringtie on each sorted BAM file. producing GTF files (one per sample)

cd /n/netscratch/edwards_lab/Lab/kelsielopez/suboscines/annotation


# In a shell, or as a script loop:
# submit as job 


nano stringtie_bam_to_gtf.sh

#!/bin/bash
#SBATCH -t 0-12:00
#SBATCH -p test
#SBATCH -c 4
#SBATCH --mem=100G
#SBATCH -o stringtie_bam_to_gtf_thaDol_%j.out
#SBATCH -e stringtie_bam_to_gtf_thaDol_%j.err
#SBATCH --mail-type=END

for bam in mapped/*_sorted.bam
do
  sample=$(basename $bam _sorted.bam)
  stringtie "$bam" -o "${sample}.gtf" -p 4
done



# then merge all the gtfs
ls *.gtf > gtf.list


#Or, if your GTF files are in another directory, adjust as needed. gtf.list is used by StringTie below.

#Merge:
stringtie --merge -o RNA.merged.gtf gtf.list


# Step 3: (optional) Compare to reference (if you have a reference GTF
# gffcompare -o gffcompare_output RNA.merged.gtf


TransDecoder.LongOrfs -t transcripts.rna.fa
TransDecoder.Predict -t transcripts.rna.fa
/n/netscratch/edwards_lab/Lab/kelsielopez/TransDecoder-TransDecoder-v5.7.1/util/gtf_genome_to_cdna_fasta.pl
/n/netscratch/edwards_lab/Lab/kelsielopez/TransDecoder-TransDecoder-v5.7.1/util/gtf_to_alignment_gff3.pl
/n/netscratch/edwards_lab/Lab/kelsielopez/TransDecoder-TransDecoder-v5.7.1/util/cdna_alignment_orf_to_genome_orf.pl


# prepare transdecoder input . Get transcript sequences and convert merged GTF to GFF3
GENOME=/n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/Doliatus_pacbio/01_repeat_masking/03_known_out/ThaDol_18-293.p_ctg.known_mask.masked.fasta
/n/netscratch/edwards_lab/Lab/kelsielopez/TransDecoder-TransDecoder-v5.7.1/util/gtf_genome_to_cdna_fasta.pl RNA.merged.gtf $GENOME > transcriptome.fasta

/n/netscratch/edwards_lab/Lab/kelsielopez/TransDecoder-TransDecoder-v5.7.1/util/gtf_to_alignment_gff3.pl RNA.merged.gtf > transcriptome.gff3


# run transdecoder
TransDecoder.LongOrfs -t transcriptome.fasta
TransDecoder.Predict -t transcriptome.fasta
# This will produce protein and GFF3 outputs for candidate ORFs


# Step 6: Map ORF predictions back to genome space
/n/netscratch/edwards_lab/Lab/kelsielopez/TransDecoder-TransDecoder-v5.7.1/util/cdna_alignment_orf_to_genome_orf.pl \
    transcriptome.fasta.transdecoder.gff3 \
    transcriptome.gff3 \
    transcriptome.fasta > transcriptome_transdecoder.gff3
    
    
#Done.  44833 / 54468 transcript orfs could be propagated to the genome

    
```

# 5. Braker 

```bash

# downloading OrthoDB vertebrate protein fasta

py_env

conda activate snakemake

# i cloned this 
git clone https://github.com/tomasbruna/orthodb-clades.git

# changed into the directory and just ran this to get it to downlaod everything. apparantly the API is not great

(snakemake) [kelsielopez@boslogin08 orthodb]$ grep 'Vert' levels.tab
7742	Vertebrata	9754863	64382	470



(snakemake) [kelsielopez@boslogin08 orthodb]$ head levels.tab
2	Bacteria	66010659	644134	17551

python3 selectClade.py orthodb/raw.fasta orthodb/levels.tab orthodb/level2species.tab Vertebrata > vertebrata_proteins.fa

# then i can go and delete the raw fasta because it is huge 

```


```bash

# need to unlock snakemake if i cancelled the job when in the middle of running
snakemake --unlock

# now run it 
[kelsielopez@boslogin07 AnnotationBRAKER]$ head braker_snakemake_slurm_runner_thaDol.sh -n 100
#!/bin/bash
#SBATCH -J thaDol_brakersnake
#SBATCH -n 1                 
#SBATCH -t 72:00:00        
#SBATCH -p shared,edwards # add partition      
#SBATCH --mem=150G           
#SBATCH -o logs/test_thaDol.%A.out # need to create logs directory first  
#SBATCH -e logs/test_thaDol.%A.err  

#module purge
#module load python
#mamba activate snakemake # need to have created a snakemake conda environment
#profile=$1

snakemake --snakefile workflow/Snakefile --profile profiles/slurm




```

# 6. TOGA

```bash
# preparing for TOGA

### downloading ref genomes from ncbi
#####################################################################################################################
#####################################################################################################################
#####################################################################################################################

#!/bin/bash
# Pipeline: Project Zebra Finch (taeGut) annotation onto a target genome using TOGA
# Author: [your name]
# Date: 2024-06

############################
## ---------- USER-DEFINED VARIABLES TO EDIT ----------------
############################

# Base work directory for the project
PROJECT_BASE="/n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/annotation/02_toga_taeGut/thaDol"
# Zebra Finch NCBI accession and nicknames
REF_ACC="GCF_048771995.1"
REF_SHORT="taeGut"
# Target genome nickname and fasta path
TARGET_SHORT="thaDol"
TARGET_FA_ORIG="/n/netscratch/edwards_lab/Lab/kelsielopez/suboscines/Doliatus_pacbio/pacbio_hifi_assembly/workflow/results/assembly/ThaDol_18-293.p_ctg.fa"
# Location of make_lastz_chains and TOGA installs
MAKE_LASTZ_CHAINS_DIR="/n/netscratch/edwards_lab/Lab/kelsielopez/make_lastz_chains"
TOGA_DIR="/n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/TOGA"
NEXTFLOW_CONFIG_DIR="${TOGA_DIR}/nextflow_config_files"

# Conda environments for lastz and toga
ENV_LASTZ="nf_lastz"
ENV_TOGA="nf_toga"
ENV_NCBI="ncbi_datasets"
ENV_PY_ENV="python_env1"


# Set thread/resource requests for slurm jobs
SLURM_CPUS_CHAIN=4
SLURM_MEM_CHAIN=64000
SLURM_TIME_CHAIN="3-00:00"

SLURM_CPUS_TOGA=16
SLURM_MEM_TOGA=100000
SLURM_TIME_TOGA="3-00:00"

############################
## ---------- AUTOMATED VARIABLES (no user edit needed) ----------
############################

REF_DIR="${PROJECT_BASE}/reference"
TARGET_DIR="${PROJECT_BASE}/target"
CHAINS_DIR="${PROJECT_BASE}/chains"
TOGA_PROJECT_DIR="${PROJECT_BASE}/toga_project"

mkdir -p "$REF_DIR" "$TARGET_DIR" "$CHAINS_DIR" "$TOGA_PROJECT_DIR"

############################
## ---- STEP 1: Download reference genome & annotation ----
############################

conda activate ${ENV_NCBI}

cd "$REF_DIR"
# Download and unzip reference genome + annotation (skip if already done)
if [ ! -f "${REF_DIR}/ncbi_dataset.zip" ]; then
    datasets download genome accession $REF_ACC --include gff3,rna,cds,protein,genome,seq-report,gtf
    unzip ncbi_dataset.zip
fi

REF_NCBI_DATA="${REF_DIR}/ncbi_dataset/data/${REF_ACC}"
REF_FA_ORIG=$(ls $REF_NCBI_DATA/*.fna | head -n1)
REF_GTF=$(ls $REF_NCBI_DATA/*.gtf | head -n1)

# Clean fasta headers
awk '{if(/^>/){sub(/ .*/, "", $0)} print $0}' $REF_FA_ORIG > "$REF_DIR/${REF_SHORT}_genomic_simple.fna"
REF_FA="$REF_DIR/${REF_SHORT}_genomic_simple.fna"

# Clean headers from TARGET genome
awk '{if(/^>/){sub(/ .*/, "", $0)} print $0}' $TARGET_FA_ORIG > "$TARGET_DIR/${TARGET_SHORT}_simple.fa"
TARGET_FA="$TARGET_DIR/${TARGET_SHORT}_simple.fa"

############################
## ---- STEP 2: Prepare ref annotation (BED12, isoforms) ----
############################

cd $REF_DIR
# GTF → genePred → BED12

conda activate ${ENV_PY_ENV}

gtfToGenePred -ignoreGroupsWithoutExons $REF_GTF ${REF_SHORT}.genePred
genePredToBed ${REF_SHORT}.genePred ${REF_SHORT}.genePred.bed12.bed
REF_BED="$REF_DIR/${REF_SHORT}.genePred.bed12.bed"

# Make isoforms file for TOGA
awk '$3=="transcript" && /transcript_id/ && /gene_id/' $REF_GTF \
 | awk '{
    match($0, /gene_id "([^"]+)"/, g)
    match($0, /transcript_id "([^"]+)"/, t)
    if(length(g) && length(t)) print g[1] "\t" t[1]
   }' | sort | uniq > "$REF_DIR/isoforms.tsv"
sed -i '1igeneID\ttranscriptID' "$REF_DIR/isoforms.tsv"
ISOFORM_FILE="$REF_DIR/isoforms.tsv"

############################
## ---- STEP 3: Make 2bit files ----
############################

cd $REF_DIR
faToTwoBit "$REF_FA" "${REF_DIR}/${REF_SHORT}.2bit"
REF_2BIT="${REF_DIR}/${REF_SHORT}.2bit"

cd $TARGET_DIR
faToTwoBit "$TARGET_FA" "${TARGET_DIR}/${TARGET_SHORT}.2bit"
TARGET_2BIT="${TARGET_DIR}/${TARGET_SHORT}.2bit"

############################
## ---- STEP 4: Make LASTZ chains (separate SLURM job) ----
############################

cat > $CHAINS_DIR/run_lastz_chains.slurm <<EOF
#!/bin/bash
#SBATCH -p shared
#SBATCH -c ${SLURM_CPUS_CHAIN}
#SBATCH -t ${SLURM_TIME_CHAIN}
#SBATCH --mem=${SLURM_MEM_CHAIN}
#SBATCH -o lastz_chains_%j.out
#SBATCH -e lastz_chains_%j.err
#SBATCH --mail-type=END

source ~/.bashrc
conda activate ${ENV_LASTZ}
cd ${MAKE_LASTZ_CHAINS_DIR}

./make_chains.py ${REF_SHORT} ${TARGET_SHORT} \\
  ${REF_FA} ${TARGET_FA} \\
  --project_dir ${CHAINS_DIR}/${REF_SHORT}_${TARGET_SHORT}_chains \\
  --job_time_req 12h --executor slurm --cluster_queue shared
EOF

echo "Submit the following to generate LASTZ chains and wait for the .final.chain.gz completion:"
echo "cd $CHAINS_DIR && sbatch run_lastz_chains.slurm"

# User must submit and wait for chain file to finish before running TOGA

############################
## ---- STEP 5: Run TOGA (separate SLURM job) ----
############################

cat > $TOGA_PROJECT_DIR/run_toga.slurm <<EOF
#!/bin/bash
#SBATCH -p edwards,shared
#SBATCH -c ${SLURM_CPUS_TOGA}
#SBATCH -t ${SLURM_TIME_TOGA}
#SBATCH --mem=${SLURM_MEM_TOGA}
#SBATCH -o toga_%j.out
#SBATCH -e toga_%j.err
#SBATCH --mail-type=END

source ~/.bashrc
conda activate ${ENV_TOGA}
cd ${TOGA_DIR}

CHAIN="${CHAINS_DIR}/${REF_SHORT}_${TARGET_SHORT}_chains/${REF_SHORT}.${TARGET_SHORT}.final.chain.gz"
RENAMED_BED="${REF_BED}"
REF_2BIT="${REF_2BIT}"
TARGET_2BIT="${TARGET_2BIT}"
ISOFORM_FILE="${ISOFORM_FILE}"

./toga.py \$CHAIN \$RENAMED_BED \\
\$REF_2BIT \$TARGET_2BIT \\
--kt \\
--nc ${NEXTFLOW_CONFIG_DIR} --cb 3,5 --cjn 500 \\
--ms \\
--isoforms \$ISOFORM_FILE \\
--project_dir ${TOGA_PROJECT_DIR}/toga_${REF_SHORT}_on_${TARGET_SHORT}
EOF

echo ""
echo "Once your chain file is generated, submit:"
echo "cd $TOGA_PROJECT_DIR && sbatch run_toga.slurm"
echo ""
echo "TOGA orthology results will be in:"
echo "${TOGA_PROJECT_DIR}/toga_${REF_SHORT}_on_${TARGET_SHORT}/orthology_classification.tsv"
```
# 7. Combine

```bash

# use EVidence modeler. 

```
# 8. Find Orthologs with other thamnophilus genomes? 
# 8.1 toga - projecting the best annotation?

#  Combine TOGA and RNAseq

This recipe describes how to combine TOGA annotation with RNAseq assembly.
NB: It assumes bed12 format of both annotations
all scripts can be found here: https://github.com/osipovarev/Annotation_scripts

```
db=entCya
TOGA=galGal6_to_entCya2.annotation.bed
RNA=?? ## stringtie->tama assembly prodiced before
```

## clean up transcript names in tama rnaseq annotation 
`rename_duplicated_id_annotation.py -a <(awk '{OFS="\t"}{$4="tama_mRNA"; print}' $RNA) > temp.rna.bed`


## assign gene names (runs overlapSelect, which can be found here: http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/)
`assign_gene_names.sh temp.rna.bed $TOGA named.query_isoforms.tsv named.$db.rna_final.bed`


## remove transcripts with identical CDS
`getUniqTranscripts.py -f <(cat $TOGA named.$db.rna_final.bed) > $db.toga_rnaseq.anno.bed 2> /dev/null`


## Prepare combined isoforms file
```
rename_duplicated_id_annotation.py -c 1 -a <(cut -f4  $db.toga_rnaseq.anno.bed | awk -F"_" '{print $1"\t"$0}' | grep tama ) > tama.isoforms.tsv

cat named.query_isoforms.tsv tama.isoforms.tsv > combined_isofomrs.tsv
```

## clean up
`rm  isoforms.tsv temp.rna.bed`

# my notes from using katya's to combine annotations. I had to change some things:
```bash
# copy over all python scripts to actual merge directory 
cp /n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/annotation/test_merge/katya/*py


/n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/annotation/04_merging



cp /n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/annotation/test_merge/katya/*py \
/n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/annotation/04_merging

# copy over the shell scripts 
cp /n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/annotation/test_merge/katya/*sh \
/n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/annotation/04_merging

# first directory to work in because i ran stirngtie for thaDol already 
/n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/annotation/04_merging/thaDol







						##### prepare stringtie/transdecoder file ####

db="thaDol"
indir="/n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/annotation/04_merging"

cd /n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/annotation/04_merging/${db}

# this file comes from the last step of transcoder after running string tie 
gff="/n/netscratch/edwards_lab/Lab/kelsielopez/suboscines/annotation/transcriptome_transdecoder.gff3"

cp ${gff} /n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/annotation/04_merging/${db}


gff="/n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/annotation/04_merging/${db}/transcriptome_transdecoder.gff3"
gff_no_prefix="/n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/annotation/04_merging/${db}/transcriptome_transdecoder_no_prefix.gff3"

# remove prefix if needed 
sed 's/^ThaDol_18-293#pri#//' ${gff} > ${gff_no_prefix}


# need to add this to the top i guess 
nano $gff_no_prefix

##gff-version 3


gff_no_prefix="/n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/annotation/04_merging/${db}/transcriptome_transdecoder_no_prefix.gff3"
gff_no_prefix_noName="/n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/annotation/04_merging/${db}/transcriptome_transdecoder_no_prefix_noName.gff3"


#py_env


# Remove Name="...":
sed 's/;Name="[^"]*"//g' "$gff_no_prefix" > ${gff_no_prefix_noName}


gff_no_prefix_noName="/n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/annotation/04_merging/${db}/transcriptome_transdecoder_no_prefix_noName.gff3"
genePred_no_prefix_noName="/n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/annotation/04_merging/${db}/transcriptome_transdecoder_no_prefix_noName.genePred"


gff3ToGenePred ${gff_no_prefix_noName} ${genePred_no_prefix_noName}

bed12_no_prefix_noName="/n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/annotation/04_merging/${db}/transcriptome_transdecoder_no_prefix_noName.bed12"

genePredToBed ${genePred_no_prefix_noName} ${bed12_no_prefix_noName}

bed_no_prefix_noName="/n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/annotation/04_merging/${db}/transcriptome_transdecoder_no_prefix_noName.bed"

cp ${bed12_no_prefix_noName} ${bed_no_prefix_noName}



RNA="/n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/annotation/04_merging/${db}/transcriptome_transdecoder_no_prefix_noName.bed"

python3 ${indir}/rename_duplicated_id_annotation.py -a <(awk '{OFS="\t"}{$4="stringtie_mRNA"; print}' $RNA) > temp.rna.bed




						##### get names from blastp for stringtie/transdecoder annotation ####

bedToGenePred temp.rna.bed temp.rna.genePred
grep stringtie temp.rna.genePred > not_assigned.gp


GENOME2BIT=/n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/annotation/02_toga_taeGut/ThaDol/target/thaDol.2bit

genePredToProt -includeStop -starForInframeStops \
  not_assigned.gp \
  $GENOME2BIT \
  not_assigned.aa.fasta

mkdir -p split_fasta results_blast
faSplit sequence not_assigned.aa.fasta 100 split_fasta/chunk


######
nano blast_stringtie_uniprot.sh

#!/bin/bash
#SBATCH --job-name=ortho_blast_uniprot
#SBATCH --output=logs/blast_%A_%a.out
#SBATCH --error=logs/blast_%A_%a.err
#SBATCH --array=1-91
#SBATCH --cpus-per-task=24
#SBATCH --time=24:00:00
#SBATCH --mem=30G
#SBATCH -p edwards,shared

FILE=$(ls split_fasta | sed -n ${SLURM_ARRAY_TASK_ID}p)

databse="/n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/annotation/test_merge/katya/uniprot_sprot"

blastp -evalue 1e-10 \
       -num_threads 24 \
       -db ${databse} \
       -query split_fasta/$FILE \
       -outfmt 6 \
       -out results_blast/${FILE%.fa}.BlastHits_out6.tsv



######

cat results_blast/*.tsv > merged.not_assigned.BlastHits_out6.tsv


UNIPROT=/n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/annotation/test_merge/katya/uniprot_sprot.fasta

python3 ${indir}/assign_genes_to_hits.py \
  -u "$UNIPROT" \
  -b <(python3 ${indir}/filter_blast_hits.py -b merged.not_assigned.BlastHits_out6.tsv -n 1) \
  -s like \
  > transc.gene.dict.csv


# this should work 

python_env1) [kelsielopez@holy8a26602 thaDol]$ 
(python_env1) [kelsielopez@holy8a26602 thaDol]$ head transc.gene.dict.csv
stringtie_mRNA_21497,SMARCA5-like_stringtie_mRNA_21497
stringtie_mRNA_18652,KLHL8-like_stringtie_mRNA_18652
stringtie_mRNA_15900,SMARCA5-like_stringtie_mRNA_15900
stringtie_mRNA_13178,KLHL8-like_stringtie_mRNA_13178
stringtie_mRNA_10710,SMARCA5-like_stringtie_mRNA_10710
stringtie_mRNA_8320,KLHL8-like_stringtie_mRNA_8320



ANNO=temp.rna.genePred

python3 ${indir}/renameToHLscaffolds.py \
  -c 1 \
  -a "$ANNO" \
  -d transc.gene.dict.csv \
  > geneNames.$ANNO
  



genePredToBed geneNames.temp.rna.genePred geneNames.temp.rna.bed



# prepare toga files

#
#
# /n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/annotation/04_merging/${db}
#

							# tae gut #
# copy raw output from toga 
run="taeGut"
# copy raw output from toga 
cp /n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/annotation/02_toga_${run}/ThaDol/toga_project/toga_${run}_on_thaDol/query_annotation.bed \
/n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/annotation/04_merging/${db}

# rename raw annotation output from toga 
mv /n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/annotation/04_merging/${db}/query_annotation.bed \
/n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/annotation/04_merging/${db}/${run}_to_thaDol_query_annotation.bed


 # copy raw output from toga 
cp /n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/annotation/02_toga_${run}/ThaDol/toga_project/toga_${run}_on_thaDol/query_isoforms.tsv \
/n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/annotation/04_merging/${db}

# rename raw isoform output from toga 
mv /n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/annotation/04_merging/${db}/query_isoforms.tsv \
/n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/annotation/04_merging/${db}/${run}_to_thaDol_query_isoforms.tsv


							# gal gal #
							
#/n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/annotation/02_toga_galGal/thaDol/toga_project/toga_galGal_on_thaDol


run="galGal"
# copy raw output from toga 
cp /n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/annotation/02_toga_${run}/thaDol/toga_project/toga_${run}_on_thaDol/query_annotation.bed \
/n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/annotation/04_merging/${db}

# rename raw annotation output from toga 
mv /n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/annotation/04_merging/${db}/query_annotation.bed \
/n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/annotation/04_merging/${db}/${run}_to_thaDol_query_annotation.bed


 # copy raw output from toga 
cp /n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/annotation/02_toga_${run}/thaDol/toga_project/toga_${run}_on_thaDol/query_isoforms.tsv \
/n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/annotation/04_merging/${db}

# rename raw isoform output from toga 
mv /n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/annotation/04_merging/${db}/query_isoforms.tsv \
/n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/annotation/04_merging/${db}/${run}_to_thaDol_query_isoforms.tsv





cd /n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/annotation/04_merging/${db}


##
db=thaDol
# TOGA=taeGut_to_thaDol_query_annotation.bed

# not sure if this is totally necessary???
# ${indir}/assign_gene_names.sh geneNames.temp.rna.bed $TOGA taeGut_to_thaDol_query_isoforms.tsv named.$db.rna_final.bed


db=thaDol

TOGA1=taeGut_to_thaDol_query_annotation.bed
TOGA2=galGal_to_thaDol_query_annotation.bed

python3 ${indir}/getUniqTranscripts.py -f <(cat $TOGA1 ${TOGA2} geneNames.temp.rna.bed) > $db.toga_galGal_taeGut_rnaseq.anno.bed 2> /dev/null


python3 ${indir}/rename_duplicated_id_annotation.py -c 1 -a <(cut -f4  $db.toga_galGal_taeGut_rnaseq.anno.bed | awk -F"_" '{print $1"\t"$0}' | grep stringtie ) > stringtie.isoforms.tsv

# ^ not sure what to do about this with both togas 

cat taeGut_to_thaDol_query_isoforms.tsv galGal_to_thaDol_query_isoforms.tsv stringtie.isoforms.tsv > combined_isoforms.tsv




bedToGenePred $db.toga_galGal_taeGut_rnaseq.anno.bed $db.toga_galGal_taeGut_rnaseq.anno.genePred

genePredToGtf file \
  $db.toga_galGal_taeGut_rnaseq.anno.genePred \
  $db.toga_galGal_taeGut_rnaseq.anno.gtf


```

# run busco on combined annotation
```bash



py_env

#/n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/annotation/02_toga_taeGut/ThaDol/target/thaDol_forLASTZ.fa

REF="/n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/annotation/02_toga_taeGut/ThaDol/target/thaDol_forLASTZ.fa"
gb="thaDol"
indir="/n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/annotation/04_merging/${db}"


gffread ${indir}/$db.toga_galGal_taeGut_rnaseq.anno.gtf -T -o- | gffread - -F -o ${indir}/$db.toga_galGal_taeGut_rnaseq.anno.gff3
gffread ${indir}/$db.toga_galGal_taeGut_rnaseq.anno.gff3 -g ${REF} -w ${indir}/$db.toga_galGal_taeGut_rnaseq.anno.fasta


# to run busco 
#${indir}/$db.toga_galGal_taeGut_rnaseq.anno.fasta


cd /n/netscratch/edwards_lab/Lab/kelsielopez/busco-5.8.3

nano combined_annotations_thaDol_busco.sh

#!/bin/bash
#SBATCH -p edwards,shared
#SBATCH -c 1
#SBATCH -t 3-00:00
#SBATCH --mem=100000
#SBATCH --mail-type=END
#SBATCH -o combined_annotations_thaDol_busco_%j.out
#SBATCH -e combined_annotations_thaDol_busco_%j.err

# enter py_env to load dependencies 


input_files=(
"/n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/annotation/04_merging/thaDol/thaDol.toga_galGal_taeGut_rnaseq.anno.fasta"
)

BUSCO_DIR="/n/netscratch/edwards_lab/Lab/kelsielopez/busco-5.8.3/bin"

for input_file in "${input_files[@]}"; do
  base_name=$(basename "$input_file" .fasta)
  echo "Running BUSCO on $base_name"
  $BUSCO_DIR/busco -i "$input_file" -l aves_odb10 -o "${base_name}_thaDol_galGal_taeGut_rnaseq_busco_aves_odb10" -m transcriptome -f
done


```



# 8.2  fastOMA

```bash
# download and do test run

conda activate env_nf_rna

cd /n/netscratch/edwards_lab/Lab/kelsielopez


git clone https://github.com/DessimozLab/FastOMA.git
cd FastOMA

nextflow run FastOMA.nf -profile docker --container_version "sha-$(git rev-list --max-count=1 --abbrev-commit HEAD)" ...



# this works !

# run test data


# enter interactive sessino using 10 cores beacuse test data requires it 

salloc -p test -t 120:00 --mem=100000 -c 12


cd FastOMA/testdata


nextflow run ../FastOMA.nf  \
         --input_folder in_folder  \
         --omamer_db in_folder/omamerdb.h5 \
         --output_folder out_folder \
         --report \
         -profile singularity





```


```bash
# take toga bed output and convert to gtf. then predict ORF using transdecorder to make a protein fasta file.
#############################################################################################
#############################################################################################
							### dry Squ ###
#############################################################################################
#############################################################################################

# need to generate protein fasta files...

cd /n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/annotation/02_toga_taeGut/drySqu/toga_project/toga_taeGut_on_drySqu

py_env

# these are failing to download but these are what i need for converting bed to gff / gtf 
conda activate nf_lastz

# this is how to convert the toga output into gtf / gff3 format 
# Convert BED12 to genePred
bedToGenePred query_annotation.bed query_annotation.genePred

# Convert genePred to GTF
genePredToGtf file query_annotation.genePred query_annotation.gtf

# Optionally: genePred to GFF3 (newer version of ucsc tools may have genePredToGtf -gff3)
# genePredToGtf -gff3 file query_annotation.genePred query_annotation.gff3
# ^ THAT DIDN"T WORK!

# back to py_env
py_env

cd /n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/annotation/02_toga_taeGut/drySqu/toga_project/toga_taeGut_on_drySqu


# this doesn't have the prefix on it. the original thaDol genome has # # prefixes

genome="/n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/annotation/02_toga_taeGut/drySqu/target/drySqu_forLASTZ.fa"

# use gffread to convert gtf to fa
gffread query_annotation.gtf -g ${genome} -w transcripts.rna.fa


# Predict ORFs and translate to protein sequences with TransDecoder:

TransDecoder.LongOrfs -t transcripts.rna.fa
TransDecoder.Predict -t transcripts.rna.fa

# retain just the transcript ID here...
awk '/^>/{print ">"substr($1,2)} !/^>/' transcripts.rna.fa.transdecoder.pep > clean_for_fastoma.fa






# prepare protein files and splice files for fastOMA


ls /n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/annotation/02_toga_taeGut/*/toga_project/*/clean_for_fastoma.fa


wc -l /n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/annotation/02_toga_taeGut/*/toga_project/*/clean_for_fastoma.fa

head /n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/annotation/02_toga_taeGut/*/toga_project/*/clean_for_fastoma.fa

cp /n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/annotation/02_toga_taeGut/drySqu/toga_project/toga_taeGut_on_drySqu/clean_for_fastoma.fa /n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/annotation/02_toga_taeGut/drySqu/toga_project/toga_taeGut_on_drySqu/drySqu.fa
cp /n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/annotation/02_toga_taeGut/sakCri/toga_project/toga_taeGut_on_sakCri/clean_for_fastoma.fa /n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/annotation/02_toga_taeGut/sakCri/toga_project/toga_taeGut_on_sakCri/sakCri.fa
cp /n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/annotation/02_toga_taeGut/ThaDol/toga_project/toga_taeGut_on_thaDol/clean_for_fastoma.fa /n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/annotation/02_toga_taeGut/ThaDol/toga_project/toga_taeGut_on_thaDol/thaDol.fa


cp /n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/annotation/02_toga_taeGut/drySqu/toga_project/toga_taeGut_on_drySqu/drySqu.fa /n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/orthology_inference/fast_oma/proteome
cp /n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/annotation/02_toga_taeGut/sakCri/toga_project/toga_taeGut_on_sakCri/sakCri.fa /n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/orthology_inference/fast_oma/proteome
cp /n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/annotation/02_toga_taeGut/ThaDol/toga_project/toga_taeGut_on_thaDol/thaDol.fa /n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/orthology_inference/fast_oma/proteome

# i need to make a splice file for fastOMA

awk '
    /^>/ {
        name = $1
        name = substr(name,2)          # Drop the >
        base = name
        # Group by gene ID: assume gene is first part before the first .
        split(name, arr, ".")
        gene = arr[1]
        a[gene] = (a[gene] ? a[gene]";"name : name)
    } 
    END {
        for (g in a) print a[g]
    }
' drySqu.fa > drySqu.splice


awk '
    /^>/ {
        name = $1
        name = substr(name,2)          # Drop the >
        base = name
        # Group by gene ID: assume gene is first part before the first .
        split(name, arr, ".")
        gene = arr[1]
        a[gene] = (a[gene] ? a[gene]";"name : name)
    } 
    END {
        for (g in a) print a[g]
    }
' thaDol.fa > thaDol.splice


awk '
    /^>/ {
        name = $1
        name = substr(name,2)          # Drop the >
        base = name
        # Group by gene ID: assume gene is first part before the first .
        split(name, arr, ".")
        gene = arr[1]
        a[gene] = (a[gene] ? a[gene]";"name : name)
    } 
    END {
        for (g in a) print a[g]
    }
' sakCri.fa > sakCri.splice


# this could have been one line
for fa in *.fa; do
  awk '/^>/{name=$1; name=substr(name,2); split(name,a,"."); gene=a[1]; map[gene]=(map[gene]?map[gene]";"name:name)} END{for (g in map) print map[g]}' "$fa" > "${fa%.fa}.splice"
done

species_tree.nwk

((thaDol,sakCri)inter1,drySqu)inter2;

cd ..
(python_env1) [kelsielopez@holy8a24402 in_folder]$ wget https://omabrowser.org/All/LUCA.h5

mv LUCA.h5 omamerdb.h5

mkdir splice 
cp proteome/*splice splice

cd FastOMA/testdata

# now run 
salloc -p test -t 120:00 --mem=100000 -c 12


nano fast_oma_test.sh


#!/bin/bash
#SBATCH -p test
#SBATCH -c 16
#SBATCH -t 0-12:00
#SBATCH --mem=100G
#SBATCH -o fast_oma_test_%j.out
#SBATCH -e fast_oma_test_%j.err
#SBATCH --mail-type=END

conda activate env_nf_rna # need to be in this because it has profile singularity

fastOMA_path="/n/netscratch/edwards_lab/Lab/kelsielopez/FastOMA"
workdir="/n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/orthology_inference/fast_oma"
nextflow_slurm_config="/n/netscratch/edwards_lab/Lab/kelsielopez/FastOMA/nextflow_slurm.config"

cd ${workdir}

nextflow run ${fastOMA_path}/FastOMA.nf  \
         --input_folder in_folder  \
         --omamer_db in_folder/omamerdb.h5 \
         --output_folder out_folder \
         --report \
         -profile singularity


```






# X.X  PhyloAcc

```bash


```
