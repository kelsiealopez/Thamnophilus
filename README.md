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

# 5. Braker 

# 6. TOGA

# 7. Combine

# 8. Find Orthologs with other thamnophilus genomes? 
# 8.1 toga - pairwise
# 8.2 orthofinder - fastOMA
