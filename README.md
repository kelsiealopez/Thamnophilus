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

# 7. Combine

# 8. Find Orthologs with other thamnophilus genomes? 
# 8.1 toga - pairwise
# 8.2 orthofinder - fastOMA
