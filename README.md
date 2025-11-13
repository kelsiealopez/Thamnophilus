# Thamnophilus

# repeat masking 

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
