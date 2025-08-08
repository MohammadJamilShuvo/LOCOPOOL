#!/bin/bash
#SBATCH --job-name=ref_index
#SBATCH --partition=cpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=256G
#SBATCH --time=02:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=your_email@example.com
#SBATCH --output=11_logs/%x_%j.out
#SBATCH --error=11_logs/%x_%j.err

cd /pfs/work9/workspace/scratch/your_username-collembola/EnvilisPoolseq

LOGDIR="11_logs"
mkdir -p "$LOGDIR"
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
LOGFILE="$LOGDIR/ref_index_$TIMESTAMP.log"

source ~/miniconda3/etc/profile.d/conda.sh
conda activate genome_env

echo "Downloading genome..." | tee -a "$LOGFILE"
datasets download genome accession GCA_034695485.1 --filename entomobrya_genome.zip

echo "Unzipping genome..." | tee -a "$LOGFILE"
unzip entomobrya_genome.zip -d ref/
mv ref/ncbi_dataset/data/*/*.fna ref/entomobrya_nevilis_genome.fa

REF=ref/entomobrya_nevilis_genome.fa

echo "Running BWA index..." | tee -a "$LOGFILE"
bwa index $REF

echo "Running samtools faidx..." | tee -a "$LOGFILE"
samtools faidx $REF

echo "Creating sequence dictionary..." | tee -a "$LOGFILE"
gatk CreateSequenceDictionary -R $REF

echo "Step 2 completed." | tee -a "$LOGFILE"
