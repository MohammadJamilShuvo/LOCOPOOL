#!/bin/bash
#SBATCH --job-name=map_reads
#SBATCH --partition=cpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=1-00:00:00
#SBATCH --output=11_logs/%x_%j.out
#SBATCH --error=11_logs/%x_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=your_email@example.com

echo "▶ SLURM job started on $(date)"

# Load Conda and activate environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate ./tools/envs/envilis_env

# Confirm tool versions
echo "▶ Using BWA version:"
bwa 2>&1 | head -n 1
echo "▶ Using samtools version:"
samtools --version | head -n 1

# Define reference genome path
REF="ref/entomobrya_nevilis_genome.fa"

# Ensure output directory exists
mkdir -p 03_alignments

# Loop over all trimmed read pairs
for R1 in 01_trimmed_reads/*_R1.trimmed.fastq.gz; do
  SAMPLE=$(basename "$R1" _R1.trimmed.fastq.gz)
  R2="01_trimmed_reads/${SAMPLE}_R2.trimmed.fastq.gz"

  # Check if R2 exists
  if [[ ! -f "$R2" ]]; then
    echo "⚠️  Missing R2 for $SAMPLE — skipping" >&2
    continue
  fi

  echo "▶ Mapping $SAMPLE at $(date)..."

  bwa mem -t 4 "$REF" "$R1" "$R2" | \
    samtools view -bS - > "03_alignments/${SAMPLE}.bam"

  if [[ $? -eq 0 ]]; then
    echo "✅ Done: $SAMPLE"
  else
    echo "❌ Error mapping $SAMPLE" >&2
  fi
done

echo "🏁 SLURM job completed on $(date)"
