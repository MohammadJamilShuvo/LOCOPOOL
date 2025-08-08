#!/bin/bash
#SBATCH --job-name=sort_index_dedup_bam
#SBATCH --partition=cpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=6:00:00
#SBATCH --output=/pfs/work9/workspace/scratch/your_username-collembola/EnvilisPoolseq/11_logs/%x_%j.out
#SBATCH --error=/pfs/work9/workspace/scratch/your_username-collembola/EnvilisPoolseq/11_logs/%x_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=your_email@example.com

# -----------------------------
# Step 4.6: Sort & Index Deduplicated BAMs
# -----------------------------

cd /pfs/work9/workspace/scratch/your_username-collembola/EnvilisPoolseq || exit 1

# Activate environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate tools/envs/envilis_env

# Check samtools availability
if ! command -v samtools &> /dev/null; then
  echo "❌ samtools not found in envilis_env. Exiting."
  exit 1
fi

# Loop over deduplicated BAMs
for BAM in 03_alignments/*.dedup.bam; do
  SAMPLE=$(basename "$BAM" .dedup.bam)
  SORTED_BAM=03_alignments/${SAMPLE}.dedup.sorted.bam
  BAI_FILE=${SORTED_BAM}.bai

  echo "▶ Sorting $SAMPLE..."
  samtools sort -@ 2 -o "$SORTED_BAM" "$BAM" 2> 03_alignments/${SAMPLE}.sort.log

  if [[ -s "$SORTED_BAM" ]]; then
    echo "✅ Sorted BAM created for $SAMPLE"
    echo "▶ Indexing $SAMPLE..."
    samtools index "$SORTED_BAM"
    if [[ -s "$BAI_FILE" ]]; then
      echo "✅ Index created for $SAMPLE"
    else
      echo "❌ Index failed for $SAMPLE" >&2
    fi
  else
    echo "❌ Sorting failed for $SAMPLE — no output BAM written" >&2
  fi
done

echo "🏁 SLURM job finished on $(date)"
