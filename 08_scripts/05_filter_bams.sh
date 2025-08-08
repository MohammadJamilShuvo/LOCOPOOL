#!/bin/bash
#SBATCH --job-name=filter_bams
#SBATCH --partition=cpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=4:00:00
#SBATCH --output=/pfs/work9/workspace/scratch/your_username-collembola/EnvilisPoolseq/11_logs/%x_%j.out
#SBATCH --error=/pfs/work9/workspace/scratch/your_username-collembola/EnvilisPoolseq/11_logs/%x_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=your_email@example.com

# -----------------------------
# Step 5: Filter Deduplicated Sorted BAM Files
# -----------------------------

cd /pfs/work9/workspace/scratch/your_username-collembola/EnvilisPoolseq || exit 1

# Activate environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate tools/envs/envilis_env

# Create output folder
mkdir -p 04_bam_filtered

# Loop through dedup.sorted.bam files
for BAM in 03_alignments/*.dedup.sorted.bam; do
  SAMPLE=$(basename "$BAM" .dedup.sorted.bam)
  OUT_BAM=04_bam_filtered/${SAMPLE}.filtered.bam

  echo "▶ Filtering $SAMPLE..."

  samtools view -@ 2 -b     -q 30     -f 0x2     -F 0x4 -F 0x100 -F 0x800     "$BAM" > "$OUT_BAM"

  if [[ -s "$OUT_BAM" ]]; then
    echo "✅ Filtered BAM created for $SAMPLE"
    samtools index "$OUT_BAM"
  else
    echo "❌ Filtering failed for $SAMPLE (empty or missing output)" >&2
  fi
done

echo "🏁 SLURM job finished on $(date)"
