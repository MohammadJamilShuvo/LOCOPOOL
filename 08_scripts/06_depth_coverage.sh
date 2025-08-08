#!/bin/bash
#SBATCH --job-name=depth_coverage
#SBATCH --partition=cpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=2:00:00
#SBATCH --output=/pfs/work9/workspace/scratch/your_username-collembola/EnvilisPoolseq/11_logs/%x_%j.out
#SBATCH --error=/pfs/work9/workspace/scratch/your_username-collembola/EnvilisPoolseq/11_logs/%x_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=your_email@example.com

# -----------------------------
# Step 6: Coverage Estimation
# -----------------------------

cd /pfs/work9/workspace/scratch/your_username-collembola/EnvilisPoolseq || exit 1

# Activate environment with samtools
source ~/miniconda3/etc/profile.d/conda.sh
conda activate tools/envs/envilis_env

# Create output directory
mkdir -p 05_depth_coverage

# Output file
OUTFILE=05_depth_coverage/mean_depth.tsv
echo -e "Sample\tMean_Depth" > "$OUTFILE"

# Loop through filtered BAMs and calculate mean depth
for BAM in 04_bam_filtered/*.filtered.bam; do
  SAMPLE=$(basename "$BAM" .filtered.bam)
  MEAN=$(samtools depth "$BAM" | awk '{sum+=$3} END {if (NR>0) print sum/NR; else print 0}')
  echo -e "${SAMPLE}\t${MEAN}" >> "$OUTFILE"
done

echo "🏁 Coverage calculation finished on $(date)"
