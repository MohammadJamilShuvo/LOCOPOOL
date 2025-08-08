#!/bin/bash
#SBATCH --job-name=mark_duplicates_verified
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
# Step 4.5: Mark PCR Duplicates
# -----------------------------

cd /pfs/work9/workspace/scratch/your_username-collembola/EnvilisPoolseq
echo "▶ SLURM job started in $(pwd) on $(date)"

# Activate working Picard conda env
source ~/miniconda3/etc/profile.d/conda.sh
conda activate picard_env

# Create writable temporary directory
TMPDIR=$(pwd)/tmp
mkdir -p "$TMPDIR"

# Check for input BAMs
shopt -s nullglob
FILES=(03_alignments/*.sorted.bam)
shopt -u nullglob

if [[ ${#FILES[@]} -eq 0 ]]; then
  echo "❌ No .sorted.bam files found in 03_alignments/. Exiting."
  exit 1
fi

# Run MarkDuplicates
for BAM in "${FILES[@]}"; do
  SAMPLE=$(basename "$BAM" .sorted.bam)
  DEDUP=03_alignments/${SAMPLE}.dedup.bam
  METRICS=03_alignments/${SAMPLE}.dup_metrics.txt

  if [[ -f "$DEDUP" && -f "$METRICS" ]]; then
    echo "⚠️  Skipping $SAMPLE — already processed"
    continue
  fi

  echo "▶ Marking duplicates for $SAMPLE..."

  picard MarkDuplicates \
    I="$BAM" \
    O="$DEDUP" \
    M="$METRICS" \
    TMP_DIR="$TMPDIR" \
    VALIDATION_STRINGENCY=LENIENT \
    REMOVE_DUPLICATES=false

  if [[ $? -eq 0 && -s "$DEDUP" ]]; then
    echo "✅ $SAMPLE processed successfully"
  else
    echo "❌ $SAMPLE failed to process or output was empty" >&2
  fi
done

echo "🏁 SLURM job finished on $(date)"
