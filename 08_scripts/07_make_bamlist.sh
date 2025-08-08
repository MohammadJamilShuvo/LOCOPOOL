#!/bin/bash
#SBATCH --job-name=make_bamlist
#SBATCH --partition=cpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:10:00
#SBATCH --output=/pfs/work9/workspace/scratch/your_username-collembola/EnvilisPoolseq/11_logs/%x_%j.out
#SBATCH --error=/pfs/work9/workspace/scratch/your_username-collembola/EnvilisPoolseq/11_logs/%x_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=your_email@example.com

# -----------------------------
# Step 7: Prepare bamlist.txt for ANGSD
# -----------------------------

cd /pfs/work9/workspace/scratch/your_username-collembola/EnvilisPoolseq || exit 1

# Activate env (not needed for this step, but safe)
source ~/miniconda3/etc/profile.d/conda.sh
conda activate tools/envs/envilis_env

# Create bamlist.txt file from filtered BAMs
ls 04_bam_filtered/*.filtered.bam > 04_bam_filtered/bamlist.txt

echo "✅ bamlist.txt created with $(wc -l < 04_bam_filtered/bamlist.txt) entries"
echo "🏁 Done on $(date)"
