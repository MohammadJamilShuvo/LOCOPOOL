#!/bin/bash
#SBATCH --job-name=trim_qc
#SBATCH --partition=cpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=256G
#SBATCH --time=1-00:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=your_email@example.com
#SBATCH --output=11_logs/%x_%j.out
#SBATCH --error=11_logs/%x_%j.err

cd /pfs/work9/workspace/scratch/your_username-collembola/EnvilisPoolseq

LOGDIR="11_logs"
mkdir -p "$LOGDIR"
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
LOGFILE="$LOGDIR/trim_qc_$TIMESTAMP.log"

source ~/miniconda3/etc/profile.d/conda.sh
conda activate genome_env

bash 08_scripts/01_trim_and_qc.sh 2>&1 | tee "$LOGFILE"
