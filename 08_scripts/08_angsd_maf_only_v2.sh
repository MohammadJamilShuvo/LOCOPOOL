#!/bin/bash
#SBATCH --job-name=angsd_maf_only
#SBATCH --partition=cpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=2-00:00:00
#SBATCH --output=11_logs/angsd_maf_only_%j.out
#SBATCH --error=11_logs/angsd_maf_only_%j.err

echo "🧬 Running ANGSD to estimate allele frequencies (MAF only, no SFS)..."

source ~/miniconda3/etc/profile.d/conda.sh
conda activate tools/envs/envilis_angsd_env

angsd -bam 04_bam_filtered/bamlist.txt \
  -ref ref/entomobrya_nevilis_genome.fa \
  -out 06_freq_results/nevilis_maf \
  -GL 1 -doMajorMinor 1 -doMaf 1 \
  -minMapQ 20 -minQ 15 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 \
  -minInd 20 \
  -SNP_pval 1e-3 \
  -P 4

echo "✅ ANGSD MAF-only analysis completed."
