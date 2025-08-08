#!/bin/bash
#SBATCH --job-name=angsd_snp_calling
#SBATCH --output=11_logs/angsd_snp_calling_%j.out
#SBATCH --error=11_logs/angsd_snp_calling_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=6:00:00
#SBATCH --partition=cpu

# Load environment and activate ANGSD
source ~/miniconda3/etc/profile.d/conda.sh
conda activate tools/envs/envilis_angsd_env

# Run ANGSD SNP calling
angsd -b 04_bam_filtered/bamlist.txt \
  -ref ref/entomobrya_nevilis_genome.fa \
  -out 06_freq_results/nevilis_snps \
  -doMaf 1 -GL 1 -doMajorMinor 1 \
  -minMapQ 30 -minQ 20 \
  -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 \
  -SNP_pval 1e-6 -minInd 20 \
  -P 4
