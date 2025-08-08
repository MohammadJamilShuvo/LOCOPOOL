# Envilis Pool-Seq Pipeline

A reproducible workflow for SNP and allele frequency estimation from pooled low-coverage whole genome sequencing (Pool-Seq) data of *Entomobrya nevilis* in a landscape genomics context.

## 🧪 Project Summary

This pipeline processes 95 pooled population samples (6–10 individuals each) collected from the Black Forest (ConFoBi project) for landscape genomics analyses.

- Input: Paired-end FASTQ files (~6× coverage per pool)
- Output: Filtered BAMs, allele frequency tables, SNP calls, and QC plots

## 📁 Directory Layout

```
EnvilisPoolseq/
├── 00_raw_reads/             # Input FASTQs
├── 01_trimmed_reads/         # Trimmed reads
├── 02_fastqc_reports/        # QC reports and summaries
├── 03_alignments/            # All BAMs (unsorted, dedup, etc.)
├── 04_bam_filtered/          # Filtered BAMs for ANGSD
├── 05_depth_coverage/        # Coverage estimates
├── 06_freq_results/          # ANGSD output
├── 07_plots/                 # QC plots
├── 08_scripts/               # All executable scripts
├── ref/                      # Reference genome
├── tools/                    # Conda environments and dependencies
└── 11_logs/                  # SLURM logs
```

## 🧬 Key Tools

- `fastp`, `jq`, `bwa`, `samtools`, `picard`, `ANGSD`, `realSFS`
- Visualization: `matplotlib`, `seaborn`
- Scheduler: SLURM

## 🛠️ Installation

```bash
git clone https://github.com/YOUR_USERNAME/envilis-poolseq-pipeline.git
cd envilis-poolseq-pipeline

# Follow PIPELINE.md to run all steps
```

## 📜 License

MIT License
