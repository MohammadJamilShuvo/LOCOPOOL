# Envilis Pool-Seq Workflow

This document provides a complete reproducible bioinformatics pipeline for pooled low-coverage WGS data using ANGSD and related tools.

## Step-by-Step Workflow

1. **Setup and Env Creation**
   - Run `tools/conda_setup.sh` to install miniconda and create all environments

2. **Reference Genome Indexing**
   - `submit_ref_index.sh` downloads and indexes the *E. nevilis* genome

3. **Trimming & QC**
   - `submit_trim_qc.sh` executes `01_trim_and_qc.sh`
   - Visual summaries via `01a_plot_fastp_summary.py`

4. **Read Mapping**
   - `03_map_reads_updated.slurm` aligns trimmed reads with BWA

5. **Mark PCR Duplicates**
   - `04_5_mark_duplicates_verified_FINAL.sh` uses Picard

6. **Sort & Index Deduplicated BAMs**
   - `04_6_sort_index_dedup_bam_FINAL.sh`

7. **Filter BAMs**
   - `05_filter_bams.sh` keeps high-quality, properly paired reads

8. **Coverage Estimation**
   - `06_depth_coverage.sh`

9. **Prepare bamlist**
   - `07_make_bamlist.sh`

10. **Allele Frequency Estimation**
    - `08_angsd_maf_only_v2.sh` (for MAF)
    - `08_angsd_snp_calling.sh` (for SNPs)

---

Each SLURM script logs outputs in `11_logs/`. Scripts are modular and can be run stepwise.
