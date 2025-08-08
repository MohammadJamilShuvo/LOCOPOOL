<p align="center">
  <img src="docs/banner.png" alt="LOCOPOOL" width="100%">
</p>

# LOCOPOOL

[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE) [![Platform: SLURM](https://img.shields.io/badge/Platform-SLURM-lightgrey)](#) [![Environment: Conda](https://img.shields.io/badge/Environment-Conda-blue)](#) [![Status: Active](https://img.shields.io/badge/Status-Active-brightgreen)](#)

LOCOPOOL is a reproducible pipeline for **low-coverage Pool-Seq analysis** of archived field samples preserved in ethanol, optimized for landscape genomics and population-level allele frequency estimation.

---

##  Quick Start

**1. Clone the repository**
```bash
git clone https://github.com/MohammadJamilShuvo/LOCOPOOL.git
cd LOCOPOOL
```
**2. Set up Conda environments**
```bash
conda env create -f envs/envilis_env.yml
conda env create -f envs/envilis_angsd_env.yml
conda env create -f envs/picard_env.yml
```
**3. Execute the pipeline**
See detailed, modular SLURM steps in [PIPELINE.md], or run the whole workflow:
```bash
chmod +x run_pipeline.sh
./run_pipeline.sh
```

---

##  Repository Structure

```
.
├── 08_scripts/          # Core processing scripts
├── envs/                # Conda environment YAMLs
├── docs/                # Documentation and banner image
├── README.md            # This file
├── PIPELINE.md          # Detailed step-by-step guide
├── LICENSE              # MIT License
├── .gitignore           # Ignores logs, data files, caches, etc.
├── run_pipeline.sh      # Automates full pipeline execution
```

---

##  Pipeline Overview

1. **Reference genome indexing**
2. **Read trimming & QC**
3. **Read alignment (BWA)**
4. **Marking duplicates & filtering**
5. **Depth coverage estimation**
6. **Allele frequency estimation using ANGSD**  
  *Optional*: SNP calling

All scripts are SLURM-ready and configured for ease of customization.

---

##  Case Study

Initially implemented for 95 pooled *Entomobrya nevilis* samples (collected in Central European forests, stored in 99% ethanol), but fully generalizable across taxa and landscapes.

---

##  Citation & License

- **License**: MIT — see [LICENSE](LICENSE) for details.  
- **Citation**: Pipeline will be archived on Zenodo with DOI for academic referencing.

---  
<p align="center"><em>Built for reproducibility. Adapt, use, and cite!</em></p>
