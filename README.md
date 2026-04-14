# Integrating genomic and tagging data in European sea bass

**Gagnaire et al. (2026)**  
*Integrating genomic and tagging data reveals spatio-temporal population structure in Northeast Atlantic European sea bass*

---

## Overview

This repository contains the R scripts used to produce the six publication-ready figures in the paper above. The analysis integrates genome-wide SNP data with electronic-tagging data to characterise the spatio-temporal population structure of European sea bass (*Dicentrarchus labrax*) across the Northeast Atlantic and Mediterranean Sea.

The repository provides two equivalent entry points:

| File | Description |
|------|-------------|
| `Gagnaire_et_al_2026_figures.Rmd` | Self-contained R Markdown document. Knit it to reproduce all six figures at once. |
| `R/Figure1_sampling_map.R` … `R/Figure6_integrated.R` | Standalone R scripts, one per figure. Run them individually in any order. |

The standalone scripts are extracted **verbatim** from the Rmd (see [Changes from Rmd](#changes-from-rmd-to-standalone-scripts) for the three minor additions made to support standalone execution).

---

## Repository structure

```
.
├── Gagnaire_et_al_2026_figures.Rmd   # Full R Markdown document (all 6 figures)
├── R/
│   ├── Figure1_sampling_map.R        # Figure 1 – Sampling map
│   ├── Figure2_PCA.R                 # Figure 2 – PCA of SNP data
│   ├── Figure3_admixture.R           # Figure 3 – ADMIXTURE ancestry proportions
│   ├── Figure4_FST_IBD.R             # Figure 4 – FST heatmap & Isolation-by-Distance
│   ├── Figure5_tagging.R             # Figure 5 – Tagging trajectories & displacement
│   └── Figure6_integrated.R          # Figure 6 – Integrated genomic & tagging comparison
├── data/                             # Input data files (see Data section below)
└── figures/                          # Output figures (created automatically on first run)
```

---

## Requirements

### R version

R ≥ 4.2.0 is recommended.

### R packages

Install all required packages with:

```r
install.packages(c(
  "ggplot2", "dplyr", "tidyr", "readr", "tibble",
  "patchwork", "cowplot", "ggrepel",
  "sf", "rnaturalearth", "rnaturalearthdata",
  "RColorBrewer", "viridis", "scales",
  "adegenet", "hierfstat", "vegan",
  "reshape2", "lubridate", "geosphere", "fields"
))
```

> **Note:** `sf` requires system libraries (GDAL, GEOS, PROJ). On Ubuntu/Debian run  
> `sudo apt-get install libgdal-dev libgeos-dev libproj-dev` before installing from CRAN.

---

## Data

The scripts expect the following files inside a `data/` directory (relative to the working directory):

| File | Description |
|------|-------------|
| `genotypes.txt` | Tab-separated genotype matrix (rows = individuals, columns = SNP loci; allele dosage coded 0 / 1 / 2, missing = `NA`). |
| `metadata.csv` | One row per individual. Required columns: `individual_id`, `site`, `region`, `longitude`, `latitude`, `year`. |
| `tagging_data.csv` | One row per detection event. Required columns: `tag_id`, `release_lon`, `release_lat`, `release_region`, `recapture_lon`, `recapture_lat`, `date` (ISO 8601: `YYYY-MM-DD`). |
| `fst_matrix.csv` | Symmetric pairwise FST matrix; row and column names are population/site identifiers matching `metadata$site`. |
| `admixture_K2.txt` | Space-separated Q-matrix from ADMIXTURE (K = 2); rows must be in the same order as `metadata`. |
| `admixture_K3.txt` | Space-separated Q-matrix from ADMIXTURE (K = 3); rows must be in the same order as `metadata`. |
| `diversity_statistics.csv` | Per-population diversity statistics (used in the Rmd preamble; not required by the standalone scripts). |

Data files are available from [TODO: insert data repository DOI/URL once published].

---

## How to reproduce the figures

### Option A – Knit the R Markdown document (all figures at once)

1. Open `Gagnaire_et_al_2026_figures.Rmd` in RStudio (or run from the command line).
2. Make sure your working directory is the **root of this repository** (i.e., the folder containing `data/`).
3. Click **Knit** in RStudio, or run:

```r
rmarkdown::render("Gagnaire_et_al_2026_figures.Rmd")
```

Output figures are written to `figures/` in both PDF and PNG format.

### Option B – Run individual R scripts

Each script in `R/` is self-contained. To run a single figure script:

```r
# From an R session with the repository root as working directory:
source("R/Figure1_sampling_map.R")
```

Or from the shell:

```bash
Rscript R/Figure1_sampling_map.R
```

#### Changing input / output directories

The two lines at the top of every standalone script control the data and output locations:

```r
data_dir    <- "data"      # folder containing the input files
figures_dir <- "figures"   # folder where PDFs and PNGs will be saved
```

Edit these paths if your files are stored elsewhere. No other code needs to change.

---

## Changes from Rmd to standalone scripts

The standalone R scripts are extracted directly from the correspondingly-named `## FIGURE N ##` section of the Rmd with **three minimal additions** to support independent execution. These additions are clearly marked with `# [CHANGE N]` comments at the top of each script:

| Tag | Addition |
|-----|----------|
| `[CHANGE 1]` | `library()` calls moved to the top of the script. |
| `[CHANGE 2]` | Data-loading lines from the shared "Data loading" chunk of the Rmd. |
| `[CHANGE 3]` | Two configurable path variables (`data_dir`, `figures_dir`) and corresponding `file.path()` wrappers around hardcoded filenames. |

All plotting code is **unchanged** from the Rmd.

---

## Citation

If you use these scripts, please cite:

> Gagnaire P.-A. et al. (2026). Integrating genomic and tagging data reveals spatio-temporal population structure in Northeast Atlantic European sea bass. *[Journal TBC]*.

---

## License

Code is released under the [MIT License](LICENSE).

---

## Contact

Pierre-Alexandre Gagnaire – [p-a.gagnaire@ifremer.fr](mailto:p-a.gagnaire@ifremer.fr)  
MARBEC, University of Montpellier / CNRS / Ifremer / IRD, Montpellier, France
