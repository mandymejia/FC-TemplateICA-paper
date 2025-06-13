# FC-TemplateICA Paper Reproducibility Repository

This repository contains code to reproduce all figures and analyses from the FC-TemplateICA paper, organized into **simulation studies** and **real data analysis** using Human Connectome Project (HCP) data.

## Repository Structure

- `simulation/`: Simulation studies comparing FC-tICA with other methods using simulated data
- `code/`: Real data analysis using preprocessed HCP resting-state fMRI data  
- Generated plots and results are saved throughout the repository.

## Quick Start

1. **Setup**: Review `0_setup.R` in both `simulation/` and `code/` directories for required R packages and file paths
2. Specific figures are commented in the code (look for `#%# Creates Fig. X.X` comments)

## Data Notes

Intermediate results and templates are included in the repository to enable reproduction without requiring access to the full HCP dataset. For simulations, we supply `TCs.RDS` but users can generate their own time courses using `simulation/00_get_TCs_HCP.R` if desired.

**Full reproduction from scratch**:

**Full reproduction from scratch**:

- **Simulation analyses**: Execute simulation scripts in order (`1_make_templates.R` → `2_run_models.R` → `3_collect_results.R` → `4_visualize.R`) from the `simulation/` directory.

- **HCP real data analysis**: Execute `code/1_analysis.R` with HCP data.

> **Note**: For either path, set flags like `first_run` and other control flags at the beginning of scripts to `TRUE` to enable full computation from scratch.
