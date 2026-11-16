# Are We Mapping Ecosystems or Models? How Reconstruction Choices Shape Food Web Inference

This repository contains the analytical pipeline, data, and manuscript sources for a study testing how different food-web reconstruction approaches alter ecological inference in deep-time ecosystems. The project uses the Early Jurassic Toarcian Oceanic Anoxic Event as a case study to compare alternative reconstruction frameworks and assess how assumptions about interaction inference change conclusions about ecosystem structure, extinction dynamics, and recovery.

The central question is whether biological conclusions are robust to the method used to reconstruct trophic interactions, or whether model choice itself drives apparent ecological patterns.

## Project overview

The project compares multiple food-web reconstruction approaches, including:

- allometric and trait-based network models
- body-size based reconstruction
- niche and random models
- PFIM-style paleo food-web inference
- downsampled and parameter-perturbed variants

For each model, it generates replicate networks, compares their topology, simulates extinctions, and evaluates how much model choice affects inferred ecological structure and time-series signals. The work is organized around both the main manuscript and supplementary analyses and is built to be reproducible from raw data and scripted analysis steps.

## Repository structure

- `code/` — main analysis scripts and project environment
  - `01_build_networks.jl` — generate synthetic food-web networks across models and time slices
  - `02_topology.jl` — compute network topology metrics
  - `03_extinctions.jl` — simulate extinction sequences under multiple mechanisms
  - `04_extinction_analysis.jl` — compare simulated results against empirical recovery transitions
  - `05_betadiv.jl` — analyse beta diversity / interaction turnover between models
  - `06_structural_differences.R` and other `*.R` scripts — downstream statistical analysis and plotting
  - `lib/` — shared helper functions and model-specific routines
  - `Project.toml` and `Manifest.toml` — Julia project environment
- `data/` — input and processed data
  - `data/raw/` — guild data, size classes, and feeding rules used to generate networks
  - `data/processed/` — generated network objects, topology summaries, extinction sequences, and derived outputs
- `notebooks/` — supplementary analyses and summary notebooks
  - `model_summaries.qmd` — expanded model descriptions and methodological notes
  - `results_expanded.qmd` — extended analysis and supplementary results
  - `figures/` and `tables/` — exported figures and summary tables
- `index.qmd` — main manuscript source
- `figures/` and `tables/` — rendered project output directories
- `_quarto.yml` — Quarto project configuration
- `requirements.txt` — Python environment used by the notebook workflow

## Research workflow

The project follows a reproducible research pipeline:

1. Build food-web networks from fossil guild data and trait rules.
2. Calculate topology and other structural metrics for each network.
3. Simulate extinction sequences under random and trait-structured scenarios.
4. Compare reconstructed models using multivariate statistics and beta-diversity metrics.
5. Evaluate whether model identity or time explains the strongest ecological signal.
6. Render the manuscript and supplementary material with Quarto.

## Requirements

To reproduce the analysis, you will need:

- Julia 1.9 or newer
- R and the R packages used by the analysis scripts
- Quarto
- Git

Optional but often useful:

- local or remote compute resources for heavier network simulations

## Setup

From the project root:

```bash
julia --project=code -e 'using Pkg; Pkg.instantiate()'
```

If any package is missing from the default registry, install it explicitly in the project environment. For example:

```julia
using Pkg
Pkg.add("https://github.com/BecksLab/pfim.jl")
Pkg.add("https://github.com/BecksLab/Extinctions.jl")
```

For the R workflow, install the packages required by the scripts in `code/*.R` if they are not already available.

## Quickstart

### 1. Generate the network ensemble

```bash
julia --project=code code/01_build_networks.jl
```

### 2. Compute topology summaries

```bash
julia --project=code code/02_topology.jl
```

### 3. Run extinction simulations

```bash
julia --project=code code/03_extinctions.jl
```

### 4. Analyse extinction outcomes

```bash
julia --project=code code/04_extinction_analysis.jl
```

### 5. Quantify model agreement / beta diversity

```bash
julia --project=code code/05_betadiv.jl
```

### 6. Run downstream R analyses and generate figures

```bash
Rscript code/06_structural_differences.R
Rscript code/07_beta_div.jl.R
Rscript code/08_dunhill_compare.R
Rscript code/09_variance_partitioning.R
Rscript code/10_ANOVA_time.R
```

The exact R scripts used in the workflow depend on which analysis you want to reproduce; several supplementary analyses are also present in the `code/` directory.

### 7. Render the manuscript and supplements

```bash
quarto render .
```

This will render the main manuscript and configured Quarto documents, including the project outputs and notebook-based supplementary material.

## Data notes

- Raw data are stored in `data/raw/`.
- Processed outputs such as network objects, topology tables, and extinction summaries are written to `data/processed/`.
- Large generated products are not always tracked in git and may need to be regenerated by rerunning the relevant analysis scripts.

## Main outputs

The workflow produces a variety of outputs, including:

- simulated network ensembles across reconstruction models and time periods
- network topology summaries
- extinction trajectories and robustness metrics
- model-comparison statistics and multivariate analyses
- manuscript figures and supplementary tables

## Reproducibility and environment

The Julia project environment is pinned via:

- `code/Project.toml`
- `code/Manifest.toml`

Keeping both files intact helps preserve the package versions used in the analysis and supports reproducible results across machines.

## Citation

If you use this repository or its methodological pipeline in your work, please cite the associated manuscript when available.

## Contact

Maintainer: Tanya Strydom

## Notes

This project is designed to support a comparative ecological-network workflow rather than a general-purpose package. The scripts assume a specific raw-data layout and a particular sequence of analysis steps, so rerunning them works best when the repository structure is preserved.
