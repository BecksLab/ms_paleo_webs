# Are We Mapping Ecosystems or Models? How Reconstruction Choices Shape Food Web Inference

This repository contains the analytical pipeline used to test how food-web reconstruction choices affect ecological inference in deep-time ecosystems. The project uses the Early Jurassic (Toarcian) Oceanic Anoxic Event as a case study to compare multiple network reconstruction approaches and evaluate how assumptions about interaction inference influence inferred ecosystem structure, extinction dynamics, and recovery.

## Overview

This project evaluates whether conclusions about paleoecosystem dynamics are robust to the method used to reconstruct trophic interactions. The workflow compares multiple food-web models, quantifies network structure and interaction turnover, and asks how model choice changes conclusions about biodiversity, extinctions, resilience, and ecosystem reorganisation.

Key features:
- Julia-based food-web generation and ecological simulations
- R-based statistical and plotting workflows
- Quarto manuscript rendering for the project outputs
- Reproducible project environment via `code/Project.toml` and `code/Manifest.toml`

## Repository layout

- `code/` — main analysis scripts and project environment
  - `01_build_networks.jl` and related scripts for model construction
  - `02_topology.jl`, `03_extinctions.jl`, `04_extinction_analysis.jl`, and others for downstream analysis
  - `S1_...`, `S2_...`, `S3_...`, `S4_...` — supplementary analyses
- `data/` — input and processed datasets
  - `data/raw/` — species, trait, and guild data
  - `data/processed/` — generated networks and analysis products (large files are not generally tracked in the repo)
- `notebooks/` — outputs and supporting analyses
- `index.qmd` and other Quarto files — manuscript source and rendering files
- `figures/`, `tables/` — generated figures and summary tables

## Requirements

- Julia 1.9+
- Quarto
- R (for the R-based analysis and plotting scripts)
- GitHub access for any required packages hosted outside the standard Julia registry

## Quickstart

1. Clone the repository and open a terminal in the project root.
2. Install the Julia project dependencies:

   ```bash
   julia --project=code -e 'using Pkg; Pkg.instantiate()'
   ```

3. If needed, install the GitHub-hosted Julia packages used by the project:

   ```julia
   using Pkg
   Pkg.add("https://github.com/BecksLab/pfim.jl")
   Pkg.add("https://github.com/BecksLab/Extinctions.jl")
   ```

4. Run key analysis steps in order (examples):

   ```bash
   julia --project=code code/01_build_networks.jl
   julia --project=code code/02_topology.jl
   julia --project=code code/03_extinctions.jl
   ```

5. Render the manuscript and supplementary materials:

   ```bash
   quarto render .
   ```

Notes:
- Scripts are generally written to be run in sequence, but individual steps can be re-run as needed.
- Some analyses are computationally intensive and may be run locally or scheduled on a cluster.

## Reproducibility

The project uses the Julia environment defined in:
- `code/Project.toml`
- `code/Manifest.toml`

The Manifest is intended to preserve package versions used in the analysis. Keeping it intact helps ensure reproducible results.

## Data availability

Large processed output files are not stored in the repository. To regenerate them, run the relevant analysis scripts in `code/`.

Examples of generated outputs include:
- `data/processed/extinction_seq.jlds`
- other network and analysis products written during the pipeline

If a dataset is required for a study or replication workflow, please contact the project maintainer or consult the associated manuscript supplements.

## Citation

If you use the code, analyses, or results from this repository, please cite the associated manuscript when available.

## Contact

Maintainer: Tanya Strydom
