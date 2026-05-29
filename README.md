# Reconstructing food webs in deep time: why model choice matters for ecological inference

This repository contains the full analytical pipeline used to evaluate how food-web reconstruction models influence ecological inference in paleoecosystems. Using the Early Jurassic (Toarcian) Oceanic Anoxic Event as a case study, multiple network reconstruction approaches are compared to assess how their assumptions shape predictions of ecosystem structure, extinction dynamics, and recovery.

## Quick overview

- Languages: Julia and R
- Primary analysis scripts: /code (numerically ordered)
- Notebooks and figures: /notebooks
- Data: /data (raw + processed; large processed files are not included)

## Quickstart (recommended)

1. Install Julia (recommended version: 1.9+).
2. Open a Julia REPL in the repository root and run:

   using Pkg
   Pkg.activate("code")
   Pkg.instantiate()

3. (One-off) Install required GitHub-only packages:

   pkg> add https://github.com/BecksLab/pfim.jl
   pkg> add https://github.com/BecksLab/Extinctions.jl

4. Run the analysis (examples):

   julia --project=code code/01_build_networks.jl
   julia --project=code code/02_topology.jl
   julia --project=code code/03_generate_extinctions.jl

5. Render the manuscript and Supplementary Materials (Quarto):

   quarto render .

Notes:
- Scripts are designed to be run in order, but individual steps can be re-run as needed.
- Long-running jobs can be run on a cluster or backgrounded locally.

## What the pipeline does

- Generates ensembles of synthetic food webs using multiple reconstruction models (ADBM, ATN, Body-mass ratio, Niche, Random, PFIM variants).
- Computes network metrics at macro, meso (motifs), and micro scales.
- Simulates primary extinction sequences and propagating secondary extinctions under multiple removal scenarios.
- Quantifies interaction beta-diversity and performs statistical comparisons (MANOVA, LDA, PERMANOVA, mixed models).

## Project structure (high level)

- /code: analysis scripts (01_... to S2_...)
  - /code/lib: reusable functions and model implementations
- /data
  - /data/raw: curated input trait and guild data
  - /data/processed: generated networks and analysis outputs (not tracked when large)
- /notebooks: supplementary outputs and tables used in the manuscript
- Quarto files (.qmd) in the repository root for manuscript rendering

## Reproducibility

- Use the Project.toml and Manifest.toml inside /code to reproduce the Julia environment exactly.
- The Manifest snapshots package versions; keep it intact to ensure identical analysis results.
- External packages hosted on GitHub (PFIM.jl and Extinctions.jl) must be installed separately (see Quickstart).

## Data availability

Large processed files are not stored in the repo. To reproduce results, run the relevant Julia scripts to regenerate the following (examples):

- data/processed/extinction_seq.jlds

If you need access to the generated datasets, contact the repository maintainers or check any supplementary data archive associated with the manuscript.

## Contributing

Contributions and bug reports are welcome. Please open issues or pull requests describing the change and the rationale.

## Citation

If using results or code from this repository, please cite the associated manuscript (when available).

## Contact

Maintainer: Tanya STrydom