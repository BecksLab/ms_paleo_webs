# Reconstructing food webs in deep time: why model choice matters for ecological inference

This repository provides the analytical pipeline for evaluating how model choice conditions ecological inference in paleoecological contexts. Using the Early Jurassic (Toarcian) Oceanic Anoxic Event as a case study, the project compares six distinct network reconstruction models to determine how their underlying assumptions (mechanistic vs. trait-based vs. structural) shape our understanding of ancient extinction dynamics.

## Scientific Context
Reconstructing paleo-food webs is a challenge because species interactions are rarely preserved. This project bridges paleoecology and modern network theory by treating reconstructions as alternative ecological hypotheses. We evaluate how these hypotheses diverge regarding:

* **Network Topology**: Macro-scale (e.g., connectance) and meso-scale (motifs) structure.

* **Extinction Dynamics**: Robustness to primary extinctions and the scale of secondary cascading effects.

* **Interaction Turnover**: How model assumptions influence "rewiring" across successive community states.

## Prerequisites

**External Julia Packages**

Two critical components of this pipeline are hosted as independent packages on GitHub and must be installed manually, as they are not currently in the general Julia registry:

PFIM.jl: The core paleo modeling framework.

`pkg> add https://github.com/BecksLab/pfim.jl`

Extinctions.jl: Logic for simulating secondary extinction cascades.

`pkg> add https://github.com/BecksLab/Extinctions.jl`

**Environment Management (Julia)**

This project uses Julia's built-in package manager to ensure reproducibility. You will find two key files in the /code directory:

Project.toml: Defines the project dependencies and their versions.

Manifest.toml: A complete "snapshot" of every package and sub-dependency used in the original analysis. (Note: This is ignored by Git in some contexts but is vital for matching the exact state of the SpeciesInteractionNetworks and Graphs packages).

To instantiate the environment:

`using Pkg; Pkg.activate("code"); Pkg.instantiate()`

## Project Structure

**/code**

The analysis is designed to be run in numerical order:

* 01_build_networks.jl: Generates synthetic food webs using five different models: ADBM, ATN (L-matrix), Body Mass Ratio, Niche, and Random. It also incorporates the downsampled and full PFIM.

* 02_topology.jl: Calculates macro, meso (motifs), and micro-scale network metrics.

* 03_generate_extinctions.jl: Simulates primary extinction sequences based on different proposed extinction mechanisms.

* 04_extinction_analysis.jl: Evaluates robustness and secondary extinction counts.

* 05_betadiv.jl: Quantifies network dissimilarity (β OS).

* 06_structural_differences.R: Performs MANOVA and LDA to visualize network space.

* 07_beta_div.jl.R: Analyses interaction turnover between modelling frameworks.

* 08_dunhill_compare.R: Uses GAMs to evaluate how well models track temporal changes and compare extinction simulations.

* S1_bodysize_distribution.jl and S1_bodysize_distribution.R: effect size analyses of body size distributions.

**/code/lib**

Internal functions and model implementations:

* `adbm.jl`, `lmatrix.jl`, `bodymass.jl`, `niche.jl`, `random.jl`: The specific model engines.

* `internals.jl`: Core topological summary functions and TSS validation logic.

* `plotting_theme.R`: Centralized aesthetics, fonts, and color palettes for all figures.

**/notebooks**

Contains the "raw" outputs for the Supplementary Materials. Tables for effect sizes, ANOVA results from beta-diversity analysis, and model coefficients are exported here from the R scripts to be dynamically pulled into the final manuscript.

**Manuscript & Quarto**

The final paper and supplementary documents are rendered using Quarto (.qmd files located in the root).

The manuscript integrates code-generated figures and tables directly.

To render the document, run quarto render from the root directory.

## Data Note

Due to file size constraints, several large processed datasets are not included in the remote repository (tracked via .gitignore). You must run the Julia scripts locally to regenerate the following:

* data/processed/extinction_seq.jlds (extinction results)
