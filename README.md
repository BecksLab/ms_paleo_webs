# Reconstructing food webs in deep time: why model choice matters for ecological inference

This repository contains the full analytical pipeline used to evaluate how food-web reconstruction models influence ecological inference in paleoecosystems. Using the Early Jurassic (Toarcian) Oceanic Anoxic Event as a case study, we compare multiple network reconstruction approaches to assess how their underlying assumptions shape predictions of ecosystem structure, extinction dynamics, and recovery.

## Scientific Context

Reconstructing paleo-food webs is inherently uncertain because species interactions are rarely preserved in the fossil record. This project treats different reconstruction models as alternative ecological hypotheses and evaluates how these hypotheses diverge in their predictions.
We focus on three key dimensions:

* **Network Topology**: Differences in macro-scale (e.g., connectance), meso-scale (motifs), and micro-scale (degree structure) properties.

* **Extinction Dynamics**: Variation in robustness to species loss and the magnitude of secondary extinction cascades.

* **Interaction Turnover**: The extent to which models differ in predicted interaction 'rewiring' across time.

## Network generation overview

For each reconstruction approach and time period, ensembles of food webs were generated under a unified framework:

* All models use a shared species pool and identical species richness.

* Each model–time combination is simulated across multiple stochastic replicates (100 per configuration).

* Networks are represented as binary directed adjacency matrices.

### Trait handling

* Body mass–based models (ADBM, ATN, body size ratio): Species body masses are stochastically assigned per replicate, generating variation in interaction structure.

* Trait-rule models (PFIM): Discrete ecological traits (e.g., guild, motility, tiering) are fixed within each time period, with variability arising from link filtering/downsampling.

* Structural models (Niche, Random): Interactions are generated from statistical rules independent of observed traits.

This design ensures that differences among networks reflect model assumptions rather than differences in species composition.


### Extinction simulations

Each network replicate is subjected to multiple extinction scenarios:

* Random removal

* Topological removal (based on generality or vulnerability)

* Trait-ordered removal (based on ecological hierarchies)

Extinctions are simulated dynamically:
* Primary removals are applied sequentially

* Secondary extinctions occur when species lose all resources

* Each scenario is repeated across replicates (~50 extinction runs per network)

### Beta diversity and structural comparison

Network dissimilarity is quantified using interaction beta diversity:
* βS: total dissimilarity

* βWN: rewiring among shared species

* βOS: species turnover component

Structural differences are further analysed using:
* MANOVA and canonical discriminant analysis (CDA)

* Linear discriminant analysis (LDA)

* Mixed-effects models and ANOVA

* Variance partitioning (PERMANOVA)

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

* `01_build_networks.jl`: Generates synthetic food webs using five different models: ADBM, ATN (L-matrix), Body Mass Ratio, Niche, and Random. It also incorporates the downsampled and full PFIM.

* `02_topology.jl`: Calculates macro, meso (motifs), and micro-scale network metrics.

* `03_generate_extinctions.jl`: Simulates primary extinction sequences based on different proposed extinction mechanisms.

* `04_extinction_analysis.jl`: Evaluates robustness and secondary extinction counts.

* `05_betadiv.jl`: Quantifies network dissimilarity (β OS).

* `06_structural_differences.R`: Performs MANOVA and LDA to visualize network space.

* `07_beta_div.jl.R`: Analyses interaction turnover between modelling frameworks.

* `08_dunhill_compare.R`: Temporal comparisons and compare extinction simulations.

* `09_variance_partitioning.R`: Variance partitioning between model and time.

* `08_ANOVA_time.R`: Two-way ANOVA and post hoc analyses to quantify how food web structure varies across model and time.

* `S1_bodysize_distribution.jl` and `S1_bodysize_distribution.R`: effect size analyses of body size distributions and body mass ratio parameters.

* `S2_pfim_downsample.jl` and `S2_pfim_downsample.R`: effect of downsampling on networks generated with PFIM.

**/code/lib**

Internal functions and model implementations:

* `adbm.jl`, `lmatrix.jl`, `bodymass.jl`, `niche.jl`, `random.jl`: The specific model engines.

* `internals.jl`: Core topological summary functions and TSS validation logic.

* `plotting_theme.R`: Centralized aesthetics, fonts, and color palettes for all figures.

**/data**

* `/raw`: input datasets (traits, guilds)

* `/processed`: generated networks and analysis outputs

⚠️ Large processed files are not tracked in the repository.

**/notebooks**

Contains the outputs for the Supplementary Materials. Tables for effect sizes, ANOVA results from beta-diversity analysis, and model coefficients are exported here from the R scripts to be dynamically pulled into the final manuscript.

**Manuscript & Quarto**

The final paper and supplementary documents are rendered using Quarto (.qmd files located in the root).

The manuscript integrates code-generated figures and tables directly.

To render the document, run quarto render from the root directory.

## Data Note

Due to file size constraints, several large processed datasets are not included in the remote repository (tracked via .gitignore). You must run the Julia scripts locally to regenerate the following:

* data/processed/extinction_seq.jlds (extinction results)
