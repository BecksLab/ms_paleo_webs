# Paleo Webs

A small little toy repo to act as a template to play around with some of the functionality afforded to us through GitHub.

The code itself is very utilitarian and simply uses a (modified) subset of the dataset from Dunhill et. al., (in prep) to construct food webs using the Paleo Food web Inference Model (PFIM, Shaw et. al., 2024) and the niche model (Williams and Martinez, 2000). We then compare the computed connectance and number of links for the same community by the different models.

## General use

This work has been developed within the `R` projects environment and it is possible to instantiate the project through `paleo_webs.Rproj`. Additionally, so as to maintain reproducibility, this project uses `{renv}`for dependency management and it is possible to download and install the relevant package versions by running `renv::restore()`.

## Folder and file structure

Data files are stored in `data/` which includes the trait data and feeding rules needed for the pfim model. Functions needed to implement the different models can be found in `lib/` and the actual workflow in `scripts/`.