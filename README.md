# Paleo Webs

A small little toy repo to act as a template to play around with some of the functionality afforded to us through GitHub.

The code itself is very utilitarian and simply uses a (modified) subset of the dataset from Dunhill et. al., (in prep) to construct food webs using the Paleo Food web Inference Model (PFIM, Shaw et. al., 2024) and the niche model (Williams and Martinez, 2000).

## General use

This work has been developed within the `R` projects environment and it is possible to instantiate the project through `paleo_webs.Rproj`. Additionally, so as to maintain reproducibility, this project uses `{renv}`for dependency management and it is possible to download and install the relevant package versions by running `renv::restore()`.

## Folder and file structure

The code workflow is contained in the `scripts/` folder, with each file being numbered sequentially. Additional functions that are needed to run some of the models are contained in the `lib/` folder, each model is contained within its own file and additional helper functions are in `lib/internals.R`.

## Associated manuscript

Link to google doc: [https://docs.google.com/document/d/1Hxh1nWtgir2b5HIpmqKjHzFzLRC9vJLTPtHNR5Z_hiE/edit?usp=sharing](https://docs.google.com/document/d/1Hxh1nWtgir2b5HIpmqKjHzFzLRC9vJLTPtHNR5Z_hiE/edit?usp=sharing)
