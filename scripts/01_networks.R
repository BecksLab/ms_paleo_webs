##Load other packages that are needed
library(igraph)
library(rgl)
library(tidyverse)

# functions needed to infer food webs
source("lib/pfim.R")
source("lib/niche.R")

# Data sources ----
# Data frame listing taxon names and life habits. 
# Taxon names must correspond to interaction data, and trait column names and
# levels must correspond to trait rules.

ex_taxonlist_g1 <- read.csv("data/community_1.csv")

#Four column df listing rules for feasible interactions.
ex_traitrules <- read.csv("data/feeding_rules.csv")

##Infer the food web using pfim model ----
pfim_web <- pfim(data = ex_taxonlist_g1,
                 col_taxon = "Guild",
                 cat_combo_list = ex_traitrules,
                 cat_trait_types = NULL, 
                 certainty_req = "all",
                 return_full_matrix = FALSE, 
                 print_dropped_taxa = TRUE
)

##Infer the food web using niche model ----

niche_web <- niche_model(nrow(ex_taxonlist_g1), 0.1)
