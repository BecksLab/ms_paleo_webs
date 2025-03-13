using DataFrames
using EcologicalNetworksDynamics
using JLD2
using SpeciesInteractionNetworks

# set seed
import Random
Random.seed!(66)

A = [0 0 0; 1 0 0; 0 1 0] # 1 <- 2 <- 3.
foodweb = Foodweb(A)

N_niche = structuralmodel(NicheModel, 10, 0.4)

A = _get_matrix(N_niche)

Foodweb(A)

adbm_networks = load_object("data/processed/networks/adbm_networks.jlds")


N = adbm_networks.network[1]

A = _get_matrix(N)
Foodweb(A)