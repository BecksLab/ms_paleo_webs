using DataFrames
using EcologicalNetworksDynamics
using JLD2
using SpeciesInteractionNetworks

# set seed
import Random
Random.seed!(66)

A = [0 0 0; 1 0 0; 0 1 0] # 1 <- 2 <- 3.
foodweb = Foodweb(A)

A = _get_matrix(N_niche)

Foodweb(A)

adbm_networks = load_object("data/processed/networks/adbm_networks.jlds")


N = adbm_networks.network[1]

N_niche = structuralmodel(NicheModel, SpeciesInteractionNetworks.richness(N), connectance(N))

fw_adbm = Foodweb(_get_matrix(N))
fw_niche = Foodweb(_get_matrix(N_niche))

m2 = default_model(Foodweb([3 => 2, 2 => 1]), BodyMass([1.5, 1.2, 1.5]))
m2.body_mass
m2.species.names