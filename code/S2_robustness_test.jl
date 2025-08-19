# Purpose: Compare the Julia PFIM metawebs with Dunhill 2024 webs
# Author: Tanya Strydom
# Copyright/License: MIT
# Contact: t.strydom@sheffield.ac.uk
# Date: 2025-08-19

# libraries
using AlgebraOfGraphics
using CairoMakie
using Combinatorics
using CSV
using DataFrames
using Extinctions
using pfim
using SpeciesInteractionNetworks

# internal functions
include("lib/internals.jl")

ints = CSV.read("data/dunhill/G1_extinction.csv", DataFrame)

# build DUNHILL network
S = unique(vcat(ints.resource, ints.consumer))
nodes = Unipartite(Symbol.(S))
edgs = Binary(zeros(Bool, (length(S), length(S))))

N_dunhill = SpeciesInteractionNetwork(nodes, edgs)

for i = 1:nrow(ints)
    interaction = (Symbol(ints.consumer[i]), Symbol(ints.resource[i]))
    N_dunhill[interaction...] = true
end

# remove cannibals
# N = remove_cannibals(N_dunhill)

robustness_vals = DataFrame(
    protect = Any[],
    mechanism = Any[],
    threshold = Any[],
    robustness = Any[],
);

# represents the % of species that have gone extinct (primary and secondary)
spread = collect(1:1:99)

# first we create `ext_reps` number of extinction series for each network
# note we do not protect the basal node

combos = [[:none; :cascade],
        [:none; :secondary],
        [:basal; :cascade],
        [:basal; :secondary]]
        
ext_reps = 20
        
for i = eachindex(combos)

    for l = 1:ext_reps
        # random extinction
        Ns = extinction(N; protect = combos[i][1], mechanism = combos[i][2])

        for j in eachindex(spread)

            rob = robustness(Ns; threshold = spread[j])
        
            D = DataFrame(
                protect = combos[i][1],
                mechanism = combos[i][2],
                threshold = spread[j],
                robustness = rob,
                )
        
            # send to results
            append!(robustness_vals, D)
        end
    end
end 

# robustness curves

layer = data(robustness_vals) * mapping(:threshold, :robustness, color = :protect => nonnumeric, linestyle = :mechanism) * (smooth())
fig = draw(layer)

figure = fig.figure