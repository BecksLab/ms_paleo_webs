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
using JLD2
using SpeciesInteractionNetworks

# internal functions
include("lib/internals.jl")

# import networks
networks = load_object("data/processed/networks.jlds")

# we only care about pfim metawebs here...
filter!(:model => x -> x == "pfim_metaweb", networks)

robustness_vals = DataFrame(
    protect = Any[],
    mechanism = Any[],
    time = Any[],
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
        
ext_reps = 3

for l in 1:nrow(networks)

    # remove cannibals
    N = remove_cannibals(networks.network[l])

    for i = eachindex(combos)
    
        for l = 1:ext_reps
            
            Ns = extinction(N; protect = combos[i][1], mechanism = combos[i][2])

            for j in eachindex(spread)

                rob = robustness(Ns; threshold = spread[j])
            
                D = DataFrame(
                    protect = combos[i][1],
                    mechanism = combos[i][2],
                    time = networks.time[l],
                    threshold = spread[j],
                    robustness = rob,
                    )
            
                # send to results
                append!(robustness_vals, D)
            end
        end
    end 
end

# robustness curves

layer = data(robustness_vals) * mapping(:threshold, :robustness, 
                                        color = :protect => nonnumeric, 
                                        linestyle = :mechanism,
                                        layout = :time) * (smooth())
fig = draw(layer)

figure = fig.figure