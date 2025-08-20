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

# set seed
import Random
Random.seed!(66)

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

auc = DataFrame(
    protect = Any[],
    mechanism = Any[],
    time = Any[],
    AUC = Any[],
);

# represents the % of species that have gone extinct (primary and secondary)
spread = collect(1:1:99)

# first we create `ext_reps` number of extinction series for each network
# note we do not protect the basal node

combos = [[:none; :cascade],
        [:none; :secondary],
        [:basal; :cascade],
        [:basal; :secondary]]
        
ext_reps = 500

for h in 1:4

    # remove cannibals
    N = remove_cannibals(networks.network[h])

    for i = eachindex(combos)
    
        for l = 1:ext_reps
            
            Ns = extinction(N; protect = combos[i][1], mechanism = combos[i][2])

             D = DataFrame(
                    protect = combos[i][1],
                    mechanism = combos[i][2],
                    time = h,
                    AUC = resilience(Ns)
                    )

            # send to results
            append!(auc, D)

            for j in eachindex(spread)

                rob = robustness(Ns; threshold = spread[j])
            
                D = DataFrame(
                    protect = combos[i][1],
                    mechanism = combos[i][2],
                    time = networks.time[h],
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
                                        color = :time, 
                                        linestyle = :mechanism,
                                        layout = :protect) * (smooth())
fig = draw(layer)

figure = fig.figure

save("figures/robustness_compare.png", figure)

# write summaries as .csv
CSV.write(
    "data/processed/robustness_test.csv",
    robustness_vals,
)

# write summaries as .csv
CSV.write(
    "data/processed/auc_test.csv",
    auc,
)