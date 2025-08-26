# Purpose: Compare the Julia PFIM metawebs with Dunhill 2024 webs
# Author: Tanya Strydom
# Copyright/License: MIT
# Contact: t.strydom@sheffield.ac.uk
# Date: 2025-08-19

# libraries
using AlgebraOfGraphics
using CairoMakie
using CSV
using DataFrames
using Extinctions
using JLD2
using pfim
using SpeciesInteractionNetworks
using StatsBase

# set seed
import Random
Random.seed!(66)

# internal functions
include("lib/internals.jl")

# get the name of all communities
matrix_names = readdir("../data/dunhill")
# select only species datasets
matrix_names = matrix_names[occursin.(r"^.*edgelist.*$", matrix_names)]

# df to store networks
networks = DataFrame(time = Any[], network = Any[]);

# build networks from Dunhill edgelists
for i in eachindex(matrix_names)

    file_name = matrix_names[i]
    # get relevant info from slug
    str_cats = split(file_name, r"_")
    # get interactions
    ints = DataFrame(CSV.File.(joinpath("../data/dunhill/", "$file_name")))

    # build network
    S = unique(vcat(ints.resource, ints.consumer))
    nodes = Unipartite(Symbol.(S))
    edgs = Binary(zeros(Bool, (length(S), length(S))))

    N = SpeciesInteractionNetwork(nodes, edgs)

    for i = 1:nrow(ints)
        interaction = (Symbol(ints.consumer[i]), Symbol(ints.resource[i]))
        N[interaction...] = true
    end

    # push to df
    d = Dict{Symbol,Any}(
        :time => str_cats[1],
        :network => N,
    )
    push!(networks, d)
    
end

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

combos = [[:none; :cascade],
        [:none; :secondary],
        [:basal; :cascade],
        [:basal; :secondary]]
        
ext_reps = 500

for h in 1:nrow(networks)

    # remove cannibals
    N = networks.network[h]

    for i = eachindex(combos)
    
        for l = 1:ext_reps

            spp = StatsBase.shuffle(species(N))

            if combos[i][1] == :basal
                filter!(x -> x != Symbol("BASAL NODE"),spp)
            end

            # pre-defined extinction sequence
            Ns = extinction(N, spp; protect = :none, mechanism = combos[i][2])
            
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

save("../figures/robustness_compare.png", figure)

# write summaries as .csv
CSV.write(
    "../data/processed/robustness_test.csv",
    robustness_vals,
)

# write summaries as .csv
CSV.write(
    "../data/processed/auc_test.csv",
    auc,
)

# robustness but we call a new extinction series for every threshold value

robustness_vals = DataFrame(
    protect = Any[],
    mechanism = Any[],
    time = Any[],
    threshold = Any[],
    robustness = Any[],
);

for h in 1:nrow(networks)

    # remove cannibals
    N = networks.network[h]

    for i = eachindex(combos)
    
        for l = 1:ext_reps

            for j in eachindex(spread)

            spp = StatsBase.shuffle(species(N))

            if combos[i][1] == :basal
                filter!(x -> x != Symbol("BASAL NODE"),spp)
            end

            # pre-defined extinction sequence
            Ns = extinction(N, spp; protect = :none, mechanism = combos[i][2])

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

# write summaries as .csv
CSV.write(
    "../data/processed/robustness_many_to_many.csv",
    robustness_vals,
)