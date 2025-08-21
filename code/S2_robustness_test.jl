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
using pfim
using SpeciesInteractionNetworks

# set seed
import Random
Random.seed!(66)

# internal functions
include("lib/internals.jl")

# get the name of all communities
matrix_names = readdir("data/raw")
# select only species datasets
matrix_names = matrix_names[occursin.(r"^.*Guilds.*$", matrix_names)]

# feeding rules
feeding_rules = DataFrame(CSV.File("data/raw/feeding_rules.csv"))

# size classes (for creating continuous body sizes)
size_classes = DataFrame(CSV.File("data/raw/size_classes.csv"))

# df to store networks
networks = DataFrame(time = Any[], network = Any[]);

for i in eachindex(matrix_names)
    
    file_name = matrix_names[i]
    # get relevant info from slug
    str_cats = split(file_name, r"_")
    # import data frame
    df = DataFrame(CSV.File.(joinpath("data/raw/", "$file_name")))
    select!(df, [:Guild, :motility, :tiering, :feeding, :size])
    rename!(df, :Guild => :species)
    N = pfim.PFIM(df, feeding_rules; downsample = false)
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

# first we create `ext_reps` number of extinction series for each network
# note we do not protect the basal node

combos = [[:none; :cascade],
        [:none; :secondary],
        [:basal; :cascade],
        [:basal; :secondary]]
        
ext_reps = 100

for h in 1:nrow(networks)

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