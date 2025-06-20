# Purpose: Compare the Julia PFIM metawebs with Dunhill 2024 webs
# Author: Tanya Strydom
# Copyright/License: MIT
# Contact: t.strydom@sheffield.ac.uk
# Date: 2025-06-10

# libraries
using CSV
using DataFrames
using pfim
using SpeciesInteractionNetworks

# internal functions
include("lib/internals.jl")

# get the name of all Dunhill communities
dunhill_comms = readdir("data/dunhill")
dunhill_comms = dunhill_comms[occursin.(r"^.*csv.*$", dunhill_comms)]

# get the name of all communities
matrix_names = readdir("data/raw")
# select only species datasets
matrix_names = matrix_names[occursin.(r"^.*Guilds.*$", matrix_names)]

# feeding rules
feeding_rules = DataFrame(CSV.File("data/raw/feeding_rules.csv"))

for j = 1:2

    # read in edgelist
    file_name = dunhill_comms[j]
    ints = CSV.read(joinpath("data/dunhill/", "$file_name"), DataFrame)

    # build DUNHILL network
    S = unique(vcat(ints.resource, ints.consumer))

    nodes = Unipartite(Symbol.(S))
    edgs = Binary(zeros(Bool, (length(S), length(S))))
    N_dunhill = SpeciesInteractionNetwork(nodes, edgs)

    for i = 1:nrow(ints)

        interaction = (Symbol(ints.consumer[i]), Symbol(ints.resource[i]))

        N_dunhill[interaction...] = true
    end

    d_dunhill = _network_summary(N_dunhill)

    # build JULIA network
    # import data frame
    df = DataFrame(CSV.File.(joinpath("data/raw/", matrix_names[j])))
    select!(df, [:Guild, :motility, :tiering, :feeding, :size])
    rename!(df, :Guild => :species)

    N_julia = pfim.PFIM(df, feeding_rules; downsample = false)

    d_julia = _network_summary(N_julia)

end

# set seed
import Random
Random.seed!(66)

# get the name of all communities
matrix_names = readdir("data/raw")
# select only species datasets
matrix_names = matrix_names[occursin.(r"^.*Guilds.*$", matrix_names)]

# feeding rules
feeding_rules = DataFrame(CSV.File("data/raw/feeding_rules.csv"))

# number of network reps
n_reps = 100
y_val = collect(1.0:0.5:50.0)

topology = DataFrame(
    time = Any[],
    n_rep = Any[],
    y_val = Any[],
    richness = Int64[],
    connectance = Float64[],
    diameter = Int64[],
    complexity = Float64[],
    trophic_level = Float64[],
    distance = Float64[],
    generality = Float64[],
    vulnerability = Float64[],
    redundancy = Float64[],
    S1 = Float64[],
    S2 = Float64[],
    S4 = Float64[],
    S5 = Float64[],
);

for j = 1:n_reps

    for i in eachindex(matrix_names)

        file_name = matrix_names[i]
        # get relevant info from slug
        str_cats = split(file_name, r"_")

        # import data frame
        df = DataFrame(CSV.File.(joinpath("data/raw/", "$file_name")))
        select!(df, [:Guild, :motility, :tiering, :feeding, :size])
        rename!(df, :Guild => :species)

        # remove BASAL_NODE for now...
        filter!(:species => x -> x != "BASAL_NODE", df)

        for k in eachindex(y_val)

            N = pfim.PFIM(df, feeding_rules; y = y_val[k], downsample = true)

            d = _network_summary(N)

            d[:time] = str_cats[1]
            d[:n_rep] = j
            d[:y_val] = y_val[k]

            push!(topology, d)
        end
    end
end

# write summaries as .csv
CSV.write("data/dunhill/topology.csv", topology)
