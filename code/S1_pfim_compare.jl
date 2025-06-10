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

for j in eachindex(dunhill_comms)

    # read in edgelist
    file_name = dunhill_comms[j]
    ints = CSV.read(joinpath("data/dunhill/", "$file_name"), DataFrame)

    # build DUNHILL network
    S = unique(vcat(ints.resource, ints.consumer))

    nodes = Unipartite(Symbol.(S))
    edgs = Binary(zeros(Bool, (length(S), length(S))))
    N_dunhill = SpeciesInteractionNetwork(nodes, edgs)

    for i in 1:nrow(ints)

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
