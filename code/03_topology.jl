using CSV
using DataFrames
using JLD2

include("lib/internals.jl")

# import networks object
networks = load_object("data/processed/networks.jlds")

topology = DataFrame(
    model = String[],
    location = String[],
    time = Any[],
    richness = Int64[],
    links = Int64[],
    connectance = Float64[],
    diameter = Int64[],
    complexity = Float64[],
    distance = Float64[],
    basal = Float64[],
    top = Float64[],
    generality = Float64[],
    vulnerability = Float64[],
    S1 = Float64[],
    S2 = Float64[],
    S4 = Float64[],
    S5 = Float64[],
);

for i in 1:nrow(networks)
    
    d = _network_summary(networks.network[i])

    d[:model] = networks.model[i]
    d[:location] = networks.location[i]
    d[:time] = networks.time[i]

    push!(topology, d)
end

# write summaries as .csv
CSV.write(
    "data/processed/topology.csv",
    topology,
)