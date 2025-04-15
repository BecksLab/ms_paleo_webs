using CSV
using DataFrames
using Extinctions
using JLD2
using SpeciesInteractionNetworks

include("lib/internals.jl")

# set seed
import Random
Random.seed!(66)

# import networks object
extinctions = load_object("data/processed/extinction_seq.jlds")

# find richness of post extinction community
df = CSV.read("data/raw/G2_Guilds.csv", DataFrame)
post_rich = nrow(df)

topology = DataFrame(
    model = String[],
    extinction_mechanism = Any[],
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
    redundancy = Float64[],
    S1 = Float64[],
    S2 = Float64[],
    S4 = Float64[],
    S5 = Float64[],
    resilience = Float64[],
    r50 = Any[],
    rep = Int64[],
);

for i = 1:nrow(extinctions)

    _ext_seq = extinctions.extinction_seq[i]

    if length(_ext_seq) > 1

        net_ind = findfirst(x -> x < post_rich, richness.(_ext_seq))

        # this is to catch networks that 
        if typeof(net_ind) == Int64 && richness(_ext_seq[net_ind]) > 0

            d = _network_summary(_ext_seq[net_ind])

            d[:model] = extinctions.model[i]
            d[:extinction_mechanism] = extinctions.extinction_mechanism[i]
            d[:resilience] = resilience(_ext_seq)
            d[:r50] = robustness(_ext_seq)
            d[:rep] = i

            push!(topology, d)
        end
    end
end

# write summaries as .csv
CSV.write("data/processed/extinction_topology.csv", topology)
