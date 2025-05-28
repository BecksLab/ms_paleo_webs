using Combinatorics
using CSV
using DataFrames
using JLD2
using SpeciesInteractionNetworks

include("lib/internals.jl")

β_div = DataFrame(
    time = Any[],
    left = String[],
    right = String[],
    β_int = Int64[],
    β_spp = Int64[],
    links_left = Int64[],
    links_right = Int64[],
)

# same comment RE elegance
networks = load_object("data/processed/networks.jlds")

# get the different communities
time_groups = unique(networks.time)

for i in eachindex(time_groups)

    df = filter(:time => x -> x == time_groups[i], networks)
    # make unique ID for making combos
    df.uniq_id = map(join,zip(df.model, df.n_rep))

    # for each combination of datasets
    for (x, y) in combinations(unique(df[:, :uniq_id]), 2)
        U = df[occursin.(x, df.uniq_id), :network][1]
        V = df[occursin.(y, df.uniq_id), :network][1]
        β_spp = SpeciesInteractionNetworks.betadiversity(βS, U, V)
        β_int = SpeciesInteractionNetworks.betadiversity(βWN, U, V)
        # collate
        b = Dict{Symbol,Any}()
        b[:time] = time_groups[i]
        b[:left] = x
        b[:right] = y
        b[:β_int] = β_int.shared
        b[:β_spp] = β_spp.shared
        b[:links_left] = links(U)
        b[:links_right] = links(V)
        # send to results
        push!(β_div, b)
    end   
end

# write outputs as .csv
CSV.write("data/processed/beta_div.csv", β_div)
