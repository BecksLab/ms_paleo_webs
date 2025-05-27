using Combinatorics
using CSV
using DataFrames
using JLD2
using SpeciesInteractionNetworks

include("lib/internals.jl")

β_div = DataFrame(
    left = String[],
    right = String[],
    β_int = Int64[],
    β_spp = Int64[],
    links_left = Int64[],
    links_right = Int64[],
)

# same comment RE elegance
networks = load_object("data/processed/networks.jlds")

# this can for sure be done more elegantly but brute force for now...
# filter!(!=("pfim_networks_basal.jlds"), matrix_names)

networks.uniq_id = map(join,zip(networks.model,networks.time,networks.n_rep))

# for each combination of datasets
for (x, y) in combinations(unique(networks[:, :uniq_id]), 2)
    U = networks[occursin.(x, networks.uniq_id), :network][1]
    V = networks[occursin.(y, networks.uniq_id), :network][1]
    β_spp = SpeciesInteractionNetworks.betadiversity(βS, U, V)
    β_int = SpeciesInteractionNetworks.betadiversity(βWN, U, V)
    # collate
    b = Dict{Symbol,Any}()
    b[:left] = x
    b[:right] = y
    b[:β_int] = β_int.shared
    b[:β_spp] = β_spp.shared
    b[:links_left] = links(U)
    b[:links_right] = links(V)
    # send to results
    push!(β_div, b)
end

# write outputs as .csv
CSV.write("data/processed/beta_div.csv", β_div)
