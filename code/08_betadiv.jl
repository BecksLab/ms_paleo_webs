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

matrix_names = readdir("../data/processed/networks")

# this can for sure be done more elegantly but brute force for now...
filter!(!=("pfim_networks_basal.jlds"), matrix_names)
filter!(!=("pfim_networks_trophic.jlds"), matrix_names)
filter!(!=("pfim_networks_size.jlds"), matrix_names)
filter!(!=("pfim_networks_metaweb.jlds"), matrix_names)
filter!(!=("pfim_networks_no_scav.jlds"), matrix_names)

# same comment RE elegance
networks = load_object(joinpath("../data/processed/networks/", matrix_names[1]))

for i in 2:5
    append!(networks, load_object(joinpath("../data/processed/networks/", matrix_names[i])))
end

for time ∈ ["pre", "during", "post"]

    df = networks[occursin.(time, networks.id), :]

    for i in eachindex(matrix_names)

        # for each combination of datasets
        for (x, y) in combinations(unique(networks[:,:model]), 2)
            
            U = df[occursin.(x, df.model), :network][1]
            V = df[occursin.(y, df.model), :network][1]
            β_spp = SpeciesInteractionNetworks.betadiversity(βS,U,V)
            β_int = SpeciesInteractionNetworks.betadiversity(βWN,U,V)
    
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
    end
end

# write outputs as .csv
CSV.write(
    "../data/processed/beta_div/beta_div.csv",
    β_div,
)
