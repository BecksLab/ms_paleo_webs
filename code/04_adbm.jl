using CSV
using DataFrames
using Distributions
using SpeciesInteractionNetworks

include("../lib/adbm/adbm.jl")
include("../lib/internals.jl")

# set seed
import Random
Random.seed!(66)

topology = topo_df();

# get the name of all communities (still need their richness)
size_names = readdir(joinpath("data", "clean", "size"))
size_names = replace.(size_names, ".csv" => "")
trait_names = readdir(joinpath("data", "clean", "trait"))
trait_names = replace.(trait_names, ".csv" => "")

for i in eachindex(size_names)
    
    size_file = size_names[i]
    trait_file = trait_names[i]
    bodymass = DataFrame(
        CSV.File.(
            joinpath("data", "clean", "size", "$size_file.csv"),
        ),
    )
    bodymass = rename(bodymass, :size => :bodymass)
    trait = DataFrame(
        CSV.File.(
            joinpath("data", "clean", "trait", "$trait_file.csv"),
        ),
    )

    df = innerjoin(trait, bodymass, on = :species)

    # create some mock abundance/biomass values
    biomass = rand(Truncated(Normal(0, 1), 0, 10), nrow(df))

    d = model_summary(df, trait_file, "adbm"; bodymass = df.bodymass, biomass = biomass)

    push!(topology, d)

end

# write summaries as .csv
CSV.write("data/processed/topology_adbm.csv", topology)
# write networks as object
save_object("data/processed/networks/adbm_networks.jlds", topology[:, ["id", "model", "network"]])