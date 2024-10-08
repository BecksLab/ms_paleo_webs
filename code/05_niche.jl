using CSV
using DataFrames
using SpeciesInteractionNetworks

include("../lib/internals.jl")

# specify connectance
Co = 0.1

# set seed
import Random
Random.seed!(66)

topology = topo_df();

# get the name of all communities (still need their richness)
matrix_names = readdir(joinpath("data", "clean", "trait"))
matrix_names = replace.(matrix_names, ".csv" => "")

for i in eachindex(matrix_names)
    
    file_name = matrix_names[i]
    df = DataFrame(
        CSV.File.(
            joinpath("data", "clean", "trait", "$file_name.csv"),
        ),
    )

    d = model_summary(df, file_name, "niche"; connectance = Co)

    push!(topology, d)

end

# write summaries as .csv
CSV.write("data/processed/topology_niche.csv", topology)
# write networks as object
save_object("data/processed/networks/niche_networks.jlds", topology[:, ["id", "model", "network"]])