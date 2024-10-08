using CSV
using DataFrames
using JLD2

include("../lib/pfim/pfim.jl")
include("../lib/internals.jl")

# set seed
import Random
Random.seed!(66)

topology = topo_df();

# get the name of all communities
matrix_names = readdir(joinpath("data", "clean", "trait"))
matrix_names = replace.(matrix_names, ".csv" => "")

for i in eachindex(matrix_names)
    
    file_name = matrix_names[i]
    df = DataFrame(
        CSV.File.(
            joinpath("data", "clean", "trait", "$file_name.csv"),
        ),
    )

    d = model_summary(df, file_name, "pfim")

    push!(topology, d)

end

# write summaries as .csv
CSV.write("data/processed/topology_pfim.csv", topology[:,setdiff(names(topology), ["network"])])
# write networks as object
save_object("data/processed/networks/pfim_networks.jlds", topology[:, ["id", "model", "network"]])