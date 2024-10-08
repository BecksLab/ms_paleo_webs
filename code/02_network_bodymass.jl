using CSV
using DataFrames

include("../lib/bodymass/bodymass.jl")
include("../lib/internals.jl")

# set seed
import Random
Random.seed!(66)

topology = topo_df();

# get the name of all communities
matrix_names = readdir(joinpath("data", "clean", "size"))
matrix_names = replace.(matrix_names, ".csv" => "")

for i in eachindex(matrix_names)
    
    file_name = matrix_names[i]
    df = DataFrame(
        CSV.File.(
            joinpath("data", "clean", "size", "$file_name.csv"),
        ),
    )

    d = model_summary(df, file_name, "bodymassratio"; bodymass = df.size)

    push!(topology, d)

end

# write summaries as .csv
CSV.write("data/processed/topology_bodymassratio.csv", topology)
# write networks as object
save_object("data/processed/networks/bodymassratio_networks.jlds", topology[:, ["id", "model", "network"]])