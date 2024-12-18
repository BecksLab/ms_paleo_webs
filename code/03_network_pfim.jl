using CSV
using DataFrames
using JLD2

include("lib/pfim/pfim.jl")
include("lib/internals.jl")

# set seed
import Random
Random.seed!(66)

topology = topo_df();

# get the name of all communities
matrix_names = readdir("../data/clean/trait")
matrix_names = replace.(matrix_names, ".csv" => "")

for i in eachindex(matrix_names)

    file_name = matrix_names[i]
    df = DataFrame(CSV.File.(joinpath("../data/clean/trait", "$file_name.csv"),))

    d = model_summary(df, file_name, "pfim")

    push!(topology, d)

end

# write summaries as .csv
CSV.write(
    "../data/processed/topology/topology_pfim.csv",
    topology[:, setdiff(names(topology), ["network"])],
)
# write networks as object
save_object(
    "../data/processed/networks/pfim_networks.jlds",
    topology[:, ["id", "model", "network"]],
)

# repeat but without downsampling (i.e., metawebs)

topology = topo_df();

for i in eachindex(matrix_names)

    file_name = matrix_names[i]
    df = DataFrame(CSV.File.(joinpath("../data/clean/trait", "$file_name.csv"),))

    d = model_summary(df, file_name, "pfim", downsample = false)

    d[:model] = "pfim_metaweb"
    push!(topology, d)

end

# write summaries as .csv
CSV.write(
    "../data/processed/topology/topology_pfim_metaweb.csv",
    topology[:, setdiff(names(topology), ["network"])],
)
# write networks as object
save_object(
    "../data/processed/networks/pfim_networks_metaweb.jlds",
    topology[:, ["id", "model", "network"]],
)

# now we try with numeric size and not categorical

size_names = readdir("../data/clean/size")
size_names = replace.(size_names, ".csv" => "")

topology = topo_df();

for i in eachindex(matrix_names)

    file_name = matrix_names[i]
    df = DataFrame(CSV.File.(joinpath("../data/clean/trait", "$file_name.csv"),))

    size_file = size_names[i]
    bodymass = DataFrame(CSV.File.(joinpath("../data/clean/size", "$size_file.csv"),))

    # remove categorical sizes
    select!(df, Not([:size]))
    # replace with continuous
    df = innerjoin(df, bodymass, on = :species)

    d = model_summary(df, file_name, "pfim", downsample = false)

    d[:model] = "pfim_size"
    push!(topology, d)

end

# write summaries as .csv
CSV.write(
    "../data/processed/topology/topology_pfim_size.csv",
    topology[:, setdiff(names(topology), ["network"])],
)
# write networks as object
save_object(
    "../data/processed/networks/pfim_networks_size.jlds",
    topology[:, ["id", "model", "network"]],
)

# remove scavengers and parasites

topology = topo_df();

for i in eachindex(matrix_names)

    file_name = matrix_names[i]
    df = DataFrame(CSV.File.(joinpath("../data/clean/trait", "$file_name.csv"),))

    # remove parasites and scavengers
    filter!(row -> row.feeding ∉ ["parasitic", "scavenger"], df)

    d = model_summary(df, file_name, "pfim", downsample = false)

    d[:model] = "pfim_no_scav"
    push!(topology, d)

end

# write summaries as .csv
CSV.write(
    "../data/processed/topology/topology_pfim_no_scav.csv",
    topology[:, setdiff(names(topology), ["network"])],
)
# write networks as object
save_object(
    "../data/processed/networks/pfim_networks_no_scav.jlds",
    topology[:, ["id", "model", "network"]],
)

# repeat but without downsampling (i.e., metawebs) with basal node

topology = topo_df();

for i in eachindex(matrix_names)

    file_name = matrix_names[i]
    df = DataFrame(CSV.File.(joinpath("../data/clean/trait", "$file_name.csv"),))

    _basal = Dict{Symbol,Any}()
    _basal[:species] = "basal"
    _basal[:feeding] = "producer"
    _basal[:motility] = "non_motile"
    _basal[:tiering] = "primary"
    _basal[:size] = "tiny"

    push!(df, _basal)

    d = model_summary(df, file_name, "pfim", downsample = false)

    d[:model] = "pfim_basal"
    push!(topology, d)

end

# write summaries as .csv
CSV.write(
    "../data/processed/topology/topology_pfim_basal.csv",
    topology[:, setdiff(names(topology), ["network"])],
)
# write networks as object
save_object(
    "../data/processed/networks/pfim_networks_basal.jlds",
    topology[:, ["id", "model", "network"]],
)

# trophic species

# get the name of all communities
matrix_names = readdir("../data/clean/trophic")
matrix_names = replace.(matrix_names, ".csv" => "")

topology = topo_df();

for i in eachindex(matrix_names)

    file_name = matrix_names[i]
    df = DataFrame(CSV.File.(joinpath("../data/clean/trophic", "$file_name.csv"),))

    d = model_summary(df, file_name, "pfim", downsample = false)

    d[:model] = "pfim_trophic"
    push!(topology, d)

end

# write summaries as .csv
CSV.write(
    "../data/processed/topology/topology_pfim_trophic.csv",
    topology[:, setdiff(names(topology), ["network"])],
)
# write networks as object
save_object(
    "../data/processed/networks/pfim_networks_trophic.jlds",
    topology[:, ["id", "model", "network"]],
)