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
matrix_names = readdir("../data/clean/trait_maximal")
matrix_names = replace.(matrix_names, ".csv" => "")

# maximal feeding rules
feeding_rules = DataFrame(CSV.File("../data/clean/feeding_rules/feeding_rules_maximal.csv"))

for i in eachindex(matrix_names)

    file_name = matrix_names[i]
    df = DataFrame(CSV.File.(joinpath("../data/clean/trait_maximal", "$file_name.csv"),))

    # remove # primary species
    filter!(row -> row.feeding ∉ ["parasitic", "scavenger", "primary_feeding"], df)

    d = model_summary(
        df,
        file_name,
        "pfim";
        feeding_rules = feeding_rules,
        downsample = false,
    )

    d[:model] = "pfim_maximal"
    push!(topology, d)

end

# write summaries as .csv
CSV.write(
    "../data/processed/topology/topology_pfim_maximal.csv",
    topology[:, setdiff(names(topology), ["network"])],
)
# write networks as object
save_object(
    "../data/processed/networks/pfim_maximal_networks.jlds",
    topology[:, ["id", "model", "network"]],
)

# minimum feeding rules

topology = topo_df();

feeding_rules = DataFrame(CSV.File("../data/clean/feeding_rules/feeding_rules_minimum.csv"))

for i in eachindex(matrix_names)

    file_name = matrix_names[i]
    df = DataFrame(CSV.File.(joinpath("../data/clean/trait_minimum", "$file_name.csv"),))

    # remove # primary species
    filter!(row -> row.feeding ∉ ["parasitic", "scavenger", "primary_feeding"], df)

    d = model_summary(
        df,
        file_name,
        "pfim";
        feeding_rules = feeding_rules,
        downsample = false,
    )

    d[:model] = "pfim_minimum"
    push!(topology, d)

end

# write summaries as .csv
CSV.write(
    "../data/processed/topology/topology_pfim_minimum.csv",
    topology[:, setdiff(names(topology), ["network"])],
)
# write networks as object
save_object(
    "../data/processed/networks/pfim_minimum_networks.jlds",
    topology[:, ["id", "model", "network"]],
)

# numerical size

size_names = readdir("../data/clean/size")
size_names = replace.(size_names, ".csv" => "")

topology = topo_df();

for i in eachindex(matrix_names)

    file_name = matrix_names[i]
    df = DataFrame(CSV.File.(joinpath("../data/clean/trait_minimum", "$file_name.csv"),))

    # remove primary node
    filter!(row -> row.feeding ∉ ["parasitic", "scavenger", "primary_feeding"], df)

    size_file = size_names[i]
    bodymass = DataFrame(CSV.File.(joinpath("../data/clean/size", "$size_file.csv"),))

    select!(df, Not([:size]))
    # replace with continuous
    df = innerjoin(df, bodymass, on = :species)

    d = model_summary(
        df,
        file_name,
        "pfim";
        feeding_rules = feeding_rules,
        downsample = false,
    )

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

feeding_rules = DataFrame(CSV.File("../data/clean/feeding_rules/feeding_rules_noscavs.csv"))

for i in eachindex(matrix_names)

    file_name = matrix_names[i]
    df = DataFrame(CSV.File.(joinpath("../data/clean/trait_minimum", "$file_name.csv"),))

    # keep parasites and scavengers
    filter!(row -> row.feeding ∉ ["primary_feeding"], df)

    d = model_summary(
        df,
        file_name,
        "pfim";
        feeding_rules = feeding_rules,
        downsample = false,
    )

    d[:model] = "pfim_with_scav"
    push!(topology, d)

end

# write summaries as .csv
CSV.write(
    "../data/processed/topology/topology_pfim_with_scav.csv",
    topology[:, setdiff(names(topology), ["network"])],
)
# write networks as object
save_object(
    "../data/processed/networks/pfim_networks_with_scav.jlds",
    topology[:, ["id", "model", "network"]],
)

# keep primary (basal) node

topology = topo_df();

# get the name of all communities
matrix_names = readdir("../data/clean/trait_minimum")
matrix_names = replace.(matrix_names, ".csv" => "")

feeding_rules = DataFrame(CSV.File("../data/clean/feeding_rules/feeding_rules_minimum.csv"))

for i in eachindex(matrix_names)

    file_name = matrix_names[i]
    df = DataFrame(CSV.File.(joinpath("../data/clean/trait_minimum", "$file_name.csv"),))

    # remove parasites and scavengers
    filter!(row -> row.feeding ∉ ["parasitic", "scavenger"], df)

    d = model_summary(
        df,
        file_name,
        "pfim";
        feeding_rules = feeding_rules,
        downsample = false,
    )

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

# trophic species - minimum rules

# get the name of all communities
matrix_names = readdir("../data/clean/trophic_minimum")
matrix_names = replace.(matrix_names, ".csv" => "")

topology = topo_df();

for i in eachindex(matrix_names)

    file_name = matrix_names[i]
    df = DataFrame(CSV.File.(joinpath("../data/clean/trophic_minimum", "$file_name.csv"),))

    # remove parasites and scavengers
    filter!(row -> row.feeding ∉ ["parasitic", "scavenger", "primary_feeding"], df)

    d = model_summary(
        df,
        file_name,
        "pfim";
        feeding_rules = feeding_rules,
        downsample = false,
    )

    d[:model] = "pfim_trophic_minimum"
    push!(topology, d)

end

# write summaries as .csv
CSV.write(
    "../data/processed/topology/topology_pfim_trophic_minimum.csv",
    topology[:, setdiff(names(topology), ["network"])],
)
# write networks as object
save_object(
    "../data/processed/networks/pfim_networks_trophic_minimum.jlds",
    topology[:, ["id", "model", "network"]],
)

# trophic species - maximal rules

# get the name of all communities
matrix_names = readdir("../data/clean/trophic_maximal")
matrix_names = replace.(matrix_names, ".csv" => "")

# maximal feeding rules
feeding_rules = DataFrame(CSV.File("../data/clean/feeding_rules/feeding_rules_maximal.csv"))

topology = topo_df();

for i in eachindex(matrix_names)

    file_name = matrix_names[i]
    df = DataFrame(CSV.File.(joinpath("../data/clean/trophic_maximal", "$file_name.csv"),))

    # remove parasites and scavengers
    filter!(row -> row.feeding ∉ ["parasitic", "scavenger", "primary_feeding"], df)

    d = model_summary(
        df,
        file_name,
        "pfim";
        feeding_rules = feeding_rules,
        downsample = false,
    )

    d[:model] = "pfim_trophic_maximal"
    push!(topology, d)

end

# write summaries as .csv
CSV.write(
    "../data/processed/topology/topology_pfim_trophic_maximal.csv",
    topology[:, setdiff(names(topology), ["network"])],
)
# write networks as object
save_object(
    "../data/processed/networks/pfim_networks_trophic_maximal.jlds",
    topology[:, ["id", "model", "network"]],
)

# downsample minimum

topology = topo_df();

feeding_rules = DataFrame(CSV.File("../data/clean/feeding_rules/feeding_rules_minimum.csv"))

for i in eachindex(matrix_names)

    file_name = matrix_names[i]
    df = DataFrame(CSV.File.(joinpath("../data/clean/trait_minimum", "$file_name.csv"),))

    # remove parasites and scavengers
    filter!(row -> row.feeding ∉ ["parasitic", "scavenger", "primary_feeding"], df)

    d = model_summary(
        df,
        file_name,
        "pfim";
        feeding_rules = feeding_rules,
        downsample = true,
    )

    d[:model] = "pfim_minimum_downsample"
    push!(topology, d)

end

# write summaries as .csv
CSV.write(
    "../data/processed/topology/topology_pfim_minimum_downsample.csv",
    topology[:, setdiff(names(topology), ["network"])],
)
# write networks as object
save_object(
    "../data/processed/networks/pfim_minimum_downsample_networks.jlds",
    topology[:, ["id", "model", "network"]],
)

# downsample maximal

topology = topo_df();

feeding_rules = DataFrame(CSV.File("../data/clean/feeding_rules/feeding_rules_maximal.csv"))

for i in eachindex(matrix_names)

    file_name = matrix_names[i]
    df = DataFrame(CSV.File.(joinpath("../data/clean/trait_maximal", "$file_name.csv"),))

    # remove parasites and scavengers
    filter!(row -> row.feeding ∉ ["parasitic", "scavenger", "primary_feeding"], df)

    d = model_summary(
        df,
        file_name,
        "pfim";
        feeding_rules = feeding_rules,
        downsample = true,
    )

    d[:model] = "pfim_maximal_downsample"
    push!(topology, d)

end

# write summaries as .csv
CSV.write(
    "../data/processed/topology/topology_pfim_maximal_downsample.csv",
    topology[:, setdiff(names(topology), ["network"])],
)
# write networks as object
save_object(
    "../data/processed/networks/pfim_maximal_downsample_networks.jlds",
    topology[:, ["id", "model", "network"]],
)
