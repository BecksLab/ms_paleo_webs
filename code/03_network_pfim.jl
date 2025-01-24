using CSV
using DataFrames
using JLD2

include("lib/pfim/pfim.jl")
include("lib/internals.jl")

# set seed
import Random
Random.seed!(66)

# get the name of all communities
matrix_names = readdir("../data/clean/trait_maximal")
matrix_names = replace.(matrix_names, ".csv" => "")

# feeding rules
max = DataFrame(CSV.File("../data/clean/feeding_rules/feeding_rules_maximal.csv"))
min = DataFrame(CSV.File("../data/clean/feeding_rules/feeding_rules_minimum.csv"))

rules = [
    "maximal" "minimum"
    max min
]

# diff datasets
model_names = ["pfim", "pfim_size", "pfim_basal", "pfim_trophic", "pfim_with_scav"]

for j in eachindex(model_names)

    # run PFIM with and without downsampling with both maximal and minimum rules
    for k in axes(rules, 1)
        feeding_rules = rules[2, k]
        feeding = rules[1, k]

        topology = topo_df()

        for i in eachindex(matrix_names)
            # get data frame
            file_name = matrix_names[i]
            df = DataFrame(
                CSV.File.(joinpath("../data/clean/trait_maximal", "$file_name.csv"),),
            )

            model = model_names[j]

            # add/remove nodes based on model/dataset
            if model == "pfim_with_scav"
                # remove only basal node
                filter!(row -> row.feeding ∉ ["primary_feeding"], df)
            elseif model == "pfim_basal"
                # remove parasites and scavengers
                filter!(row -> row.feeding ∉ ["parasitic", "scavenger"], df)
            else
                # remove scavenger, parasitic, primary species
                filter!(
                    row -> row.feeding ∉ ["parasitic", "scavenger", "primary_feeding"],
                    df,
                )
            end

            for downsample ∈ [true, false]

                d = model_summary(
                    df,
                    file_name,
                    "pfim";
                    feeding_rules = feeding_rules,
                    downsample = downsample,
                )

                if downsample == true
                    d[:model] = join([model, feeding, "downsample"], "_")
                    push!(topology, d)
                else
                    d[:model] = join([model, feeding], "_")
                    push!(topology, d)
                end
            end
        end
        # write summaries as .csv
        CSV.write(
            join([
                "../data/processed/topology/",
                join(["topology", "$model", "$feeding"], "_"),
                ".csv",
            ]),
            topology[:, setdiff(names(topology), ["network"])],
        )
        # write networks as object
        save_object(
            joinpath(
                "../data/processed/networks/",
                join(["$model", "$feeding", "_networks.jlds"], "_"),
            ),
            topology[:, ["id", "model", "network"]],
        )
    end
end
