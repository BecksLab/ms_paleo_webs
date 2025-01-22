using CSV
using DataFrames
using Extinctions
using JLD2
using SpeciesInteractionNetworks

include("lib/internals.jl")

# set seed
import Random
Random.seed!(66)

# get the traits data
traits = DataFrame(CSV.File("../data/clean/extinction/traits.csv"))
traits_trophic = DataFrame(CSV.File("../data/clean/extinction/traits_trophic.csv"))

# specify hierarchies
hierarchies = [
    # trait class
    [:size, :motility, :tiering],
    # hierarchy for trait class
    [
        ["tiny", "small", "medium", "large"],
        ["non_motile", "facultative", "slow", "fast"],
        ["infaunal", "epifaunal", "nektonic", "pelagic"],
    ],
]

# Extinction sequence

# for results
extinction_results = topo_df()
# add addition column
extinction_results[!, :time] = String[]
extinction_results[!, :extinction_mechanism] = String[]
extinction_results[!, :links] = Int64[]
# remove unused columns
select!(extinction_results, Not(:network))

matrix_names = readdir("../data/processed/networks")

for i in eachindex(matrix_names)

    # import predicted network for specific model
    file_name = matrix_names[i]
    df = load_object("../data/processed/networks/$file_name")

    # now we loop through the different time periods
    for j in ["pre", "during"]
        # select only the relevant community
        pre_comm = df[occursin.(j, df.id), :]
        # get post extinction richness
        post_rich = richness(df[occursin.("post", df.id), :network][1])

        # random extinction
        extinction_series = extinction(pre_comm.network[1], post_rich)

        # summarise extinction network
        D = _network_summary(extinction_series[end])

        # append additional info
        D[:model] = pre_comm.model[1]
        D[:extinction_mechanism] = "random"
        D[:id] = pre_comm.id[1]
        D[:time] = j

        # send to results
        push!(extinction_results, D)

        # need to simulate for both 'orders' of traits/hierarchies
        # where true is descending and false is ascending
        for descending in [true, false]

            # numeric extinctions
            for k in ["generality", "vulnerability"]

                # turn into function
                f = getfield(Main, Symbol(k))

                extinction_list = extinctionsequence(f(pre_comm.network[1]); descending)

                # generate extinction sequence
                extinction_series =
                    extinction(pre_comm.network[1], extinction_list, post_rich)

                # summarise extinction network
                D = _network_summary(extinction_series[end])

                # append additional info
                D[:model] = pre_comm.model[1]
                D[:extinction_mechanism] =
                    join([k, ifelse(descending, "descending", "ascending")], "_")
                D[:id] = pre_comm.id[1]
                D[:time] = j

                push!(extinction_results, D)
            end

            # categorical extinctions
            # note we can't do this with the niche model...
            if pre_comm.model[1] ∉ ["niche"]
                for k in eachindex(hierarchies[1])

                    # select the correct traits matrix depending on nodes (species vs trophic)
                    if occursin.("trophic", file_name)
                        trait_data = traits_trophic[:, [:species, hierarchies[1][k]]]
                    else
                        trait_data = traits[:, [:species, hierarchies[1][k]]]
                    end

                    rename!(trait_data, hierarchies[1][k] => :trait)

                    # only include species that have the desired traits
                    filter!(:trait => x -> x ∈ hierarchies[2][k], trait_data)

                    extinction_list =
                        extinctionsequence(hierarchies[2][k], trait_data; descending)

                    # generate extinction sequence
                    extinction_series =
                        extinction(pre_comm.network[1], extinction_list, post_rich)

                    # summarise extinction network
                    D = _network_summary(extinction_series[end])

                    # append additional info
                    D[:model] = pre_comm.model[1]
                    D[:extinction_mechanism] = join(
                        [
                            string(hierarchies[1][k]),
                            ifelse(descending, "descending", "ascending"),
                        ],
                        "_",
                    )
                    D[:id] = pre_comm.id[1]
                    D[:time] = j

                    push!(extinction_results, D)
                end
            end
        end
    end
end

# write summaries as .csv
CSV.write(
    "../data/processed/extinctions/extinctions.csv",
    extinction_results[:, setdiff(names(extinction_results), ["network"])],
)
