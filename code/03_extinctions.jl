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
traits = readdir("data/raw")
# select only species datasets
traits = traits[occursin.(r"^.*Guilds.*$", traits)]

# import networks
networks = load_object("data/processed/networks.jlds")
# select only pre extinction networks
subset!(networks, :time => ByRow(x -> x == "G1"))

# specify hierarchies
hierarchies = [
    # trait class
    [:size, :motility, :tiering, :calcification],
    # hierarchy for trait class
    [
        ["tiny", "small", "medium", "large", "very_large", "gigantic"],
        ["non_motile", "motile"],
        ["suspension", "herbivore", "detritivore", "predator"],
        ["heavy", "moderate", "light"],
    ],
]

# data frame for results
extinction_results = DataFrame(
    model = String[],
    extinction_mechanism = Any[],
    n_rep = Int64[],
    extinction_seq = Any[],
);

# Extinction sequence

df = CSV.read("data/raw/G1_Guilds.csv", DataFrame)

select!(
    df,
    [:Guild, :motility_simple, :tiering_simple, :feeding_simple, :size, :calcification],
)

rename!(df, :Guild => :species)
rename!(df, :motility_simple => :motility)
rename!(df, :tiering_simple => :tiering)
rename!(df, :feeding_simple => :feeding)

for j = 1:nrow(networks)

    # select correct network
    N = networks.network[j]

    N = add_basal(N)

    # random extinction
    extinction_series = extinction(N; protect = :basal)

    # send results to data frame
    d = Dict{Symbol,Any}(
        :model => networks.model[j],
        :extinction_mechanism => "random",
        :n_rep => networks.n_rep[j],
        :extinction_seq => extinction_series,
    )

    push!(extinction_results, d)

    # need to simulate for both 'orders' of traits/hierarchies
    # where true is descending and false is ascending
    for descending in [true, false]

        # numeric extinctions
        for k in ["generality", "vulnerability"]

            extinction_series = extinction(N, k, descending)

            # send results to data frame
            d = Dict{Symbol,Any}(
                :model => networks.model[j],
                :extinction_mechanism => join(
                    [string(k), ifelse(descending, "descending", "ascending")],
                    "_",
                ),
                :n_rep => networks.n_rep[j],
                :extinction_seq => extinction_series,
            )

            push!(extinction_results, d)

        end

        # categorical extinctions
        # note we can't do this with the niche model...
        if networks.model[j] ∉ ["niche"]
            for k in eachindex(hierarchies[1])

                # select the correct traits matrix depending on nodes (species vs trophic)
                trait_data = df[:, [:species, hierarchies[1][k]]]
                rename!(trait_data, hierarchies[1][k] => :trait)

                # only include species that have the desired traits
                filter!(:trait => x -> x ∈ hierarchies[2][k], trait_data)

                extinction_list = extinctionsequence(
                    hierarchies[2][k],
                    trait_data;
                    descending = descending,
                )
                # generate extinction sequence
                extinction_series = extinction(N, extinction_list; protect = :basal)

                # send results to data frame
                d = Dict{Symbol,Any}(
                    :model => networks.model[j],
                    :extinction_mechanism => join(
                        [
                            string(hierarchies[1][k]),
                            ifelse(descending, "descending", "ascending"),
                        ],
                        "_",
                    ),
                    :n_rep => networks.n_rep[j],
                    :extinction_seq => extinction_series,
                )

                push!(extinction_results, d)
            end
        end
    end
end

# write networks as object
save_object("data/processed/extinction_seq.jlds", extinction_results)
