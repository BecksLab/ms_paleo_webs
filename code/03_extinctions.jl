using CSV
using DataFrames
using Extinctions
using JLD2
using ProgressMeter
using SpeciesInteractionNetworks

include("lib/internals.jl")

# Set seed for reproducible extinction sequences
import Random
Random.seed!(66)

# Data Preparation: Identify trait datasets
traits = readdir("../data/raw")
traits = traits[occursin.(r"^.*Guilds.*$", traits)]

# Load previously generated networks and isolate the "G1" time period for baseline testing
networks = load_object("../data/processed/networks.jlds")
subset!(networks, :time => ByRow(x -> x == "G1"))

# Define Trait Hierarchies for Categorical Extinctions
# This defines the order in which species are removed based on their ecological traits
hierarchies = [
    # Trait classes available in the data
    [:size, :motility, :tiering, :calcification],
    # Specific hierarchies (e.g., small to large, or infaunal to pelagic)
    [
        ["tiny", "small", "medium", "large", "very_large", "gigantic"],
        ["non_motile", "motile"],
        ["infaunal", "epifaunal", "pelagic"],
        ["heavy", "moderate", "light"],
    ],
]

# Initialize storage for extinction results
extinction_results = DataFrame(
    model = String[],
    extinction_mechanism = Any[],
    n_rep = Int64[],
    extinction_seq = Any[],
);

# Load and clean the G1 trait data specifically for mapping node attributes
df = CSV.read("../data/raw/G1_Guilds.csv", DataFrame)
select!(
    df,
    [:Guild, :motility_simple, :tiering_simple, :feeding_simple, :size, :calcification],
)
rename!(df, :Guild => :species)
rename!(df, :motility_simple => :motility)
rename!(df, :tiering_simple => :tiering)
rename!(df, :feeding_simple => :feeding)

# Exclude BASAL_NODE to treat it as a protected source of energy/refuge
filter!(:species => x -> x != "BASAL_NODE", df)

# Number of extinction simulations per network
ext_reps = 50
# get each model - we will save each extinction cluster by model to keep file size low
models = unique(networks.model)

@showprogress "Running extinctions" for model_name in models

    # subset by model
    model_networks = subset(networks, :model => ByRow(==(model_name)))

    # fresh df
    extinction_results = DataFrame(
        model = String[],
        extinction_mechanism = String[],
        n_rep = Int[],
        extinction_seq = Any[],
    )

    for j in 1:nrow(model_networks)
        for l in 1:ext_reps

            N = add_basal(model_networks.network[j])

            # Random extinction
            extinction_series = extinction(N; protect = :basal)
            push!(extinction_results, (
                model = model_name,
                extinction_mechanism = "random",
                n_rep = model_networks.n_rep[j],
                extinction_seq = extinction_series,
            ))

            # Topological + trait extinctions
            for descending in (true, false)

                for k in ("generality", "vulnerability")
                    extinction_series = extinction(N, k, descending; protect = :basal)

                    push!(extinction_results, (
                        model = model_name,
                        extinction_mechanism = "$(k)_$(descending ? "descending" : "ascending")",
                        n_rep = model_networks.n_rep[j],
                        extinction_seq = extinction_series,
                    ))
                end

                for k in eachindex(hierarchies[1])

                    trait_data = df[:, [:species, hierarchies[1][k]]]
                    rename!(trait_data, hierarchies[1][k] => :trait)
                    filter!(:trait => x -> x ∈ hierarchies[2][k], trait_data)

                    extinction_list = extinctionsequence(
                        hierarchies[2][k],
                        trait_data;
                        descending = descending,
                    )

                    extinction_series = extinction(N, extinction_list; protect = :basal)

                    push!(extinction_results, (
                        model = model_name,
                        extinction_mechanism = "$(hierarchies[1][k])_$(descending ? "descending" : "ascending")",
                        n_rep = model_networks.n_rep[j],
                        extinction_seq = extinction_series,
                    ))
                end
            end
        end
    end

    # Save one file for this model
    save_object(
        "../data/processed/extinction_seq/$(model_name).jld2",
        extinction_results,
    )

    # Release memory before the next model
    extinction_results = nothing
    GC.gc()
end