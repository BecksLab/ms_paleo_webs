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

# Main Extinction Loop
@showprogress "Running extinctions" for j = 1:nrow(networks)

    for l = 1:ext_reps

        # Retrieve the specific network and ensure it has a basal node for stability
        N = networks.network[j]
        N = add_basal(N)

        # Mechanism 1: Random Extinction
        # Species are removed in a completely stochastic order
        extinction_series = extinction(N; protect = :basal)

        # Log Random Results
        d = Dict{Symbol,Any}(
            :model => networks.model[j],
            :extinction_mechanism => "random",
            :n_rep => networks.n_rep[j],
            :extinction_seq => extinction_series,
        )
        push!(extinction_results, d)

        # Mechanism 2 & 3: Topological and Categorical Extinctions
        # Iterates through both descending and ascending directions
        for descending in [true, false]

            # Numeric (Topological) extinctions based on network metrics
            # Generality (number of resources) or Vulnerability (number of consumers)
            for k in ["generality", "vulnerability"]

                extinction_series = extinction(N, k, descending; protect = :basal)

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

            # Categorical extinctions based on biological traits
            for k in eachindex(hierarchies[1])

                    # Extract species list and the current trait being tested
                    trait_data = df[:, [:species, hierarchies[1][k]]]
                    rename!(trait_data, hierarchies[1][k] => :trait)

                    # Filter to ensure only species with defined traits are targeted
                    filter!(:trait => x -> x ∈ hierarchies[2][k], trait_data)

                    # Determine the sequence of species removal based on hierarchy
                    extinction_list = extinctionsequence(
                        hierarchies[2][k],
                        trait_data;
                        descending = descending,
                    )
                    
                    # Simulate the cascade
                    extinction_series = extinction(N, extinction_list; protect = :basal)

                    # Log Trait-based Results
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

# Save all extinction sequences for subsequent analysis (e.g., Robustness or R50 metrics)
save_object("../data/processed/extinction_seq.jlds", extinction_results)