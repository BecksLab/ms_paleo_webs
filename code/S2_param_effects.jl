using CSV
using DataFrames
using Distributions
using JLD2
using pfim
using SpeciesInteractionNetworks

# helper functions
include("lib/internals.jl")
include("lib/adbm.jl")
include("lib/bodymass.jl")
include("lib/bodysize_classes.jl")

# reproducibility
import Random
Random.seed!(66)

# Load data
matrix_names = readdir("../data/raw")
matrix_names = matrix_names[occursin.(r"^.*Guilds.*$", matrix_names)]
feeding_rules = DataFrame(CSV.File("../data/raw/feeding_rules.csv"))
size_classes = DataFrame(CSV.File("../data/raw/size_classes.csv"))

# Initialize central DataFrame
networks = DataFrame(model = String[], time = Any[], network = Any[], n_rep = Any[], bodysize_method = String[])

# Number of repetitions
n_reps = 100

# BMR params sims
# Initialize central DataFrame
networks = DataFrame(time = Any[], network = Any[], n_rep = Any[], param_set = String[])

# Number of repetitions
n_reps = 100

# Define LTM params
param_sets = Dict(
    "baseline" => (α=1.41, β=3.73, γ=-1.90),
    "narrow_niche" => (α=1.41, β=2.5, γ=-1.90),
    "wide_niche" => (α=1.41, β=5.0, γ=-1.90),
    "shifted_optimum_low" => (α=1.0, β=3.73, γ=-1.90),
    "shifted_optimum_high" => (α=1.8, β=3.73, γ=-1.90),
    "strict_threshold" => (α=1.41, β=3.73, γ=-2.5),
    "relaxed_threshold" => (α=1.41, β=3.73, γ=-1.0)
)

# Main simulation loop
for j in 1:n_reps
    
    # Step 1: Generate Synthetic Body Mass Data
    # Assigns a random mass within a specific range based on the 'size' category
    y = collect(String, size_classes.size)
    bodysize = (
        y ->
            y == "tiny" ? rand(Uniform(0.1, 10.0)) :
            y == "small" ? rand(Uniform(10.0, 50.0)) :
            y == "medium" ? rand(Uniform(50.0, 100.0)) :
            y == "large" ? rand(Uniform(100.0, 300.0)) :
            y == "very_large" ? rand(Uniform(300.0, 500.0)) :
            y == "gigantic" ? rand(Uniform(500.0, 700.0)) : y
    ).(y)

    # Update size_classes with the mass values generated for this specific repetition
    size_classes[!, :bodymass] = bodysize

    for file_name in matrix_names
            str_cats = split(file_name, r"_")

            # Load and clean community data
            df = DataFrame(CSV.File(joinpath("../data/raw/", file_name)))
            select!(df, [:Guild, :motility, :tiering, :feeding, :size])
            rename!(df, :Guild => :species)
            filter!(:species => x -> x != "BASAL_NODE", df)

            is_producer = map(==("primary"), string.(df.tiering))
            bodymass = Vector{Float64}(innerjoin(df, size_classes, on=[:species, :size]).bodymass)

            for (pname, pvals) in param_sets
                N = bmratio(df.species, bodymass;
                α=pvals.α, β=pvals.β, γ=pvals.γ)

                if richness(N) > 0
                    push!(networks, (
                        time = str_cats[1],
                        network = N,
                        n_rep = j,
                        param_set = pname
                    ))
                end
            end
        end
end

# Pre-allocate a DataFrame to store topological features [cite: 15]
# This defines the schema for the structural analysis of each network
topology = DataFrame(
    time = Any[],
    n_rep = Any[],
    param_set = String[],
    richness = Int64[],      # Number of nodes (species/guilds)
    connectance = Float64[],   # Realized fraction of possible links
    diameter = Int64[],        # Longest shortest path between any two nodes
    complexity = Float64[],    # Often related to link density or eigenvalue properties
    trophic_level = Float64[], # Mean vertical position in the food web
    distance = Float64[],      # Mean path length between nodes
    generality = Float64[],    # Average number of resources per consumer
    vulnerability = Float64[], # Average number of consumers per resource
    redundancy = Float64[],    # Overlap in ecological roles/links
    S1 = Float64[],            # Structural motif 1
    S2 = Float64[],            # Structural motif 2
    S4 = Float64[],            # Structural motif 4 
    S5 = Float64[],            # Structural motif 5 
);

# Analysis Loop: Process every generated network from the previous step 
for i = 1:nrow(networks)

    # Calculate metrics using the internal helper function 
    # This likely calls functions from SpeciesInteractionNetworks.jl internally
    d = _network_summary(networks.network[i])

    # Re-attach metadata from the original 'networks' dataframe 
    d[:n_rep] = networks.n_rep[i]
    d[:time] = networks.time[i]
    d[:param_set] = networks.param_set[i]

    # Append the resulting dictionary of metrics to the topology table 
    push!(topology, d)
end

# Data Export: Save the final table as a CSV for easy use in plotting or R
CSV.write("../data/processed/topology_bmr_params.csv", topology)