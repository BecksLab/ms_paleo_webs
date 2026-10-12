using CSV
using DataFrames
using Distributions
using FoodWebTools
using JLD2
using pfim
using SpeciesInteractionNetworks

# helper functions: Imports internal utilities and model-specific logic
include("lib/internals.jl");
include("lib/adbm.jl");
include("lib/bodymass.jl");
include("lib/lmatrix.jl");
include("lib/niche.jl");
include("lib/random.jl");

# Set seed for reproducibility of random body mass and network generation
import Random
Random.seed!(66)

# Data Ingestion: Identify relevant community files in the raw data folder
matrix_names = readdir("../data/raw")
# Filter to keep only datasets containing "Guilds" in the filename
matrix_names = matrix_names[occursin.(r"^.*Guilds.*$", matrix_names)]

# Load external rule sets for feeding interactions and body size distributions
feeding_rules = DataFrame(CSV.File("../data/raw/feeding_rules.csv"))
size_classes = DataFrame(CSV.File("../data/raw/size_classes.csv"))

# Initialize a central DataFrame to store results of all simulated networks
networks = DataFrame(model = String[], time = Any[], network = Any[], n_rep = Any[], targ_Co = Any[]);

# --- Distributions ---

size_bounds = Dict(
    "primary" => (0.01, 0.1),
    "tiny" => (0.1, 10.0),
    "small" => (10.0, 50.0),
    "medium" => (50.0, 100.0),
    "large" => (100.0, 300.0),
    "very_large" => (300.0, 500.0),
    "gigantic" => (500.0, 700.0)
)

σ = 0.5

# Simulation Loop: Run 100 repetitions to account for stochasticity in body sizes and models
n_reps = 100

for j = 1:n_reps

    # Step 1: Generate Synthetic Body Mass Data
    # Assigns a random mass within a specific range based on the 'size' category
    y = collect(String, size_classes.size)
    bodysize = [
                begin
                    lo, hi = size_bounds[s]
                    μ = log((lo + hi)/2)
                    rand(truncated(LogNormal(μ, σ), lo, hi))
                end
                for s in y
            ]

    # Update size_classes with the mass values generated for this specific repetition
    size_classes[!, :bodymass] = bodysize

    # Step 2: Process each community file individually
    for i in eachindex(matrix_names)

        file_name = matrix_names[i]
        # Extract metadata (like time period) from the file name string
        str_cats = split(file_name, r"_")

        # Load and clean the community dataframe
        df = DataFrame(CSV.File.(joinpath("../data/raw/", "$file_name")))
        select!(df, [:Guild, :motility, :tiering, :feeding, :size])
        rename!(df, :Guild => :species)

        # Exclude basal nodes from the species list
        filter!(:species => x -> x != "BASAL_NODE", df)

        # Identify primary producers based on the 'tiering' column
        is_producer = map(==("primary"), string.(df.tiering))

        # Join with size_classes to map synthetic masses to specific species in this community
        bodymass = Vector{Float64}(innerjoin(df, size_classes, on = [:species, :size]).bodymass)

        # Estimate biomass using Metabolic Theory scaling (M^-3/4)
        biomass = bodymass .^ (-3 / 4)

        # Assign a random connectance value for models requiring a fixed link density
        connectance = rand(Uniform(0.01, 0.03))

        # Step 3: Iterate through different ecological network models
        for model ∈ [
            "adbm",
            "bodymassratio",
            "lmatrix",
            "niche",
            "pfim_metaweb",
            "pfim_downsample",
            "random",
        ]

            # Model Branching Logic
            if model == "bodymassratio"
                _N = bmratio(df.species, bodymass)
                adj_mat = Matrix(downsample(_N.edges.edges, :random; target_co=connectance, max_iter=300))
                spp_list = Symbol.(df.species)
                edges = Binary(.!iszero.(adj_mat))
                nodes = Unipartite(spp_list)
                N = SpeciesInteractionNetwork(nodes, edges)
            elseif model == "pfim_downsample"
                # PFIM with random downsampling applied
                _N = PFIM(df, feeding_rules; return_type = :matrix)
                adj_mat = Matrix(downsample(_N, :random; target_co=connectance, max_iter=300))
                spp_list = Symbol.(df.species)
                edges = Binary(.!iszero.(adj_mat))
                nodes = Unipartite(spp_list)
                N = SpeciesInteractionNetwork(nodes, edges)
            elseif model == "niche"
                N = nichemodel(df.species, connectance)
            elseif model == "random"
                # Simple random network based on total possible links and connectance
                N = randommodel(df.species, connectance)
            elseif model == "lmatrix"
                _N = lmatrix(df.species, bodymass, is_producer)
                adj_mat = Matrix(downsample(_N.edges.edges, :random; target_co=connectance, max_iter=300))
                spp_list = Symbol.(df.species)
                edges = Binary(.!iszero.(adj_mat))
                nodes = Unipartite(spp_list)
                N = SpeciesInteractionNetwork(nodes, edges)
            else
                # Default case: Allometric Diet Breadth Model (ADBM)
                model == "adbm"
                parameters = adbm_parameters(df, bodymass)
                _N = adbmmodel(df, parameters, biomass)
                adj_mat = Matrix(downsample(_N.edges.edges, :random; target_co=connectance, max_iter=300))
                spp_list = Symbol.(df.species)
                edges = Binary(.!iszero.(adj_mat))
                nodes = Unipartite(spp_list)
                N = SpeciesInteractionNetwork(nodes, edges)
            end

            # Package result metadata
            d = Dict{Symbol,Any}(
                :model => model,
                :time => str_cats[1],
                :network => N,
                :n_rep => j,
                :targ_Co => connectance
            )

            # Store the network only if it successfully generated species/interactions
            if richness(N) > 0
                push!(networks, d)
            end
        end
    end
end

# Step 4: Final Serialisation
# Save the completed DataFrame containing all repetitions and models to a JLD2 file
save_object("../data/processed/connectance_networks.jlds", networks)

# Pre-allocate a DataFrame to store topological features [cite: 15]
# This defines the schema for the structural analysis of each network
topology = DataFrame(
    model = String[],
    time = Any[],
    n_rep = Any[],
    targ_Co = Any[],
    richness = Int64[],        # Number of nodes (species/guilds)
    connectance = Float64[],   # Realized fraction of possible links
    diameter = Int64[],        # Longest shortest path between any two nodes
    complexity = Float64[],    # Often related to link density or eigenvalue properties
    trophic_level = Float64[], # Mean vertical position in the food web
    distance = Float64[],      # Mean path length between nodes
    generality = Float64[],    # Average number of resources per consumer
    vulnerability = Float64[], # Average number of consumers per resource
    redundancy = Float64[],    # Overlap in ecological roles/links
    NDTI = Float64[],          # normalised difference trophic index
    NDCI = Float64[],          # normalised difference competition index
);

# Analysis Loop: Process every generated network from the previous step 
for i = 1:nrow(networks)

    # Calculate metrics using the internal helper function 
    # This likely calls functions from SpeciesInteractionNetworks.jl internally
    d = _network_summary(networks.network[i])

    # Re-attach metadata from the original 'networks' dataframe 
    d[:model] = networks.model[i]
    d[:n_rep] = networks.n_rep[i]
    d[:time] = networks.time[i]
    d[:targ_Co] = networks.targ_Co[i]

    # Append the resulting dictionary of metrics to the topology table 
    push!(topology, d)
end

# Data Export: Save the final table as a CSV for easy use in plotting or R
CSV.write("../data/processed/topology_connectance.csv", topology)