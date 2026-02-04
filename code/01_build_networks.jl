using CSV
using DataFrames
using Distributions
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
networks = DataFrame(model = String[], time = Any[], network = Any[], n_rep = Any[]);

# Simulation Loop: Run 100 repetitions to account for stochasticity in body sizes and models
n_reps = 100

for j = 1:n_reps

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
        connectance = rand(Uniform(0.07, 0.15))

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
                N = bmratio(df.species, bodymass)
            elseif model == "pfim_metaweb"
                # Probabilistic Feeding Interaction Model (Metaweb version)
                N = pfim.PFIM(df, feeding_rules; downsample = false)
            elseif model == "pfim_downsample"
                # PFIM with downsampling applied
                N = pfim.PFIM(df, feeding_rules; y = 30.0, downsample = true)
            elseif model == "niche"
                N = nichemodel(df.species, connectance)
            elseif model == "random"
                # Simple random network based on total possible links and connectance
                links = floor(Int, connectance * (nrow(df)^2))
                N = randommodel(df.species, links)
            elseif model == "lmatrix"
                N = lmatrix(df.species, bodymass, is_producer)
            else
                # Default case: Allometric Diet Breadth Model (ADBM)
                model == "adbm"
                parameters = adbm_parameters(df, bodymass)
                N = adbmmodel(df, parameters, biomass)
            end

            # Package result metadata
            d = Dict{Symbol,Any}(
                :model => model,
                :time => str_cats[1],
                :network => N,
                :n_rep => j,
            )

            # Store the network only if it successfully generated species/interactions
            if richness(N) > 0
                push!(networks, d)
            end
        end
    end
end

# Step 4: Final Serialization
# Save the completed DataFrame containing all repetitions and models to a JLD2 file
save_object("../data/processed/networks.jlds", networks)