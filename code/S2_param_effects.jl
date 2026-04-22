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
networks = DataFrame(
    model = String[],
    time = Any[],
    network = Any[],
    n_rep = Int[],
    param_set = String[]
)

# Number of repetitions
n_reps = 100

# Define params sets (BMR and ADBM)
bmr_param_sets = Dict(
    "baseline" => (α=1.41, β=3.73, γ=-1.90),
    "narrow_niche" => (α=1.41, β=2.5, γ=-1.90),
    "wide_niche" => (α=1.41, β=5.0, γ=-1.90),
    "shifted_optimum_low" => (α=1.0, β=3.73, γ=-1.90),
    "shifted_optimum_high" => (α=1.8, β=3.73, γ=-1.90),
    "strict_threshold" => (α=1.41, β=3.73, γ=-2.5),
    "relaxed_threshold" => (α=1.41, β=3.73, γ=-1.0)
)

adbm_baseline = (
    e=1.0, a_adbm=0.0189, ai=-0.491, aj=-0.465,
    b=0.401, h_adbm=1.0, hi=1.0, hj=1.0,
    n=1.0, ni=-0.75, Hmethod=:ratio, Nmethod=:original
)

adbm_param_sets = Dict(
    "baseline" => (),
    "high_attack" => (a_adbm=0.05,),
    "low_attack" => (a_adbm=0.005,),
    "strict_ratio" => (b=0.3,),
    "relaxed_ratio" => (b=0.6,),
    "high_handling" => (h_adbm=2.0,),
    "low_handling" => (h_adbm=0.5,),
    "power_handling" => (Hmethod=:power,),
    "biomass_driven" => (Nmethod=:biomass,)
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

        df = DataFrame(CSV.File(joinpath("../data/raw/", file_name)))
        select!(df, [:Guild, :motility, :tiering, :feeding, :size])
        rename!(df, :Guild => :species)
        filter!(:species => x -> x != "BASAL_NODE", df)

        bodymass = Vector{Float64}(innerjoin(df, size_classes, on=[:species, :size]).bodymass)

        # -------------------------
        # BMR MODEL
        # -------------------------
        for (pname, pvals) in bmr_param_sets

            N = bmratio(df.species, bodymass;
                α=pvals.α, β=pvals.β, γ=pvals.γ)

            if richness(N) > 0
                push!(networks, (
                    model = "BMR",
                    time = str_cats[1],
                    network = N,
                    n_rep = j,
                    param_set = pname
                ))
            end
        end

        # -------------------------
        # ADBM MODEL
        # -------------------------
        for (pname, pvals) in adbm_param_sets

            params = merge(adbm_baseline, pvals)

            adbm_params = adbm_parameters(
                df,
                bodymass;
                params...
            )

            biomass = bodymass  # or experiment later

            N = adbmmodel(df, adbm_params, biomass)

            if richness(N) > 0
                push!(networks, (
                    model = "ADBM",
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
    model = String[],
    time = Any[],
    n_rep = Int[],
    param_set = String[],
    richness = Int[],
    connectance = Float64[],
    diameter = Int[],
    complexity = Float64[],
    trophic_level = Float64[],
    distance = Float64[],
    generality = Float64[],
    vulnerability = Float64[],
    redundancy = Float64[],
    S1 = Float64[],
    S2 = Float64[],
    S4 = Float64[],
    S5 = Float64[],
)

# Analysis Loop: Process every generated network from the previous step 
for i = 1:nrow(networks)

    d = _network_summary(networks.network[i])

    d[:model] = networks.model[i]
    d[:n_rep] = networks.n_rep[i]
    d[:time] = networks.time[i]
    d[:param_set] = networks.param_set[i]

    push!(topology, d)
end

# Data Export: Save the final table as a CSV for easy use in plotting or R
CSV.write("../data/processed/topology_param_effects.csv", topology)