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
include("lib/lmatrix.jl")
include("lib/niche.jl")
include("lib/random.jl")

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

# Define body size sampling methods
body_methods = ["uniform", "lognormal", "truncated_lognormal"]

# Main simulation loop
for j in 1:n_reps
    y = collect(String, size_classes.size)

    for method in body_methods
        # Generate synthetic body masses
        bodysize = [sample_body_size(sc, method=method) for sc in y]
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
            biomass = bodymass .^ (-3 / 4)
            connectance = rand(Uniform(0.07, 0.15))

            for model ∈ ["adbm", "bodymassratio", "lmatrix", "niche", "pfim_metaweb", "pfim_downsample", "random"]
                N = nothing
                if model == "bodymassratio"
                    N = bmratio(df.species, bodymass)
                elseif model == "pfim_metaweb"
                    N = pfim.PFIM(df, feeding_rules; downsample=false)
                elseif model == "pfim_downsample"
                    N = pfim.PFIM(df, feeding_rules; y=30.0, downsample=true)
                elseif model == "niche"
                    N = nichemodel(df.species, connectance)
                elseif model == "random"
                    links = floor(Int, connectance * (nrow(df)^2))
                    N = randommodel(df.species, links)
                elseif model == "lmatrix"
                    N = lmatrix(df.species, bodymass, is_producer)
                else
                    parameters = adbm_parameters(df, bodymass)
                    N = adbmmodel(df, parameters, biomass)
                end

                if richness(N) > 0
                    push!(networks, (
                        model = model,
                        time = str_cats[1],
                        network = N,
                        n_rep = j,
                        bodysize_method = method
                    ))
                end
            end
        end
    end
end

# Save results
save_object("../data/processed/networks.jlds", networks)