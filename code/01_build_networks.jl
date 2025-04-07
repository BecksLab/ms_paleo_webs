using CSV
using DataFrames
using Distributions
using JLD2
using pfim
using SpeciesInteractionNetworks

# helper functions
include("lib/internals.jl");

# additional functions for building networks
include("lib/adbm.jl");
include("lib/bodymass.jl");
include("lib/lmatrix.jl");
include("lib/random.jl");

# set seed
import Random
Random.seed!(66)

# get the name of all communities
matrix_names = readdir("data/raw")
# select only species datasets
matrix_names = matrix_names[occursin.(r"^.*Guilds.*$", matrix_names)]

# feeding rules
feeding_rules = DataFrame(CSV.File("data/raw/feeding_rules.csv"))

# specify connectance for niche model
# TODO could possibly have this be 'dynamic' based on the Co of other networks...
connectance = 0.1

# df to store networks
networks = DataFrame(
    model = String[], 
    time = Any[],
    network = Any[],
);

for i in eachindex(matrix_names)

    file_name = matrix_names[i]
    # get relevant info from slug
    str_cats = split(file_name, r"_")

    # import data frame
    df = DataFrame(CSV.File.(joinpath("data/raw/", "$file_name")))

    select!(df, [:Guild, :motility_detail, :tiering, :feeding, :size])

    rename!(df, :Guild => :species)
    rename!(df, :motility_detail => :motility)

    # specify if producer (basal node)
    is_producer = map(==("primary"), string.(df.tiering))
    
    # create some synthetic bodysize data (based on some distributions)
    y = collect(String, df.size)

    bodymass = (y -> y == "tiny" ? rand(Uniform(0.1, 0.3)) :
                y == "small" ? rand(Uniform(0.2, 0.5)) :
                y == "medium" ? rand(Uniform(0.4, 0.7)) :
                y == "large" ? rand(Uniform(0.6, 0.9)) :
                y == "very_large" ? rand(Uniform(0.8, 1.1)) :
                y == "gigantic" ? rand(Uniform(1.0, 1.3)) :
                y == "primary" ? rand(Uniform(0.1, 0.3)) :
                y).(y)

    # create some mock abundance/biomass values using a *very* basic scaling law
    biomass = bodymass .^ (-3 / 4)
    
    for model ∈ ["adbm", "bodymassratio", "lmatrix", "niche", "pfim", "random"]
                
        if model == "bodymassratio"
            N = bmratio(df.species, bodymass)
            N = randomdraws(N) # from probabilistic to binary
        elseif model == "pfim"
            N = pfim.PFIM(df, feeding_rules)
        elseif model == "niche"
            N = structuralmodel(NicheModel, nrow(df), connectance)
        elseif model == "random"
            links = floor(Int, connectance * (nrow(df)^2))
            N = randommodel(df.species, links)
        elseif model == "lmatrix"
            N = lmatrix(df.species, bodymass, is_producer)
            N = randomdraws(N) # from probabilistic to binary
        else model == "adbm"
            parameters = adbm_parameters(df, bodymass)
            N = adbmmodel(df, parameters, biomass)
        end

        d = Dict{Symbol,Any}(
            :model => model,
            :time => str_cats[1],
            :network => N,)
            
        # only push if network exists
        if richness(N) > 0
            push!(networks, d)
        end
    end
end

# write networks as object
save_object(
    "data/processed/networks.jlds",
    networks,
)