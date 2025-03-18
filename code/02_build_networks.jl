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
matrix_names = matrix_names[occursin.(r"^.*species.*$", matrix_names)]

# feeding rules
feeding_rules = DataFrame(CSV.File("data/raw/feeding_rules.csv"))

# specify connectance for niche model
# TODO could possibly have this be 'dynamic' based on the Co of other networks...
connectance = 0.1

# df to store networks
networks = DataFrame(
    model = String[], 
    location = String[],
    time = Any[],
    network = Any[],
);

for i in eachindex(matrix_names)

    file_name = matrix_names[i]
    # get relevant info from slug
    str_cats = split(file_name, r"_")

    # import data frame
    df = DataFrame(CSV.File.(joinpath("data/raw/", "$file_name")))

    for time = 1:5

        # select correct time period
        time_df = filter(df -> occursin.("$time", df.time_pre_during_post), df)
        
        # add primary node
        push!(time_df, ["primary" "primary" "primary" "primary" "primary" "$time"])

        # specify if producer (basal node)
        is_producer = map(==("primary"), string.(time_df.tiering))

        # create some synthetic bodysize data (based on some distributions)
        y = collect(String, time_df.size)

        bodymass = (y -> y == "tiny" ? rand(Uniform(0.1, 0.3)) :
                    y == "small" ? rand(Uniform(0.2, 0.5)) :
                    y == "medium" ? rand(Uniform(0.4, 0.7)) :
                    y == "large" ? rand(Uniform(0.6, 0.9)) :
                    y == "primary" ? rand(Uniform(0.1, 0.3)) :
                    y).(y)

        # create some mock abundance/biomass values using a *very* basic scaling law
        biomass = bodymass .^ (-3 / 4)

        # only build network for community (richness > 1)
        if nrow(time_df) > 1

            for model ∈ ["adbm", "bodymassratio", "lmatrix", "niche", "pfim", "random"]
                
                if model == "bodymassratio"
                    N = bmratio(time_df.species, bodymass)
                    N = randomdraws(N) # from probabilistic to binary
                elseif model == "pfim"
                    N = pfim.PFIM(time_df, feeding_rules)
                elseif model == "niche"
                    N = structuralmodel(NicheModel, nrow(time_df), connectance)
                elseif model == "random"
                    links = floor(Int, connectance * (nrow(time_df)^2))
                    N = randommodel(time_df.species, links)
                elseif model == "lmatrix"
                    N = lmatrix(time_df.species, bodymass, is_producer)
                    N = randomdraws(N) # from probabilistic to binary
                else model == "adbm"
                    parameters = adbm_parameters(time_df, bodymass)
                    N = adbmmodel(time_df, parameters, biomass)
                end

                d = Dict{Symbol,Any}(
                    :model => model,
                    :time => time,
                    :location => str_cats[1],
                    :network => N,
                    )
            
                push!(networks, d)

            end
        end
    end
end

# write networks as object
save_object(
    "data/processed/networks.jlds",
    networks,
)