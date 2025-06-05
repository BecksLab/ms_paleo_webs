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

# size classes (for creating continuous body sizes)
size_classes = DataFrame(CSV.File("data/raw/size_classes.csv"))

# df to store networks
networks = DataFrame(model = String[], time = Any[], network = Any[], n_rep = Any[]);

# number of network reps
n_reps = 100

for j = 1:n_reps

        # create some synthetic bodysize data (based on some distributions)
        # spp body sizes will be constant across communities but vary by rep
        y = collect(String, size_classes.size)

        bodysize =
            (
                y ->
                    y == "tiny" ? rand(Uniform(0.1, 10.0)) :
                    y == "small" ? rand(Uniform(10.0, 50.0)) :
                    y == "medium" ? rand(Uniform(50.0, 100.0)) :
                    y == "large" ? rand(Uniform(100.0, 300.0)) :
                    y == "very_large" ? rand(Uniform(300.0, 500.0)) :
                    y == "gigantic" ? rand(Uniform(500.0, 700.0)) : y
            ).(y)

        # add to size classes df
        size_classes[!,:bodymass] = bodysize
    
        for i in eachindex(matrix_names)
    
            file_name = matrix_names[i]
            # get relevant info from slug
            str_cats = split(file_name, r"_")
    
            # import data frame
            df = DataFrame(CSV.File.(joinpath("data/raw/", "$file_name")))
            select!(df, [:Guild, :motility, :tiering, :feeding, :size])
            rename!(df, :Guild => :species)
    
            # remove BASAL_NODE for now...
            filter!(:species => x -> x != "BASAL_NODE", df)
    
            # specify if producer (basal node)
            is_producer = map(==("primary"), string.(df.tiering))

            # get the bodysizes of species only present in df
            bodymass = Vector{Float64}(innerjoin(df, size_classes, on = [:species, :size]).bodymass)

            # create some mock abundance/biomass values using a *very* basic scaling law
            biomass = bodymass .^ (-3 / 4)

        # specify connectance for niche/random model
        # TODO could possibly have this be 'dynamic' based on the Co of other networks...
        connectance = rand(Uniform(0.07, 0.15))

        for model ∈ ["adbm", "bodymassratio", "lmatrix", "niche", "pfim_metaweb", "pfim_downsample", "random"]

            if model == "bodymassratio"
                N = bmratio(df.species, bodymass)
                N = randomdraws(N) # from probabilistic to binary
            elseif model == "pfim_metaweb"
                N = pfim.PFIM(df, feeding_rules; downsample = false)
            elseif model == "pfim_downsample"
                N = pfim.PFIM(df, feeding_rules; downsample = true)
            elseif model == "niche"
                N = structuralmodel(NicheModel, nrow(df), connectance)
            elseif model == "random"
                links = floor(Int, connectance * (nrow(df)^2))
                N = randommodel(df.species, links)
            elseif model == "lmatrix"
                N = lmatrix(df.species, bodymass, is_producer)
            else
                model == "adbm"
                parameters = adbm_parameters(df, bodymass)
                N = adbmmodel(df, parameters, biomass)
            end

            d = Dict{Symbol,Any}(
                :model => model,
                :time => str_cats[1],
                :network => N,
                :n_rep => j,
            )

            # only push if network exists
            if richness(N) > 0
                push!(networks, d)
            end
        end
    end
end

# write networks as object
save_object("data/processed/networks.jlds", networks)
