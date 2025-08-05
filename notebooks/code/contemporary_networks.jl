using CSV
using DataFrames
using Distributions
using JLD2
using pfim
using ProgressMeter
using SpeciesInteractionNetworks

# set seed
import Random
Random.seed!(66)

# import model libs
include("../../code/lib/bodymass.jl")
include("../../code/lib/adbm.jl")
include("../../code/lib/lmatrix.jl")
include("../../code/lib/random.jl")
include("../../code/lib/internals.jl")

# import datasets
#nz_summaries = DataFrame(CSV.File("../notebooks/data/processed/nz_summary.csv"))
nz_networks = load_object("../data/raw/nz_networks.jlds")
comm_data = DataFrame(CSV.File("../data/raw/taxa_dry_weight_abundance.csv"))
rename!(comm_data, :dw => :bodymass)
rename!(comm_data, Symbol("no.m2") => :abundance)
rename!(comm_data, :taxa => :species)

# summarise NZ networks

nz_summaries = DataFrame(
    network = Any[],
    model = Any[],
    id = Any[],
    richness = Int64[],
    connectance = Float64[],
    diameter = Int64[],
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
);

for i in eachindex(nz_networks)

    d = _network_summary(nz_networks[i].network)

    
    d[:model] = "NA"
    d[:network] = nz_networks[i].network
    d[:id] = nz_networks[i].id

    push!(nz_summaries, d)
end

# write summaries as .csv
CSV.write(
    "../data/processed/nz_summary.csv",
    nz_summaries[:, setdiff(names(nz_summaries), ["model", "network"])],
)

# side quest: find 'producer' species
producer = DataFrame(species = String15[], tiering = String15[], id = String15[])

for i in eachindex(nz_networks)
    gen = generality(nz_networks[i].network)

    tiering = collect(values(gen))

    p = DataFrame(
        species = collect(keys(gen)),
        tiering = ifelse.(tiering .== 0, "producer", "non"),
        id = nz_networks[i].id,
    )

    producer = vcat(p, producer)
end

# master df
topology = DataFrame(
    network = Any[],
    model = String[],
    id = Any[],
    richness = Int64[],
    connectance = Float64[],
    diameter = Int64[],
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
);

# connectance (for Niche model)
Co = 0.1;

# number of networks to generate
n_reps = 1000;

@showprogress for j = 1:n_reps
    for i = 1:nrow(nz_summaries)

        site = nz_summaries.id[i]
        spp_rich = nz_summaries.richness[i]

        prods = filter(:id => x -> x == site, producer)
        comm = filter(Symbol("food.web") => x -> x == site, comm_data)

        df = innerjoin(comm, prods, on = :species)

        is_producer = map(==("producer"), string.(df.tiering))

        bodymass = abs.(rand.(Uniform.(df.bodymass .- 0.00001, df.bodymass .+ 0.00001)))
        biomass = df.abundance

        for model ∈ [
            "adbm",
            "bodymassratio",
            "lmatrix",
            "niche",
            "random",
        ]

            if model == "bodymassratio"
                N = bmratio(df.species, bodymass)
            elseif model == "niche"
                N = structuralmodel(NicheModel, nrow(df), Co)
            elseif model == "random"
                links = floor(Int, Co * (nrow(df)^2))
                N = randommodel(df.species, links)
            elseif model == "lmatrix"
                N = lmatrix(df.species, bodymass, is_producer)
            else
                model == "adbm"
                parameters = adbm_parameters(df, bodymass)
                N = adbmmodel(df, parameters, biomass)
            end

            d = _network_summary(N)

            d[:model] = model
            d[:network] = N
            d[:id] = site

            # only push if network exists
            if richness(N) > 0
                push!(topology, d)
            end
        end

    end
end


# write summaries as .csv
CSV.write(
    "../data/processed/topology_models.csv",
    topology[:, setdiff(names(topology), ["network"])],
)
