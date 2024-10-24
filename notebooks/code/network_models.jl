using CSV
using DataFrames
using JLD2
using SpeciesInteractionNetworks

# set seed
import Random
Random.seed!(66)

# import model libs
include("../../code/lib/bodymass/bodymass.jl")
include("../../code/lib/adbm/adbm.jl")
include("../../code/lib/internals.jl")

# import datasets
nz_summaries = DataFrame(CSV.File("../notebooks/data/processed/nz_summary.csv"))
nz_networks = load_object("../notebooks/data/raw/nz_networks.jlds")
comm_data = DataFrame(CSV.File("../notebooks/data/raw/taxa_dry_weight_abundance.csv"))
rename!(comm_data, :dw => :bodymass)
rename!(comm_data, Symbol("no.m2") => :abundance)
rename!(comm_data, :taxa  => :species)

# side quest: find 'producer' species

producer = 
    DataFrame(
            species = String15[],
            tiering = String15[],
            id = String15[]
            )

for i in eachindex(nz_networks)
    gen = generality(nz_networks[i].network)

    tiering = collect(values(gen))

    p = DataFrame(
            species = collect(keys(gen)),
            tiering = ifelse.(tiering .== 0, "producer", "non"),
            id = nz_networks[i].id
        )

    producer = vcat(p, producer)
end



# master df
topology = topo_df();

# connectance (for Niche model)
Co = 0.1

for i in 1:nrow(nz_summaries)

    site = nz_summaries.id[i]
    spp_rich = nz_summaries.richness[i]

    prods = filter(:id => x -> x == site, producer)
    comm = filter(Symbol("food.web") => x -> x == site, comm_data)

    df = innerjoin(comm, prods, on = :species)

    # adbm
    d = model_summary(df, site, "adbm"; bodymass = df.bodymass, biomass = df.abundance)
    push!(topology, d)
    # niche
    d = model_summary(df, site, "niche"; connectance = Co)
    push!(topology, d)
    # bodymass
    d = model_summary(df, site, "bodymassratio"; bodymass = df.bodymass)
    push!(topology, d)

end

# write summaries as .csv
CSV.write("../notebooks/data/processed/topology_models.csv",
    topology[:, setdiff(names(topology), ["network"])])