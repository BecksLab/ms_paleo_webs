using CSV
using DataFrames
using JLD2
using ProgressMeter
using SpeciesInteractionNetworks

# set seed
import Random
Random.seed!(66)

# import model libs
include("../../code/lib/bodymass/bodymass.jl")
include("../../code/lib/adbm/adbm.jl")
include("../../code/lib/internals.jl")

# import datasets
#nz_summaries = DataFrame(CSV.File("../notebooks/data/processed/nz_summary.csv"))
nz_networks = load_object("../data/raw/nz_networks.jlds")
comm_data = DataFrame(CSV.File("../data/raw/taxa_dry_weight_abundance.csv"))
rename!(comm_data, :dw => :bodymass)
rename!(comm_data, Symbol("no.m2") => :abundance)
rename!(comm_data, :taxa  => :species)

# summarise NZ networks

nz_summaries = topo_df();

for i in eachindex(nz_networks)

   d = _network_summary(nz_networks[i].network)

   D = Dict{Symbol,Any}()
   D[:id] = nz_networks[i].id
   D[:model] = "NA"
   D[:richness] = d[:richness]
   D[:connectance] = d[:connectance]
   D[:complexity] = d[:complexity]
   D[:deficiency] = d[:deficiency]
   D[:distance] = d[:distance]
   D[:basal] = d[:basal]
   D[:top] = d[:top]
   D[:generality] = d[:generality]
   D[:vulnerability] = d[:vulnerability]
   D[:S1] = d[:S1]
   D[:S2] = d[:S2]
   D[:S4] = d[:S4]
   D[:S5] = d[:S5]
   D[:network] = nz_networks[i].network

    push!(nz_summaries, D)
end

# write summaries as .csv
CSV.write("../data/processed/nz_summary.csv",
            nz_summaries[:, setdiff(names(nz_summaries), ["model", "network"])])

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
Co = 0.1;

# number of networks to generate
n_reps = 1000;

@showprogress for _ = 1:n_reps
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
end


# write summaries as .csv
CSV.write("../data/processed/topology_models.csv",
    topology[:, setdiff(names(topology), ["network"])])