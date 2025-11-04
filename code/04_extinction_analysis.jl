using AlgebraOfGraphics
using CairoMakie
using CSV
using DataFrames
using Extinctions
using JLD2
using ProgressMeter
using SpeciesInteractionNetworks

include("lib/internals.jl")

# set seed
import Random
Random.seed!(66)

# import networks object
extinctions = load_object("../data/processed/extinction_seq.jlds")

# also import the 'known' networks for tss calcs
networks = load_object("../data/processed/networks.jlds")
# we only care about the post extinction community
filter!(:time => x -> x == "G2", networks)

# find richness of post extinction community (but also only species shared between G1 and G2)
post_df = CSV.read("../data/raw/G2_Guilds.csv", DataFrame)
pre_df = CSV.read("../data/raw/G1_Guilds.csv", DataFrame)

# this is the list of species shared
shared_spp = intersect(post_df.Guild, pre_df.Guild)
# this is the list of species that entered the community post extinction
spp_originated = setdiff(post_df.Guild, pre_df.Guild)

# this is the 'target richness' of the extinction sequence
post_rich = length(shared_spp)

topology = DataFrame(
    model = String[],
    extinction_mechanism = Any[],
    n_rep = Any[],
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
    resilience = Float64[],
    rep = Int64[],
);

tss = DataFrame(
    model = String[],
    extinction_mechanism = Any[],
    n_rep = Any[],
    tss = Float64[],
)

β_div = DataFrame(
    model = String[],
    extinction_mechanism = Any[],
    n_rep = Any[],
    β_int_shared = Int64[],
    β_int_all = Int64[],
    β_spp = Int64[],
    links_left = Int64[],
    links_right = Int64[],
)

@showprogress "Getting topology" for i = 1:nrow(extinctions)

    _ext_seq = extinctions.extinction_seq[i]

    if length(_ext_seq) > 1

        net_ind = findfirst(x -> x < post_rich, richness.(_ext_seq)) - 1

        # this is to catch networks that 
        if typeof(net_ind) == Int64 && richness(_ext_seq[net_ind]) > 0

            # tss score

            # get the real network
            N_real = filter(
                [:model, :n_rep] =>
                    (x, y) -> x == extinctions.model[i] && y == extinctions.n_rep[i],
                networks,
            )

            real_net = add_basal(N_real.network[1])
            # subset so we remove species that originated
            # get spp in real and filter that list against the originated spp list
            real_spp = species(real_net)
            filter!(x -> x ∉ Symbol.(spp_originated), real_spp)

            real_N = subgraph(real_net, real_spp)
            sim_N = _ext_seq[net_ind]

            _tss = TSS(real_N, sim_N, _ext_seq[1])

            t = Dict{Symbol,Any}(
                :model => extinctions.model[i],
                :extinction_mechanism => extinctions.extinction_mechanism[i],
                :n_rep => extinctions.n_rep[i],
                :tss => _tss,
            )

            push!(tss, t)

            # network topology
            d = _network_summary(_ext_seq[net_ind])

            d[:model] = extinctions.model[i]
            d[:n_rep] = extinctions.n_rep[i]
            d[:extinction_mechanism] = extinctions.extinction_mechanism[i]
            d[:resilience] = resilience(_ext_seq)
            d[:rep] = i

            push!(topology, d)

            # beta div - as an alternative to tss
            β_spp = SpeciesInteractionNetworks.betadiversity(βS, real_N, sim_N)
            β_int_all = SpeciesInteractionNetworks.betadiversity(βWN, real_N, sim_N)
            β_int_shared = SpeciesInteractionNetworks.betadiversity(βOS, real_N, sim_N)
            # collate
            b = Dict{Symbol,Any}()
            b[:model] = extinctions.model[i]
            b[:extinction_mechanism] = extinctions.extinction_mechanism[i]
            b[:n_rep] = extinctions.n_rep[i]
            b[:β_int_all] = β_int_all.shared
            b[:β_int_shared] = β_int_shared.shared
            b[:β_spp] = β_spp.shared
            b[:links_left] = links(real_N)
            b[:links_right] = links(sim_N)

            # send to results
            push!(β_div, b)

        end
    end
end

# write summaries as .csv
CSV.write("../data/processed/extinction_topology.csv", topology)
CSV.write("../data/processed/extinction_tss.csv", tss)
CSV.write("../data/processed/extinction_betadiv.csv", β_div)


# robustness curves

# using only random extinctions
networks = load_object("../data/processed/networks.jlds")
# we only care about the metaweb pfim for now
filter!(:model => x -> x == "pfim_metaweb", networks)

# here we want to generate robustness curves from the spread 1:99
spread = collect(1:1:99)

n_reps = 500

robustness_vals = DataFrame(
    time = Any[],
    threshold = Any[],
    robustness = Any[],
);

for _ = 1:n_reps

    @showprogress "Calculating robustness" for i in 1:4

        N = remove_cannibals(networks.network[i])
        N = add_basal(N)

        for j in eachindex(spread)
    
            rob = robustness(N;
                        threshold = spread[j],
                        mechanism = :cascade)
    
            D = DataFrame(
                time = networks.time[i],
                threshold = spread[j],
                robustness = rob,
                )
    
            # send to results
            append!(robustness_vals, D)
        end
    end    
end

# robustness curves

layer = data(robustness_vals) * mapping(:threshold, :robustness, color = :time => nonnumeric) * (smooth())
fig = draw(layer)

figure = fig.figure
save("../figures/robustness_curve.png", figure)

# write summaries as .csv
CSV.write(
    "../data/processed/robustness.csv",
    robustness_vals,
)
