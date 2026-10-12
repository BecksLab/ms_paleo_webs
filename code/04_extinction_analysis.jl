using CSV
using DataFrames
using Extinctions
using JLD2
using ProgressMeter
using SpeciesInteractionNetworks

include("lib/internals.jl")

# Set seed for reproducibility in stochastic robustness simulations
import Random
Random.seed!(66)

# Data Ingestion: Load simulated sequences and full network objects
networks = load_object("../data/processed/networks.jlds")

# Filter for G2 (post-extinction) to serve as the "ground truth" for comparison
filter!(:time => x -> x == "G2", networks)

# Load raw community lists to identify which species survived and which are new arrivals
post_df = CSV.read("../data/raw/G2_Guilds.csv", DataFrame)
pre_df = CSV.read("../data/raw/G1_Guilds.csv", DataFrame)

# Identify species overlap between the two time periods
shared_spp = intersect(post_df.Guild, pre_df.Guild)
spp_originated = setdiff(post_df.Guild, pre_df.Guild) # New arrivals (ignored for simulation comparison)

# The "target richness" used to find the point in the simulation that matches the real world
post_rich = length(shared_spp)

# Comparative Loop: Matches simulations to real data snapshots
# again we do this per a model
models = replace.(readdir("../data/processed/extinction_seq"), ".jld2" => "")

@showprogress "Processing models" for model_name in models

    extinctions = load_object("../data/processed/extinction_seq/$(model_name).jld2")

    # Matching G2 networks for this model only
    model_networks = subset(networks, :model => ByRow(==(model_name)))

    topology = DataFrame(
        model = String[],
        extinction_mechanism = String[],
        n_rep = Int[],
        richness = Int[],
        connectance = Float64[],
        diameter = Int[],
        complexity = Float64[],
        trophic_level = Float64[],
        distance = Float64[],
        generality = Float64[],
        vulnerability = Float64[],
        redundancy = Float64[],
        NDTI = Float64[],
        NDCI = Float64[],
        resilience = Float64[],
        rep = Int[],
    )

    tss = DataFrame(
        model = String[],
        extinction_mechanism = String[],
        n_rep = Int[],
        tss_link = Float64[],
        tss_node = Float64[],
    )

    β_div = DataFrame(
        model = String[],
        extinction_mechanism = String[],
        n_rep = Int[],
        β_div = Float64[],
        β_type = String[],
    )

    for i in 1:nrow(extinctions)

        _ext_seq = extinctions.extinction_seq[i]

        if length(_ext_seq) <= 1
            continue
        end

        net_ind = findfirst(x -> x < post_rich, richness.(_ext_seq)) - 1

        if !(net_ind isa Int) || richness(_ext_seq[net_ind]) == 0
            continue
        end

        N_real = subset(
            model_networks,
            [:n_rep] => ByRow(==(extinctions.n_rep[i])),
        )

        real_net = add_basal(N_real.network[1])

        real_spp = species(real_net)
        filter!(x -> x ∉ Symbol.(spp_originated), real_spp)
        real_N = subgraph(real_net, real_spp)

        sim_N = _ext_seq[net_ind]

        # ---------- TSS ----------
        tss_link, tss_node = TSS(real_N, sim_N, _ext_seq[1])

        push!(tss, (
            model = model_name,
            extinction_mechanism = extinctions.extinction_mechanism[i],
            n_rep = extinctions.n_rep[i],
            tss_link = tss_link,
            tss_node = tss_node,
        ))

        # ---------- Topology ----------
        d = _network_summary(sim_N)

        push!(topology, (
            model = model_name,
            extinction_mechanism = extinctions.extinction_mechanism[i],
            n_rep = extinctions.n_rep[i],
            richness = d[:richness],
            connectance = d[:connectance],
            diameter = d[:diameter],
            complexity = d[:complexity],
            trophic_level = d[:trophic_level],
            distance = d[:distance],
            generality = d[:generality],
            vulnerability = d[:vulnerability],
            redundancy = d[:redundancy],
            NDTI = d[:NDTI],
            NDCI = d[:NDCI],
            resilience = resilience(_ext_seq),
            rep = i,
        ))

        # ---------- Beta diversity ----------
        for β in (βS, βWN, βOS)

            β_vals = SpeciesInteractionNetworks.betadiversity(β, real_N, sim_N)

            a = β_vals.shared
            b = β_vals.right
            c = β_vals.left

            _β = (a + b + c) / ((2a + b + c) / 2) - 1

            push!(β_div, (
                model = model_name,
                extinction_mechanism = extinctions.extinction_mechanism[i],
                n_rep = extinctions.n_rep[i],
                β_div = _β,
                β_type = string(β),
            ))
        end
    end

    CSV.write("../data/processed/extinction_topology/$(model_name).csv", topology)
    CSV.write("../data/processed/extinction_tss/$(model_name).csv", tss)
    CSV.write("../data/processed/extinction_betadiv/$(model_name).csv", β_div)

    GC.gc()
end

# 4. Robustness Curve Generation
# Focuses on the pfim_metaweb model to evaluate network resistance to random cascading loss
networks = load_object("../data/processed/networks.jlds")
filter!(:time => x -> x == "G2", networks)

models = unique(networks.model)

spread = collect(1:99)
n_reps = 500

@showprogress "Calculating robustness" for model_name in models

    model_networks = subset(networks, :model => ByRow(==(model_name)))

    robustness_vals = DataFrame(
        model = String[],
        time = String[],
        n_rep = Int[],
        threshold = Int[],
        robustness = Float64[],
    )

    for _ in 1:n_reps
        for i in 1:nrow(model_networks)

            N = remove_cannibals(model_networks.network[i])
            N = add_basal(N)

            for threshold in spread

                rob = robustness(
                    N;
                    threshold = threshold,
                    mechanism = :cascade,
                )

                push!(robustness_vals, (
                    model = model_name,
                    time = model_networks.time[i],
                    n_rep = model_networks.n_rep[i],
                    threshold = threshold,
                    robustness = rob,
                ))
            end
        end
    end

    CSV.write("../data/processed/robustness/$(model_name).csv", robustness_vals)

    GC.gc()
end