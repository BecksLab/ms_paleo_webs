using AlgebraOfGraphics
using CairoMakie
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
extinctions = load_object("../data/processed/extinction_seq.jlds")
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

# Initialize storage for Topology, TSS (Skill), and Beta Diversity (Similarity)
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
    tss_link = Float64[],
    tss_node = Float64[],
)

β_div = DataFrame(
    model = String[],
    extinction_mechanism = Any[],
    n_rep = Any[],
    β_div = Float64[],
    β_type = String[],
)

# Comparative Loop: Matches simulations to real data snapshots
@showprogress "Getting topology" for i = 1:nrow(extinctions)

    _ext_seq = extinctions.extinction_seq[i]

    if length(_ext_seq) > 1
        # Locate the specific step in the extinction sequence where richness matches the G2 target
        net_ind = findfirst(x -> x < post_rich, richness.(_ext_seq)) - 1

        # Check for valid index and non-empty network
        if typeof(net_ind) == Int64 && richness(_ext_seq[net_ind]) > 0

            # 1. TSS Calculation (Predictive Skill)
            # Find the corresponding real-world network object
            N_real = filter(
                [:model, :n_rep] =>
                    (x, y) -> x == extinctions.model[i] && y == extinctions.n_rep[i],
                networks,
            )

            real_net = add_basal(N_real.network[1])
            
            # Remove "originated" species from the real network because the simulation
            # (which starts at G1) cannot predict new arrivals.
            real_spp = species(real_net)
            filter!(x -> x ∉ Symbol.(spp_originated), real_spp)
            real_N = subgraph(real_net, real_spp)
            
            sim_N = _ext_seq[net_ind]

            # Compare simulation results to reality
            tss_link, tss_node = TSS(real_N, sim_N, _ext_seq[1])

            t = Dict{Symbol,Any}(
                :model => extinctions.model[i],
                :extinction_mechanism => extinctions.extinction_mechanism[i],
                :n_rep => extinctions.n_rep[i],
                :tss_link => tss_link,
                :tss_node => tss_node,
            )
            push!(tss, t)

            # 2. Network Topology and Resilience Analysis
            d = _network_summary(_ext_seq[net_ind])
            d[:model] = extinctions.model[i]
            d[:n_rep] = extinctions.n_rep[i]
            d[:extinction_mechanism] = extinctions.extinction_mechanism[i]
            d[:resilience] = resilience(_ext_seq) # Measure how well the system resists collapse
            d[:rep] = i
            push!(topology, d)

            # 3. Beta Diversity (Structure Turnover)
            # Calculates shared/unique links and nodes between real and simulated networks
            for β ∈ [βS, βWN, βOS]
                β_vals = SpeciesInteractionNetworks.betadiversity(β, real_N, sim_N)
                a = β_vals.shared
                b = β_vals.right
                c = β_vals.left

                # Calculate Sorensen-style turnover
                _β = (a + b + c)/((2a + b + c)/2) - 1

                b_dict = Dict{Symbol,Any}()
                b_dict[:model] = extinctions.model[i]
                b_dict[:extinction_mechanism] = extinctions.extinction_mechanism[i]
                b_dict[:n_rep] = i
                b_dict[:β_div] = _β
                b_dict[:β_type] = "$β"
                push!(β_div, b_dict)
            end
        end
    end
end

# Export analysis results
CSV.write("../data/processed/extinction_topology.csv", topology)
CSV.write("../data/processed/extinction_tss.csv", tss)
CSV.write("../data/processed/extinction_betadiv.csv", β_div)

# 4. Robustness Curve Generation
# Focuses on the pfim_metaweb model to evaluate network resistance to random cascading loss
networks = load_object("../data/processed/networks.jlds")
filter!(:model => x -> x == "pfim_metaweb", networks)

# Test survival at varying thresholds of node loss (1% to 99%)
spread = collect(1:1:99)
n_reps = 500

robustness_vals = DataFrame(
    time = Any[],
    threshold = Any[],
    robustness = Any[],
);

for _ = 1:n_reps
    @showprogress "Calculating robustness" for i in 1:4
        # Standardize network by removing self-loops (cannibalism) and adding energy source
        N = remove_cannibals(networks.network[i])
        N = add_basal(N)

        for j in eachindex(spread)
            # Simulate a full secondary extinction cascade at the current threshold
            rob = robustness(N; threshold = spread[j], mechanism = :cascade)
    
            D = DataFrame(
                time = networks.time[i],
                threshold = spread[j],
                robustness = rob,
            )
            append!(robustness_vals, D)
        end
    end    
end

# Visualization: Generate a smoothed robustness curve plot
layer = data(robustness_vals) * mapping(:threshold, :robustness, color = :time => nonnumeric) * (smooth())
fig = draw(layer)
save("../figures/robustness_curve.png", fig.figure)

# Final data export
CSV.write("../data/processed/robustness.csv", robustness_vals)