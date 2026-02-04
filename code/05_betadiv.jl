using Combinatorics
using CSV
using DataFrames
using JLD2
using SpeciesInteractionNetworks

# Load internal utilities
include("lib/internals.jl")

# Initialize storage for beta diversity results
# This will track variation over time and between different model pairings
β_div = DataFrame(
    time = Any[],
    left = String[],   # ID of the first network in the comparison
    right = Any[],     # ID of the second network in the comparison
    β_div = Float64[], # Calculated beta diversity value
    β_type = String[], # Type of beta diversity (e.g., βS, βWN, βOS)
)

# Data Ingestion: Load the full suite of generated networks
networks = load_object("../data/processed/networks.jlds")

# Filtering: Remove "niche" and "random" models to focus the comparison 
# on biologically-constrained models (ADBM, PFIM, etc.)
filter!(:model => x -> x ∉ ["niche", "random"], networks)

# Analysis Groups: Identify unique time periods (e.g., G1, G2)
time_groups = unique(networks.time)

for i in eachindex(time_groups)

    # Isolate all networks belonging to the current time period
    df = filter(:time => x -> x == time_groups[i], networks)
    
    # Create a unique identifier combining model name and repetition number 
    # to facilitate unique pairwise combinations
    df.uniq_id = map(join, zip(df.model, df.n_rep))

    # Combinatorial Loop: Generate every possible pair of unique networks within this time group
    for (x, y) in combinations(unique(df[:, :uniq_id]), 2)

        # Retrieve the network objects for both IDs in the pair
        U = df[occursin.(x, df.uniq_id), :network][1]
        V = df[occursin.(y, df.uniq_id), :network][1]
        
        # Iterate through standard interaction beta diversity metrics:
        # βS:  Overall dissimilarity
        # βWN: Dissimilarity due to re-wiring (interactions between shared species)
        # βOS: Dissimilarity due to species turnover
        for β ∈ [βS, βWN, βOS]
                
                # Calculate raw components using SpeciesInteractionNetworks.jl
                β_vals = SpeciesInteractionNetworks.betadiversity(β, U, V)

                a = β_vals.shared # Interactions present in both networks
                b = β_vals.right  # Interactions unique to network V
                c = β_vals.left   # Interactions unique to network U

                # Calculate the Sorensen-based dissimilarity index
                # (a+b+c)/((2a+b+c)/2) - 1 is a variation of common dissimilarity formulas
                _β = (a + b + c)/((2a + b + c)/2) - 1

                # Package metadata and results
                res_dict = Dict{Symbol,Any}()
                res_dict[:time] = time_groups[i]
                res_dict[:left] = x
                res_dict[:right] = y
                res_dict[:β_div] = _β
                res_dict[:β_type] = "$β"

                # Append to the final results table
                push!(β_div, res_dict)
            end
    end
end

# Data Export: Save for final statistical comparison of model agreement
CSV.write("../data/processed/beta_div.csv", β_div)