using Combinatorics
using CSV
using DataFrames
using JLD2
using SpeciesInteractionNetworks

include("lib/internals.jl")

β_div = DataFrame(
    time = Any[],
    left = String[],
    right = Any[],
    β_div = Float64[],
    β_type = String[],
)

# same comment RE elegance
networks = load_object("../data/processed/networks.jlds")
# remove random and niche networks...
filter!(:model => x -> x ∉ ["niche", "random"], networks)

# get the different communities
time_groups = unique(networks.time)

for i in eachindex(time_groups)

    df = filter(:time => x -> x == time_groups[i], networks)
    # make unique ID for making combos
    df.uniq_id = map(join, zip(df.model, df.n_rep))

    # for each combination of datasets
    for (x, y) in combinations(unique(df[:, :uniq_id]), 2)

        U = df[occursin.(x, df.uniq_id), :network][1]
        V = df[occursin.(y, df.uniq_id), :network][1]
        
        for β ∈ [βS, βWN, βOS]
                
                β_vals = SpeciesInteractionNetworks.betadiversity(β, U, V)

                a = β_vals.shared
                b = β_vals.right
                c = β_vals.left

                _β = (a + b + c)/((2a + b + c)/2) - 1

                b = Dict{Symbol,Any}()
                b[:time] = time_groups[i]
                b[:left] = x
                b[:right] = y
                b[:β_div] = _β
                b[:β_type] = "$β"

                # send to results
                push!(β_div, b)

            end
    end
end

# write outputs as .csv
CSV.write("../data/processed/beta_div.csv", β_div)
