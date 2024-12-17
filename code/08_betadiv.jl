using Combinatorics
using CSV
using DataFrames
using JLD2
using SpeciesInteractionNetworks

include("lib/internals.jl")

β_div = DataFrame(
    model = String[],
    left = String[],
    right = String[],
    β_int = Int64[],
    β_spp = Int64[],
)

matrix_names = readdir("../data/processed/networks")

for i in eachindex(matrix_names)

    # import predicted network for specific model
    file_name = matrix_names[i]
    df = load_object("../data/processed/networks/$file_name")

    # for each combination of datasets
    for (x, y) in combinations(["pre", "during", "post"], 2)
        
        U = df[occursin.(x, df.id), :network][1]
        V = df[occursin.(y, df.id), :network][1]
        β_int = SpeciesInteractionNetworks.betadiversity(βS,U,V)
        β_spp = SpeciesInteractionNetworks.betadiversity(βWN,U,V)

        # cllate
        b = Dict{Symbol,Any}()
        b[:model] = df.model[1]
        b[:left] = x
        b[:right] = y
        b[:β_int] = β_int.shared
        b[:β_spp] = β_spp.shared

        # send to results
        push!(β_div, b)
    end  
end