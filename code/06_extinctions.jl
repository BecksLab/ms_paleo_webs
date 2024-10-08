using CSV
using DataFrames

include("../lib/extinctions.jl")

# set seed
import Random
Random.seed!(66)

# create trait dataframe
# get the name of all communities
matrix_names = readdir(joinpath("data", "clean", "trait"))
matrix_names = replace.(matrix_names, ".csv" => "")

# get the traits data
traits = DataFrame(
    species = Any[],
    feeding = Any[],
    motility = Any[],
    tiering = Any[],
    size = Any[],
)

for i in eachindex(matrix_names)    
    file_name = matrix_names[i]
    df = DataFrame(
        CSV.File.(
            joinpath("data", "clean", "trait", "$file_name.csv"),
        ),
    )
    append!(traits, df)
end

traits = unique(traits)

# Extinction sequence
matrix_names = readdir(joinpath("data", "processed", "networks"))

for i in eachindex(matrix_names)

    # Import predicted network 
    file_name = matrix_names[i]
    df = load_object("data/processed/networks/$file_name")

    # select only the pre extinction community
    pre_comm = df[occursin.("pre", df.id), :]
    post_rich = richness(df[occursin.("post", df.id), :network][1])

    # generate extinction sequence

    N = extinction(pre_comm.network[1], post_rich)

    return _network_summary(N)
end