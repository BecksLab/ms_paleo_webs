using CSV
using DataFrames
using JLD2

include("lib/extinctions.jl")
include("lib/internals.jl")

# set seed
import Random
Random.seed!(66)

# create trait dataframe
# get the name of all communities
matrix_names = readdir("../data/clean/trait")
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
            joinpath("../data/clean/trait", "$file_name.csv"),
        ),
    )
    append!(traits, df)
end

traits = unique(traits)

# Extinction sequence

# for results
extinction_results = topo_df()
# add addition column
extinction_results[!,:time] = String[]
extinction_results[!,:extinction_mechanism] = String[]
extinction_results[!,:links] = Int64[]
# remove unused columns
select!(extinction_results, Not(:network))

matrix_names = readdir("../data/processed/networks")

for i in eachindex(matrix_names)

    # import predicted network for specific model
    file_name = matrix_names[i]
    df = load_object("../data/processed/networks/$file_name")

    # now we loop throuhg the different time periods
    for j in ["pre", "during"]
            # select only the pre extinction community
            pre_comm = df[occursin.(j, df.id), :]
            # get post extinction richness
            post_rich = richness(df[occursin.("post", df.id), :network][1])

            # random extinction
            extinction_series = extinction(pre_comm.network[1], post_rich)
    
            # summarise extinction network
            D = _network_summary(extinction_series[end])
    
            # append additional info
            D[:model] = pre_comm.model[1]
            D[:extinction_mechanism] = "random"
            D[:id] = pre_comm.id[1]
            D[:time] = j

            # send to results
            push!(extinction_results, D)

            # numeric extinctions
            for k in ["generality", "vulnerability"]

                f = getfield(Main, Symbol(k))

                extinction_list = extinction_sequence(f(pre_comm.network[1]))

                # generate extinction sequence
                extinction_series = extinction(pre_comm.network[1], extinction_list, post_rich)
    
                # summarise extinction network
                D = _network_summary(extinction_series[end])
    
                # append additional info
                D[:model] = pre_comm.model[1]
                D[:extinction_mechanism] = k
                D[:id] = pre_comm.id[1]
                D[:time] = j
    
                push!(extinction_results, D)
            end
    end
end

# write summaries as .csv
CSV.write("../data/processed/extinctions/extinctions.csv", extinction_results[:,setdiff(names(extinction_results), ["network"])])