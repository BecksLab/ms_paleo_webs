using DataFrames
using DataFramesMeta
using SpeciesInteractionNetworks
using StatsBase

"""
extinction(N::SpeciesInteractionNetwork, end_richness::Int64)

    Function to simulate random, cascading extinctions of an initial network `N` until 
    the richness is less than or equal to that specified by `end_richness`.
"""
function extinction(N::SpeciesInteractionNetwork{<:Partiteness,<:Binary}, 
                    end_richness::Int64)
    if richness(N) <= end_richness
        throw(ArgumentError("Richness of final community is less than starting community"))
    end

    network_series = Vector{SpeciesInteractionNetwork{<:Partiteness,<:Binary}}(undef, richness(N)+1)
    network_series[1] = deepcopy(N)
    final_network = deepcopy(N)
    extinction_list = StatsBase.shuffle(SpeciesInteractionNetworks.species(N))

    for (i, sp_to_remove) in enumerate(extinction_list)
            species_to_keep = filter(sp -> sp != sp_to_remove, SpeciesInteractionNetworks.species(network_series[i]))
            K = subgraph(N, species_to_keep)
            K = simplify(K)
            network_series[i+1] = K
        if richness(K) <= end_richness
            final_network = K
            break
        else
            continue
        end
    end
    return final_network
end

"""
extinction(N::SpeciesInteractionNetwork, extinction_list::Vector{String}, end_richness::Int64)

    Function to simulate cascading extinctions of an initial network `N` until the richness
    is less than or equal to that specified by `end_richness`. The order of species removal
    (extinction) is specified by `extinction_list`
"""
function extinction(N::SpeciesInteractionNetwork{<:Partiteness,<:Binary}, 
                    extinction_list::Vector{String}, 
                    end_richness::Int64)
    if richness(N) <= end_richness
        throw(ArgumentError("Richness of final community is less than starting community"))
    end
    if !issubset(species(N), extinction_list)
        throw(ArgumentError("Species in the network do not match those specified in `extinction_list`"))
    end

    network_series = Vector{SpeciesInteractionNetwork{<:Partiteness,<:Binary}}(undef, richness(N)+1)
    network_series[1] = deepcopy(N)
    final_network = deepcopy(N)

    for (i, sp_to_remove) in enumerate(extinction_list)
            species_to_keep = filter(sp -> sp != sp_to_remove, SpeciesInteractionNetworks.species(network_series[i]))
            K = subgraph(N, species_to_keep)
            K = simplify(K)
            network_series[i+1] = K
        if richness(K) <= end_richness
            final_network = K
            break
        else
            continue
        end
    end
    return final_network
end


"""
extinction_sequence(hierarchy::Vector{Any}, trait_data::DataFrame)

    Determine the order of species extinction for categorical traits. Using a specified hierarchy
"""
function extinction_sequence(
    hierarchy::Vector{String},
    trait_data::DataFrame
    )
    # data checks
    if !issubset(String.(trait_data.trait), hierarchy)
        error("Not all traits in `traits_data` are listed in `hierarchy`")
    end

    df = @rorderby trait_data findfirst(==(:trait), hierarchy)
    return String.(df.species)
end

"""
extinction_sequence(hierarchy::Vector{Any}, trait_data::DataFrame)

    Determine the order of species extinction for numeric traits.
"""
function extinction_sequence(
    trait_dict::Dict{String, Int64};
    ascending::Bool=false
    )
    return keys(sort(trait_dict; byvalue = true, rev=ascending))
end