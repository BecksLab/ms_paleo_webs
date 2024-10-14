using DataFrames
using DataFramesMeta
using SpeciesInteractionNetworks
using StatsBase

"""
_species_removal(N::SpeciesInteractionNetwork, end_richness::Int64)

    Internal function that does the species removal simulations
"""
function _species_removal(
    network_series::Vector{SpeciesInteractionNetwork{<:Partiteness,<:Binary}},
    extinction_list::Vector{Symbol},
    end_richness::Int64,
)

    for (i, sp_to_remove) in enumerate(extinction_list)
        N = network_series[i]
        species_to_keep =
            filter(sp -> sp != sp_to_remove, SpeciesInteractionNetworks.species(N))
        K = subgraph(N, species_to_keep)
        K = simplify(K)
        push!(network_series, K)
        if richness(K) <= end_richness
            break
        else
            continue
        end
    end
    return network_series
end

"""
extinction(N::SpeciesInteractionNetwork, end_richness::Int64)

    Function to simulate random, cascading extinctions of an initial network `N` until 
    the richness is less than or equal to that specified by `end_richness`.
"""
function extinction(
    N::SpeciesInteractionNetwork{<:Partiteness,<:Binary},
    end_richness::Int64,
)
    if richness(N) <= end_richness
        throw(ArgumentError("Richness of staring community is less than final community"))
    end

    extinction_list = StatsBase.shuffle(species(N))
    network_series = Vector{SpeciesInteractionNetwork{<:Partiteness,<:Binary}}()
    # push initial network
    push!(network_series, deepcopy(N))

    return _species_removal(network_series, extinction_list, end_richness)
end

"""
extinction(N::SpeciesInteractionNetwork, extinction_list::Vector{String}, end_richness::Int64)

    Function to simulate cascading extinctions of an initial network `N` until the richness
    is less than or equal to that specified by `end_richness`. The order of species removal
    (extinction) is specified by `extinction_list`
"""
function extinction(
    N::SpeciesInteractionNetwork{<:Partiteness,<:Binary},
    extinction_list::Vector{Symbol},
    end_richness::Int64,
)
    if richness(N) <= end_richness
        throw(ArgumentError("Richness of final community is less than starting community"))
    end
    if !issubset(species(N), extinction_list)
        throw(
            ArgumentError(
                "Species in the network do not match those specified in `extinction_list`",
            ),
        )
    end

    network_series = Vector{SpeciesInteractionNetwork{<:Partiteness,<:Binary}}()
    # push initial network
    push!(network_series, deepcopy(N))

    return _species_removal(network_series, extinction_list, end_richness)
end


"""
extinction_sequence(hierarchy::Vector{Any}, trait_data::DataFrame)

    Determine the order of species extinction for categorical traits. Using a specified hierarchy
"""
function extinction_sequence(hierarchy::Vector{String}, trait_data::DataFrame)
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
function extinction_sequence(trait_dict::Dict{Symbol,Int64}; ascending::Bool = false)
    return collect(keys(sort(trait_dict; byvalue = true, rev = ascending)))
end
