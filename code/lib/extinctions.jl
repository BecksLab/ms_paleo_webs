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
        # end if target richness reached
        if richness(K) == end_richness
            push!(network_series, K)
            break
        # if richenss below target then we break without pushing
        elseif richness(K) < end_richness
            break
        # continue removing species
        else
            push!(network_series, K)
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
extinction_sequence(hierarchy::Vector{Any}, trait_data::DataFrame; descending::Bool = false)

    Determine the order of species extinction for categorical traits. Using a specified hierarchy
"""
function extinction_sequence(hierarchy::Vector{String}, trait_data::DataFrame; descending::Bool = false)
    # data checks
    if !issubset(String.(trait_data.trait), hierarchy)
        error("Not all traits in `traits_data` are listed in `hierarchy`")
    end

    # reverse order of trait hierarchy if descending is true
    if descending == true
        hierarchy = reverse(hierarchy)
    end

    df = @rorderby trait_data findfirst(==(:trait), hierarchy)
    return Symbol.(df.species)
end

"""
extinction_sequence(hierarchy::Vector{Any}, trait_data::DataFrame; descending::Bool = false)

    Determine the order of species extinction for numeric traits.
"""
function extinction_sequence(trait_dict::Dict{Symbol,Int64}; descending::Bool = false)
    return collect(keys(sort(trait_dict; byvalue = true, rev = descending)))
end
