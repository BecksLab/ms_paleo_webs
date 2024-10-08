using DataFramesMeta
using SpeciesInteractionNetworks
using StatsBase

"""
extinction(N::SpeciesInteractionNetwork{<:Partiteness,<:Binary}, end_richness::Int64)

    Function to simulate secondary extinctions using an initial network N unitl the richness is 
    less than or equal to that specified by `end_richness`.
"""

function extinction(N::SpeciesInteractionNetwork{<:Partiteness,<:Binary}, end_richness::Int64)
    if richness(N) <= end_richness
        throw(ArgumentError("Richness of final community is less than starting community"))
    end

    network_series = Vector{SpeciesInteractionNetwork{<:Partiteness,<:Binary}}(undef, richness(N)+1)
    network_series[1] = deepcopy(N)
    final_network = deepcopy(N)
    # order of species to remove
    extinction_sequence = StatsBase.shuffle(SpeciesInteractionNetworks.species(N))

    for (i, sp_to_remove) in enumerate(extinction_sequence)
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
_extinction_sequence(N::SpeciesInteractionNetwork{<:Partiteness,<:Binary}, end_richness::Int64)

    Determine the order of species extinction based on the desired mechanism.
"""

function _extinction_sequence(
    N::SpeciesInteractionNetwork{<:Partiteness,<:Binary},
    mechanism::String;
    traits::DataFrame = DataFrame(A = [1, 2])
    )

    # data checks
    if mechanism ∉ ["random", "size_descend", "size_ascend", "tiering_descend", "tiering_ascend", "motility_fast_non", "motility_non_fast", "generality_ascend", "generality_descend", "vulnerability_ascend", "vulnerability_descend"]
        error("$(mechanim) is not a recognised extinction mechanim")
    end
    if mechanism ∈ ["size_descend", "size_ascend", "tiering_descend", "tiering_ascend", "motility_fast_non", "motility_non_fast"]
        for (i, v) in enumerate(["species", "motility", "tiering", "feeding", "size"])
            if v ∉ names(traits)
                error("Missing $(v) variable as a column in DataFrame, add or rename")
            end
        end
        if SpeciesInteractionNetworks.species(N) ∉ names(traits)
            error("Not all species in network N are listed in the traits data")
        end
    end

    if mechanism == "random"
        spp_list = StatsBase.shuffle(SpeciesInteractionNetworks.species(N))
    elseif mechanism == "size_descend"
        order = ["very_large", "large", "medium", "small", "tiny"]
        @rorderby traits findfirst(==(:size), order)
        spp_list = Symbol.(df.species)
    elseif mechanism == "size_ascend"
        order = ["tiny", "small", "medium", "large", "very_large"]
        @rorderby traits findfirst(==(:size), order)
        spp_list = Symbol.(df.species)
    elseif mechanism == "tiering_descend"
    elseif mechanism == "tiering_ascend"
    elseif mechanism == "motility_fast_non"
    elseif mechanism == "motility_non_fast"
    elseif mechanism == "generality_ascend"
        gen = SpeciesInteractionNetworks.generality(N)
        spp_list = keys(sort(gen; byvalue = true))
    elseif mechanism == "generality_descend"
        spp_list = gen = SpeciesInteractionNetworks.generality(N)
        keys(sort(gen; byvalue = true, rev=true))
    elseif mechanism == "vulnerability_ascend"
        vul = vulnerability(N)
        spp_list = keys(sort(vul; byvalue = true))
    else mechanism == "vulnerability_descend"
        vul = vulnerability(N)
        spp_list = keys(sort(vul; byvalue = true, rev=true))
    end
    return spp_list
end