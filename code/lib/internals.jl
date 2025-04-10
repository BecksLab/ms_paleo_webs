# General sundry internal functions

using SpeciesInteractionNetworks
using Statistics

"""
_network_summary(N::SpeciesInteractionNetwork{<:Partiteness, <:Binary})

    returns the 'summary statistics' for a network
"""
function _network_summary(N::SpeciesInteractionNetwork{<:Partiteness,<:Binary})

    A = _get_matrix(N)

    _gen = SpeciesInteractionNetworks.generality(N)
    gen = collect(values(_gen))
    vul = collect(values(SpeciesInteractionNetworks.vulnerability(N)))
    ind_maxgen = findmax(gen)[2]

    D = Dict{Symbol,Any}(
        :richness => richness(N),
        :links => links(N),
        :connectance => SpeciesInteractionNetworks.connectance(N),
        :diameter => _diameter(N),
        :complexity => complexity(N),
        :distance => distancetobase(N, collect(keys(_gen))[ind_maxgen]),
        :basal => sum(vec(sum(A, dims = 2) .== 0)),
        :top => sum(vec(sum(A, dims = 1) .== 0)),
        :generality => std(gen),
        :vulnerability => std(vul),
        :redundancy => (links(N) - (richness(N) - 1)),
        :S1 => length(findmotif(motifs(Unipartite, 3)[1], N)),
        :S2 => length(findmotif(motifs(Unipartite, 3)[2], N)),
        :S4 => length(findmotif(motifs(Unipartite, 3)[4], N)),
        :S5 => length(findmotif(motifs(Unipartite, 3)[5], N)),
    )

    return D
end

"""
    maxrank(N::SpeciesInteractionNetwork{<:Partiteness, <:Binary}}

Returns the maximum possible rank of a Network
"""
function maxrank(N::SpeciesInteractionNetwork{<:Partiteness,<:Binary})
    return minimum(size(N))
end

"""
_get_matrix(N::SpeciesInteractionNetwork{<:Partiteness, <:Binary})

    Internal function to return a matrix of interactions from a
    SpeciesInteractionNetwork
"""
function _get_matrix(N::SpeciesInteractionNetwork{<:Partiteness,<:Binary})

    species = richness(N)
    n = zeros(Int64, (species, species))
    for i in axes(n, 1)
        for j in axes(n, 2)
            if N.edges[i, j] == true
                n[i, j] = 1
            end
        end
    end

    return n
end

"""
_diameter(N::SpeciesInteractionNetwork{<:Partiteness,<:Binary})

    Calculates the diameter of a food web. Where diameter is the longest 
    shortest path between two nodes
"""
function _diameter(N::SpeciesInteractionNetwork{<:Partiteness,<:Binary})

    # extract species names
    spp = SpeciesInteractionNetworks.species(N)
    # empty vector for storing shortest path for each spp
    shortpath = zeros(Int64, length(spp))

    # get shortest path
    for i in eachindex(spp)
        shortpath[i] = length(shortestpath(N, spp[i]))
    end

    # return max shortest path
    return findmax(shortpath)[1]
end

_parser(x) = parse(Int, x)

"""
add_basal(N::SpeciesInteractionNetwork{<:Partiteness,<:Binary})

    Add a single basal node to networks for which there are multiple
    producer (generality of zero) species
"""
function add_basal(N::SpeciesInteractionNetwork{<:Partiteness,<:Binary})

    # find producer spp (gen == 0)
    gen = SpeciesInteractionNetworks.generality(N)
    filter!(v -> last(v) == 0, gen)

    # only add basal node if the number of species with gen 0 > 1
    if length(gen) > 1

        # get current interactions
        intxns = interactions(N)

        # create producer-BASAL_SPP interactions tuple
        producer_spp = collect(keys(gen))
        prod_basal = tuple.(producer_spp, :BASAL_NODE, true)

        # combine current and producer-BASAL_SPP interactions
        interactions_all = vcat(intxns, prod_basal)

        # build new network architecture

        # add BASAL_NODE to spp list
        spp = SpeciesInteractionNetworks.species(N)
        push!(spp,:BASAL_NODE)
        nodes2 = Unipartite(spp)

        # create empty edgelist
        edges2 = Binary(zeros(Bool, (richness(nodes2,1), richness(nodes2, 2))))

        # new empty network
        network2 = SpeciesInteractionNetwork(nodes2, edges2)
    
        # add interactions using interactions_all
        for i in eachindex(interactions_all)
            network2[(interactions_all[i][1], interactions_all[i][2])...] = true
        end
        return network2
    elseif isempty(gen)
        @warn ("No producer species, original network returned")
        return N
    else
        @warn ("No BASAL_NODE added as network is already rooted by a single node, original network returned")
        return N
    end
end