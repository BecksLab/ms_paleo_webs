# General sundry internal functions

using LinearAlgebra
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

    L = links(N)
    S = richness(N)
    l_s = L / S

    D = Dict{Symbol,Any}(
        :richness => richness(N),
        :connectance => SpeciesInteractionNetworks.connectance(N),
        :diameter => _diameter(N),
        :complexity => complexity(N),
        :distance => distancetobase(N, collect(keys(_gen))[ind_maxgen]),
        :generality => std(gen / l_s),
        :vulnerability => std(vul / l_s),
        :redundancy => (L - (S - 1)),
        :S1 => length(findmotif(motifs(Unipartite, 3)[1], remove_cannibals(N)))/(richness(N)^2),
        :S2 => length(findmotif(motifs(Unipartite, 3)[2], remove_cannibals(N)))/(richness(N)^2),
        :S4 => length(findmotif(motifs(Unipartite, 3)[4], remove_cannibals(N)))/(richness(N)^2),
        :S5 => length(findmotif(motifs(Unipartite, 3)[5], remove_cannibals(N)))/(richness(N)^2),
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
    n = zeros(Bool, (species, species))
    for i in axes(n, 1)
        for j in axes(n, 2)
            if N.edges[i, j] == true
                n[i, j] = true
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
        push!(spp, :BASAL_NODE)
        nodes2 = Unipartite(spp)

        # create empty edgelist
        edges2 = Binary(zeros(Bool, (richness(nodes2, 1), richness(nodes2, 2))))

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
        @warn (
            "No BASAL_NODE added as network is already rooted by a single node, original network returned"
        )
        return N
    end
end

"""
TSS(N_real::SpeciesInteractionNetwork{<:Partiteness,<:Binary}, N_sim::SpeciesInteractionNetwork{<:Partiteness,<:Binary})

    Calculates the true test skill statistic of the networks generated 
    via simulated extinctions and those constructed from the known post
    extinction community.

    Note that this assumes that the post extinction community is a 
    subset of the pre extinction community.

    #### References

    Gupta, A., Furrer, R., & Petchey, O. L. (2022). Simultaneously 
    estimating food web connectance and21 structure with uncertainty. 
    Ecology and Evolution, 12(3), e8643. https://doi.org/10.1002/ece3.8643215

"""
function TSS(N_real::SpeciesInteractionNetwork{<:Partiteness,<:Binary}, N_sim::SpeciesInteractionNetwork{<:Partiteness,<:Binary})

    # is post extinction community a subset of the pre extinction community
    species(N_real)

    # get pairwise interactions
    sim = interactions(N_sim)
    real = interactions(N_real)

    # what is in left not in right
    # interactions in simulated not in real
    fp = length(setdiff(sim ,real))
    # interactions in real not in simulated
    fn = length(setdiff(real ,sim))
    
    # what is shared between real and simulated
    tp = length(intersect(sim, real))
    
    # all potential links in real and simulated (all species^2 = num potential links)
    link_pot = length(union(species(N_real), species(N_sim)))^2
    # subtract all other accounted for links to get those missed
    tn = link_pot - (fp + fn + tp)
    
    # eq. to get tss score
    tss = ((tp*tn) - (fp*fn))/((tp + fn)*(fp + tn))

    return tss
end

"""
remove_cannibals(N::SpeciesInteractionNetwork{<:Partiteness,<:Binary})

    Identifies and sets cannibalistic link to zero
"""
function remove_cannibals(N::SpeciesInteractionNetwork{<:Partiteness,<:Binary})

    A = _get_matrix(N)

    A[diagind(A)] .= false

    nodes = Unipartite(species(N))
    edges = Binary(A)
    network = SpeciesInteractionNetwork(nodes, edges)

    return network
end