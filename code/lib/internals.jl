# General sundry internal functions

using LinearAlgebra
using SpeciesInteractionNetworks
using Statistics

# import other scripts with functions
include("diameter.jl")

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

    tls = trophic_level(N)

    D = Dict{Symbol,Any}(
        :richness => richness(N),
        :connectance => SpeciesInteractionNetworks.connectance(N),
        :diameter => diameter(A),
        :complexity => complexity(N),
        :trophic_level => findmax(collect(values(tls)))[2],
        :distance => distancetobase(N, collect(keys(_gen))[ind_maxgen]),
        :generality => std(gen / l_s),
        :vulnerability => std(vul / l_s),
        :redundancy => (L - (S - 1))/S,
        :S1 =>
            length(
                findmotif(motifs(Unipartite, 3)[1], remove_cannibals(N)),
            )/(richness(N)^2),
        :S2 =>
            length(
                findmotif(motifs(Unipartite, 3)[2], remove_cannibals(N)),
            )/(richness(N)^2),
        :S4 =>
            length(
                findmotif(motifs(Unipartite, 3)[4], remove_cannibals(N)),
            )/(richness(N)^2),
        :S5 =>
            length(
                findmotif(motifs(Unipartite, 3)[5], remove_cannibals(N)),
            )/(richness(N)^2),
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
        intxns = SpeciesInteractionNetworks.interactions(N)

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
function TSS(
    N_real::SpeciesInteractionNetwork{<:Partiteness,<:Binary},
    N_sim::SpeciesInteractionNetwork{<:Partiteness,<:Binary},
    N_original::SpeciesInteractionNetwork{<:Partiteness,<:Binary};
)

    # get pairwise interactions
    sim = SpeciesInteractionNetworks.species(N_sim)
    real = SpeciesInteractionNetworks.species(N_real)
    original = SpeciesInteractionNetworks.species(N_original)

    # NODE tss

    # what is in left not in right
    # spp that were in original and not in sim && then not in real
    # spp that should go extinct
    extinct = setdiff(original, real)
    # make sure those spp are not in sim
    tn = length(setdiff(extinct, sim))
    # spp in simulated not in real (species that should've been removed)
    fp = length(setdiff(sim, real))
    # species in real not in simulated (ie the 'wrong' species removed)
    fn = length(setdiff(real, sim))
    # what is shared between real and simulated
    tp = length(intersect(sim, real))

    # eq. to get tss score
    tss_node = (tp/(tp+fn)) + (tn/(tn+fp)) - 1

    # LINK tss

    # first we need to subset the networks to only include spp that are shared (to not punish for spp turnover)
    shared = intersect(sim, real)
    # both pre and post extinction
    sim_mat = _get_matrix(subgraph(N_sim, shared))
    real_mat = _get_matrix(subgraph(N_real, shared))

    # what is in left not in right
    # spp that were in original and not in sim && then not in real
    # spp that should go extinct
    tp = sum((sim_mat .== 1) .& (real_mat .== 1))
    tn = sum((sim_mat .== 0) .& (real_mat .== 0))
    fp = sum((sim_mat .== 1) .& (real_mat .== 0))
    fn = sum((sim_mat .== 0) .& (real_mat .== 1))

    # eq. to get tss score
    tss_link = (tp/(tp+fn)) + (tn/(tn+fp)) - 1


    return (tss_link, tss_node)
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

"""
trophic_level(N::SpeciesInteractionNetwork)

    Calculates the trophic level of all species in a network using the average 
    shortest path from the prey of species 𝑖 to a basal species

    Williams, Richard J., and Neo D. Martinez. 2004. “Limits to Trophic Levels 
    and Omnivory in Complex Food Webs: Theory and Data.” The American Naturalist 
    163 (3): 458–68. https://doi.org/10.1086/381964.
"""
function trophic_level(N::SpeciesInteractionNetwork)

    A = _get_matrix(N) # Ensure A is dense for inversion.
    S = size(A, 1) # Species richness.
    out_degree = sum(A; dims = 2)
    D = -(A ./ out_degree) # Diet matrix.
    D[isnan.(D)] .= 0.0
    D[diagind(D)] .= 1.0 .- D[diagind(D)]
    # Solve with the inverse matrix.
    inverse = iszero(det(D)) ? pinv : inv
    tls = inverse(D) * ones(S)

    # create dictionary
    Dict(zip(species(N),tls))

end
