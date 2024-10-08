using DataFrames
using SpeciesInteractionNetworks

include(joinpath("rules.jl"))

# create struct to contain traits for a species
mutable struct PFIMspecies
    species::Symbol
    feeding_trait::feeding
    size_trait::sizes
    motility_trait::motility
    tiering_trait::tier
end

"""
_pfim_community(data::DataFrame)

    Internal function that takes a DataFrame and transforms it into an
    array containing `PFIMspecies``. This is done so that the 'metadata' 
    is infomred by a species as opposed to a matrix that needs to be
    indexed...

"""
function _pfim_community(data::DataFrame)


    PFIMcommunity = Array{PFIMspecies,1}(undef, nrow(data))

    for i in eachindex(PFIMcommunity)
        _PFIMsp = PFIMspecies(
            Symbol(data.species[i]),
            getfield(Main, Symbol(data.feeding[i]))(),
            getfield(Main, Symbol(data.size[i]))(),
            getfield(Main, Symbol(data.motility[i]))(),
            getfield(Main, Symbol(data.tiering[i]))(),
        )
        PFIMcommunity[i] = _PFIMsp
    end
    return PFIMcommunity
end

"""
    _pfim_link(consumer::PFIMspecies, resource::PFIMspecies)

    Internal function that actually determines the link feasibility for
    a given pair of PFIMspecies. Retuns either 1 (link present) or 0
    (link absent)

"""
function _pfim_link(consumer::PFIMspecies, resource::PFIMspecies)

    # sum outcomes from the four rules
    trait_sum =
        feeding_rules(consumer.feeding_trait, resource.feeding_trait) +
        motility_rules(consumer.motility_trait, resource.motility_trait) +
        tiering_rules(consumer.tiering_trait, resource.tiering_trait) +
        size_rules(consumer.size_trait, resource.size_trait)

    # if all conditions are met (sums to 4) then link is present
    if trait_sum == 4
        link = 1
    else
        link = 0
    end
    return link
end

"""
    _PFIM_network(PFIMcommunity::Vector{PFIMspecies})

    Internal function that constructs a network for a given 
    PFIMcommunity by returning the oucome from _pfim_link().

"""
function _PFIM_network(PFIMcommunity::Vector{PFIMspecies})

    S = length(PFIMcommunity)
    int_matrix = zeros(Bool, (S, S))

    # populate matrix
    for i in eachindex(PFIMcommunity)
        for j in eachindex(PFIMcommunity)
            int_matrix[i, j] = _pfim_link(PFIMcommunity[i], PFIMcommunity[j])
        end
    end

    # data check
    if sum(int_matrix) == 0
            error("No viable interactions for this community")
        end

    # create SpeciesInteractionNetwork
    nodes = Unipartite(getproperty.(PFIMcommunity, :species))
    edges = Binary(int_matrix)
    network = SpeciesInteractionNetwork(nodes, edges)
    return network, int_matrix
end


"""
   PFIM(data::DataFrame; y::Float64 = 2.5)

    Takes a data frame and impliments the feeding rules to determine the
    feasibility of links between species. As well as applying the link
    distribution downsampling approach.
    
    #### References
    
    Shaw, Jack O., Alexander M. Dunhill, Andrew P. Beckerman, Jennifer A.
    Dunne, and Pincelli M. Hull. 2024. “A Framework for Reconstructing 
    Ancient Food Webs Using Functional Trait Data.” 
    https://doi.org/10.1101/2024.01.30.578036.
"""
function PFIM(data::DataFrame; y::Float64 = 2.5, downsample::Bool = true)

    # data checks
    for (i, v) in enumerate(["species", "motility", "tiering", "feeding", "size"])
        if v ∉ names(data)
            error("Missing $(v) variable as a column in DataFrame, add or rename")
        end
    end

    S = nrow(data)
    PFIMcommunity = _pfim_community(data)
    network, matrix = _PFIM_network(PFIMcommunity)
    link_dist = zeros(Float64, S)

    # downsampling protocol
    if downsample == true
        # get link distribution
        for i in eachindex(data.species)
            sp = data.species[i]
            r = generality(network, Symbol(sp))
            E = exp(log(S) * (y - 1) / y)
            link_dist[i] = exp(r / E)
        end

        # create probabilistic int matrix
        prob_matrix = zeros(AbstractFloat, (S, S))
        for i in axes(matrix, 1)
            for j in axes(matrix, 2)
                if matrix[i, j] == true
                    prob_matrix[i, j] = link_dist[i]
                end
            end
        end
        # make probabanilistic
        prob_matrix = prob_matrix ./ maximum(prob_matrix)

        edges = Probabilistic(prob_matrix)
        nodes = Unipartite(Symbol.(data.species))
        N = SpeciesInteractionNetwork(nodes, edges)
        return randomdraws(N)
    end

    return network
    
end
