using DataFrames
using SpeciesInteractionNetworks

"""
    _PFIM_network(trait_data::DataFrame, feeding_rules::DataFrame)

    Internal function that constructs a network for a community
    when only categorical traits are given.

"""
function _PFIM_network(trait_data::DataFrame, feeding_rules::DataFrame)

    S = nrow(trait_data)
    int_matrix = zeros(Bool, (S, S))
    
    for cons in 1:nrow(trait_data)
        for res in 1:nrow(trait_data)
            consumer = trait_data[cons,:]
            resource = trait_data[res,:]
    
            # keep record if rule is met or not
            tally = 0
    
            for i in Symbol.(unique(collect(feeding_rules.trait_type_resource)))
                consumer_trait = consumer[i]
                resource_trait = resource[i]
                resources = filter(:trait_consumer => x -> x == consumer_trait, feeding_rules).trait_resource
                if resource_trait ∈ resources
                    tally = tally + 1
                end
            end
    
            # only add link if all 4 rules are met
            if tally == 4
                int_matrix[cons, res] = 1
            end
    
        end
    end
    
    nodes = Unipartite(Symbol.(trait_data.species))
    edges = Binary(int_matrix)
    network = SpeciesInteractionNetwork(nodes, edges)
    return network, int_matrix
end

"""
    _downsample(network, matrix, y)

    Internal function to downsample a network based on species link 
    distributions.

    #### References
    Roopnarine, Peter D. 2006. “Extinction Cascades and Catastrophe in
    Ancient Food Webs.” Paleobiology 32 (1): 1-19. 
    https://www.jstor.org/stable/4096814.

"""
function _downsample(network, matrix, y)
    
    link_dist = zeros(Float64, richness(network))
    spp = species(network)
    S = richness(network)

            # get link distribution
            for i in eachindex(spp)
                sp = spp[i]
                r = generality(network, sp)
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
            # make probabalistic
            prob_matrix = prob_matrix ./ maximum(prob_matrix)
    
            edges = Probabilistic(prob_matrix)
            nodes = Unipartite(spp)
            N = SpeciesInteractionNetwork(nodes, edges)
            return randomdraws(N)
end

"""
   PFIM(trait_data::DataFrame, feeding_rules::DataFrame; y::Float64 = 2.5, downsample::Bool = true)

    Takes a data frame and implements the feeding rules to determine the
    feasibility of links between species. As well as applying the link
    distribution downsampling approach.
    
    #### References
    
    Shaw, Jack O., Alexander M. Dunhill, Andrew P. Beckerman, Jennifer A.
    Dunne, and Pincelli M. Hull. 2024. “A Framework for Reconstructing 
    Ancient Food Webs Using Functional Trait Data.” 
    https://doi.org/10.1101/2024.01.30.578036.
"""
function PFIM(trait_data::DataFrame, feeding_rules::DataFrame; y::Float64 = 2.5, downsample::Bool = true)

    # data checks
    for (i, v) in enumerate(["species", "motility", "tiering", "feeding", "size"])
        if v ∉ names(trait_data)
            error("Missing $(v) variable as a column in DataFrame, add or rename")
        end
    end

    network, matrix = _PFIM_network(trait_data, feeding_rules)

    # downsampling protocol
    if downsample == true
        network = _downsample(network, matrix, y)
    end

    return network

end
