using Graphs

"""
randommodel(S::Int64, L::Int64)

    Return a network of randomly assembled interactions according to
    the Erdős-Rényi model.

    This is essentially a wrapper for the `erdos_renyi` function from
    Graphs.jl and just packes it into a SpeciesInteraction network

    #### References

    Erdős, Paul, and Alfréd Rényi. 1959. “On Random Graphs I.” 
    Publicationes Mathematicae. https://doi.org/10.5486/PMD.1959.6.3-4.12.

    Graphs.jl TODO
"""
function randommodel(species::Any, L::Int64)

    S = length(species)
    N = erdos_renyi(S, L)

    # empty matrix
    edges = zeros(Bool, (S, S))

    for i in 1:S
        if length(N.fadjlist[i]) > 0
            for j in eachindex(N.fadjlist[i])
                edges[i,j] = 1
            end
        end
    end

    edges = Binary(edges)
    nodes = Unipartite(Symbol.(species))
    return SpeciesInteractionNetworks.SpeciesInteractionNetwork(nodes, edges)
end