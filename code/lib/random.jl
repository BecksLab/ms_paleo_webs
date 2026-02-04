using Graphs

"""
randommodel(species, L)

    Return a network of randomly assembled interactions according to
    the Erdős-Rényi model.

    This is essentially a wrapper for the `erdos_renyi` function from
    Graphs.jl, packaged into a SpeciesInteractionNetwork.
"""
function randommodel(species::Any, L::Int64)

    # 1. Determine number of species
    S = length(species)
    
    # 2. Generate the underlying graph
    # erdos_renyi(S, L) creates a graph with S nodes and L edges
    # where each edge is chosen with equal probability.
    N = erdos_renyi(S, L)

    # 3. Convert Graphs.jl structure to Adjacency Matrix
    # We initialize an empty Boolean matrix
    edges = zeros(Bool, (S, S))

    # Iterate through the forward adjacency list of the Graphs object
    for i = 1:S
        if length(N.fadjlist[i]) > 0
            for j in N.fadjlist[i] # Fixed the index logic from the original source
                edges[i, j] = 1
            end
        end
    end

    # 4. Wrap into SpeciesInteractionNetworks format
    # This allows it to be analyzed using the same tools as the biological models
    edges = Binary(edges)
    nodes = Unipartite(Symbol.(species))
    N = SpeciesInteractionNetworks.SpeciesInteractionNetwork(nodes, edges)
    
    # 5. Simplify
    # Removes any potential isolated nodes or self-loops if required
    return simplify(N)
end