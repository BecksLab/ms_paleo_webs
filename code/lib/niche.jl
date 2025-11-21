"""
randommodel(S::Int64, L::Int64)

    Wrapper of the structuralmodel(NicheModel) that allows us to assign 'real' species 
    as opposed to 1:n named species
"""
function nichemodel(species::Any, connectance::Float64)

    S = length(species)
    N = structuralmodel(NicheModel, S, connectance)

    A = _get_matrix(N)

    edges = Binary(A)
    nodes = Unipartite(Symbol.(species))
    N = SpeciesInteractionNetworks.SpeciesInteractionNetwork(nodes, edges)
    return simplify(N)
end