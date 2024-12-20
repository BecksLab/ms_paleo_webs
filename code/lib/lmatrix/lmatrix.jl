using DataFrames
using SpeciesInteractionNetworks

"""
  lmatrix(
    species::Any,
    bodymass::Vector{Float64},
    is_producer::Vector{Bool};
    Ropt::Float64,
    γ::Float64,
    threshold::Float64,
)

    A port of the ATNr L matrix function from Gauzens et al. (2023). 
    Interactions are determined by allometric rules and a Ricker function 
    defined by `Ropt` and `γ` and returns a probabilistic 
    `SpeciesInteractionNetwork`. `Ropt`, `γ`, and `threshold` use the ATNr
    defaults.

    Requires a species list, vector of bodymass, and if the species is a
    producer or not.
    
    #### References
    
    Rohr, R.P., H. Scherer, P. Kehrli, C. Mazza, and L-F. Bersier. 2010. 
    “Modeling Food Webs: Exploring Unexplained Structure Using Latent 
    Traits.” The American Naturalist 176 (2): 170–77. 
    https://doi.org/10.1086/653667.

    Yeakel, J.D., M.M. Pires, L. Rudolf, N.J. Dominy, P.L. Koch, P.R. 
    Guimarães, and T. Gross. 2014. “Collapse of an Ecological Network in 
    Ancient Egypt.” PNAS 111 (40): 14472–77. 
    https://doi.org/10.1073/pnas.1408471111.

"""
function lmatrix(
    species::Any,
    bodymass::Vector{Float64},
    is_producer::Vector{Bool};
    Ropt::Float64 = 100.0,
    γ::Float64 = 2.0,
    threshold::Float64 = 0.01,
)

    S = length(species)

    prob_matrix = zeros(AbstractFloat, (S, S))
    
    for i = 1:S
        # only assess if consumer is not a producer
        if !is_producer[i]
            for j = 1:S
                l = bodymass[i] / (bodymass[j] * Ropt)
                L = (l * exp(1 - l)) ^γ
                if L ≥ threshold
                    prob_matrix[i, j] = L
                end
            end
        end
    end

    edges = Probabilistic(prob_matrix)
    nodes = Unipartite(Symbol.(species))
    return SpeciesInteractionNetwork(nodes, edges)
end
