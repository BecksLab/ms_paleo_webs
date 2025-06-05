using DataFrames
using SpeciesInteractionNetworks

"""
  bmratio(
    species::Any,
    bodymass::Vector{Float64};
    α::Float64,
    β::Float64,
    γ::Float64,
)

    Implements the Bodymass-ratio model as per Rohr et al. (2010), 
    however does not include the latent trait variables. α, β, and γ 
    parameters are set to those used by Yeakel (2014) and so only
    requires a species list and their bodymass and returns a 
    `SpeciesInteractionNetwork`.
    
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
function bmratio(
    species::Any,
    bodymass::Vector{Float64};
    α::Float64 = 1.41,
    β::Float64 = 3.73,
    γ::Float64 = 1.87,
    )

    S = length(species)

    prob_matrix = zeros(AbstractFloat, (S, S))

    for i = 1:S
        for j = 1:S
            MR = bodymass[i] / bodymass[j]
            p = exp(α + β * log(MR) + γ * (log(MR))^2)

           if p / (1 - p) >= 0.0
                prob_matrix[i, j] = p / (1 + p)
            end
        end
    end

    # make probabilistic
    # prob_matrix = prob_matrix ./ maximum(prob_matrix)

    edges = Probabilistic(prob_matrix)
    nodes = Unipartite(Symbol.(species))
    return SpeciesInteractionNetwork(nodes, edges)
end
