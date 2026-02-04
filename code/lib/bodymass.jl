using DataFrames
using SpeciesInteractionNetworks

"""
  bmratio(species, bodymass; α, β, γ)

    Determines links based on the log-ratio of consumer to prey body mass.
    Returns a Unipartite Binary SpeciesInteractionNetwork.
    
    α, β, and γ are the statistical parameters that define the 
    location and breadth of the feeding niche.
"""
function bmratio(
    species::Any,
    bodymass::Vector{Float64};
    α::Float64 = 1.41,  # Center of the feeding niche (optimum ratio)
    β::Float64 = 3.73,  # Width/Breadth of the feeding niche
    γ::Float64 = -1.90, # Scaling constant
)

    S = length(species)
    # Convert body mass to log10 for allometric scaling
    logM = log10.(bodymass)
    
    # Initialize adjacency matrix
    prob_matrix = zeros(Bool, (S, S))

    for i = 1:S
        # Potential consumer loop
        for j = 1:S
            
            # 1. Calculate the log-ratio of masses
            # This represents the relative size difference in log-space
            x_ij = logM[j] - logM[i]

            # 2. Probability Function (Rohr/Yeakel Logic)
            # This calculates the 'energy' or 'fit' of the interaction.
            # It creates a peak probability at the ratio defined by α.
            # Note: The original model often includes latent traits (ε); 
            # here we assume ε = 0.
            
            # Formula: - ( (x_ij - α) / β )^2
            p_val = -((x_ij - α) / β)^2

            # 3. Deterministic Threshold
            # In this version, if the 'fit' is greater than γ, a link is formed.
            # This creates a "diet envelope" around the optimal size ratio.
            if p_val > γ
                prob_matrix[i, j] = true
            end
        end
    end

    # Wrap as a Network object
    return SpeciesInteractionNetwork(Unipartite(Symbol.(species)), Binary(prob_matrix))
end