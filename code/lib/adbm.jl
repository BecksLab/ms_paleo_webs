# For implementing the ADBM model

"""
  adbm_parameters(...)
  Returns the bioenergetic constants needed for the ADBM. 
  Defaults are based on empirical fits from Petchey et al. (2008).
"""
function adbm_parameters(
    df::DataFrame,
    bodymass::Vector{Float64};
    e::Float64 = 1.0,        # Assimilation efficiency
    a_adbm::Float64 = 0.0189, # Scaling constant for attack rate
    ai::Float64 = -0.491,    # Allometric exponent for prey size
    aj::Float64 = -0.465,    # Allometric exponent for consumer size
    b::Float64 = 0.401,      # Scaling constant for handling time
    h_adbm::Float64 = 1.0,
    hi::Float64 = 1.0,
    hj::Float64 = 1.0,
    n::Float64 = 1.0,
    ni::Float64 = -0.75,
    Hmethod::Symbol = :ratio,
    Nmethod::Symbol = :original,
)

    parameters = Dict{Symbol,Any}(
        :e => e,
        :a_adbm => a_adbm,
        :ai => ai,
        :aj => aj,
        :b => b,
        :h_adbm => h_adbm,
        :hi => hi,
        :hj => hj,
        :n => n,
        :ni => ni,
    )

    # ... [Method checks for H and N omitted] ...

    # Flag producers to ensure they don't "hunt" in the simulation
    parameters[:is_producer] = df.feeding .== "producer"
    parameters[:M] = bodymass
    return parameters
end

"""
  _get_adbm_terms(S, parameters, biomass)
  Calculates the intermediate matrices for Energy (E), Attack Rate (λ), and Handling Time (H).
"""
function _get_adbm_terms(S::Int64, parameters::Dict{Symbol,Any}, biomass::Vector{Float64})
    M = parameters[:M]
    
    # λ_ij: Attack rate of consumer j on prey i
    # Calculated as a function of the body masses of both species
    λ = [parameters[:a_adbm] * (M[i]^parameters[:ai]) * (M[j]^parameters[:aj]) for i = 1:S, j = 1:S]
    
    # E_i: Energy content of prey i (assumed proportional to its body mass)
    E = [parameters[:e] * M[i] for i = 1:S]
    
    # H_ij: Handling time (how long it takes consumer j to process prey i)
    H = [parameters[:b] * (M[i]^parameters[:hi]) * (M[j]^parameters[:hj]) for i = 1:S, j = 1:S]

    # N_i: Prey density (resource availability)
    N = [parameters[:n] * (M[i]^parameters[:ni]) for i = 1:S]

    return Dict(:E => E, :λ => λ, :H => H, :N => N)
end

"""
  adbmmodel(df, parameters, biomass)
  The main function that builds the adjacency matrix. 
  It ranks potential prey by profitability (Energy / Handling Time) and 
  adds them to the diet until the intake rate is maximized.
"""
function adbmmodel(df::DataFrame, parameters::Dict{Symbol,Any}, biomass::Vector{Float64})
    S = nrow(df)
    adbmMAT = zeros(Bool, (S, S))
    
    # Calculate bioenergetic terms
    adbmTerms = _get_adbm_terms(S, parameters, biomass)
    E = adbmTerms[:E]
    λ = adbmTerms[:λ]
    H = adbmTerms[:H]
    N = adbmTerms[:N]

    for j = 1:S
        # Only evaluate consumers
        if !parameters[:is_producer][j]
            
            # 1. Rank all potential prey by profitability (E_i / H_ij)
            profitability = [E[i] / H[i, j] for i = 1:S]
            # Create a sequence of indices from most to least profitable
            idx = sortperm(profitability, rev = true)

            # 2. Iteratively add prey to the diet
            # Optimal Foraging Theory: Intake rate increases as you add profitable prey,
            # but eventually drops as you add low-quality items that take too long to handle.
            current_intake = 0.0
            for i in idx
                # Calculate potential intake rate with the new item included
                # Formula: (Sum of Energy * Attack * Density) / (1 + Sum of Handling * Attack * Density)
                numerator = sum([λ[k, j] * E[k] * N[k] for k in idx[1:findfirst(x -> x == i, idx)]])
                denominator = 1 + sum([λ[k, j] * H[k, j] * N[k] for k in idx[1:findfirst(x -> x == i, idx)]])
                new_intake = numerator / denominator

                # If the new intake rate is better (or equal), add the link
                if new_intake >= current_intake
                    adbmMAT[i, j] = true
                    current_intake = new_intake
                else
                    # Once intake rate starts to drop, the "optimal diet" is found
                    break
                end
            end
        end
    end

    # Return as a SpeciesInteractionNetwork object
    return SpeciesInteractionNetwork(Unipartite(Symbol.(df.species)), Binary(adbmMAT))
end