using Distributions, Random

# Define bounds for each size class
size_bounds = Dict(
    "primary" => (0.001, 0.09),
    "tiny" => (0.1, 10.0),
    "small" => (10.0, 50.0),
    "medium" => (50.0, 100.0),
    "large" => (100.0, 300.0),
    "very_large" => (300.0, 500.0),
    "gigantic" => (500.0, 700.0)
)

# Function to sample body sizes
function sample_body_size(size_class::String; method::String="uniform")
    lower, upper = size_bounds[size_class]
    
    if method == "uniform"
        return rand(Uniform(lower, upper))
    elseif method == "lognormal"
        # Estimate mu and sigma from bounds to keep ~95% within bounds
        μ = log((lower + upper)/2)  # approximate median
        σ = 0.5                     # adjust to control spread
        val = rand(LogNormal(μ, σ))
        # truncate to bounds
        return clamp(val, lower, upper)
    elseif method == "truncated_lognormal"
        # Use truncated lognormal
        μ = log((lower + upper)/2)
        σ = 0.5
        dist = truncated(LogNormal(μ, σ), lower, upper)
        return rand(dist)
    else
        error("Unknown sampling method: $method")
    end
end