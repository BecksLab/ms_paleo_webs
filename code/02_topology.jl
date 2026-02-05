using CSV
using DataFrames
using JLD2

# Load internal helper functions, specifically the `_network_summary` utility
include("lib/internals.jl")

# Data Ingestion: Load the serialized JLD2 file containing the synthetic networks [cite: 15]
networks = load_object("../data/processed/networks.jlds")

# Pre-allocate a DataFrame to store topological features [cite: 15]
# This defines the schema for the structural analysis of each network
topology = DataFrame(
    model = String[],
    time = Any[],
    n_rep = Any[],
    richness = Int64[],      # Number of nodes (species/guilds)
    connectance = Float64[],   # Realized fraction of possible links
    diameter = Int64[],        # Longest shortest path between any two nodes
    complexity = Float64[],    # Often related to link density or eigenvalue properties
    trophic_level = Float64[], # Mean vertical position in the food web
    distance = Float64[],      # Mean path length between nodes
    generality = Float64[],    # Average number of resources per consumer
    vulnerability = Float64[], # Average number of consumers per resource
    redundancy = Float64[],    # Overlap in ecological roles/links
    S1 = Float64[],            # Structural motif 1
    S2 = Float64[],            # Structural motif 2
    S4 = Float64[],            # Structural motif 4 
    S5 = Float64[],            # Structural motif 5 
);

# Analysis Loop: Process every generated network from the previous step 
for i = 1:nrow(networks)

    # Calculate metrics using the internal helper function 
    # This likely calls functions from SpeciesInteractionNetworks.jl internally
    d = _network_summary(networks.network[i])

    # Re-attach metadata from the original 'networks' dataframe 
    d[:model] = networks.model[i]
    d[:n_rep] = networks.n_rep[i]
    d[:time] = networks.time[i]

    # Append the resulting dictionary of metrics to the topology table 
    push!(topology, d)
end

# Data Export: Save the final table as a CSV for easy use in plotting or R/Python 
CSV.write("../data/processed/topology.csv", topology)