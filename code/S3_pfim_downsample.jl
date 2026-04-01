using CSV
using DataFrames
using Distributions
using JLD2
using pfim
using Plots
using SpeciesInteractionNetworks
using Statistics
using StatsPlots
import Random

# Reproducibility
Random.seed!(66)

# internals
include("lib/internals.jl");

# -----------------------------
# LOAD DATA
# -----------------------------
matrix_names = readdir("../data/raw")
matrix_names = matrix_names[occursin.(r"^.*Guilds.*$", matrix_names)]

feeding_rules = DataFrame(CSV.File("../data/raw/feeding_rules.csv"))
size_classes = DataFrame(CSV.File("../data/raw/size_classes.csv"))

# -----------------------------
# PARAMETERS
# -----------------------------
n_reps = 100

# Range of downsampling values to test
downsample_vals = [5.0, 10.0, 20.0, 30.0, 50.0, 100.0]

# -----------------------------
# OUTPUT STORAGE
# -----------------------------
networks = DataFrame(
    time = String[],
    n_rep = Int[],
    downsample_y = Float64[],
    network = Any[],
    richness = Int[],
    links = Int[],
    connectance = Float64[]
)

# -----------------------------
# SIMULATION LOOP
# -----------------------------
for j = 1:n_reps

    # --- Loop through communities ---
    for file_name in matrix_names

        str_cats = split(file_name, r"_")

        df = DataFrame(CSV.File(joinpath("../data/raw/", file_name)))
        select!(df, [:Guild, :motility, :tiering, :feeding, :size])
        rename!(df, :Guild => :species)

        filter!(:species => x -> x != "BASAL_NODE", df)

        # --- First: generate baseline (no downsampling) ---
N_meta = pfim.PFIM(df, feeding_rules; downsample = false)

if richness(N_meta) > 0
    push!(networks, (
        str_cats[1],
        j,
        0.0,  # 👈 use 0 to indicate "no downsampling"
        N_meta,
        richness(N_meta),
        links(N_meta),
        connectance(N_meta)
    ))
end

# --- Then: downsampled versions ---
for y_val in downsample_vals

    N = pfim.PFIM(df, feeding_rules; y = y_val, downsample = true)

    if richness(N) > 0
        push!(networks, (
            str_cats[1],
            j,
            y_val,
            N,
            richness(N),
            links(N),
            connectance(N)
        ))
    end
end
    end
end

# -----------------------------
# SAVE OUTPUT
# -----------------------------
save_object("../data/processed/pfim_downsampling_sensitivity.jlds", networks)

function edge_set(N)
    A = _get_matrix(N)
    return Set(findall(!iszero, A))
end

function jaccard_similarity(N1, N2)
    E1 = edge_set(N1)
    E2 = edge_set(N2)
    return length(intersect(E1, E2)) / length(union(E1, E2))
end

# Store results
similarity_df = DataFrame(
    downsample_y = Float64[],
    similarity = Float64[]
)

for g in groupby(networks, [:time, :n_rep])

    # baseline network (y = 0)
    base = filter(:downsample_y => ==(0.0), g)

    if nrow(base) == 0
        continue
    end

    N_base = base.network[1]

    for row in eachrow(g)

        if row.downsample_y == 0.0
            continue
        end

        sim = jaccard_similarity(N_base, row.network)

        push!(similarity_df, (row.downsample_y, sim))
    end
end


CSV.write("../data/processed/pfim_downsampling_similarity.csv", similarity_df)