include(joinpath("types.jl"))

# Define rules for each trait type

# feeding rules
feeding_rules(consumer::T1, resource::T2) where {T1<:feeding,T2<:feeding} = 0
feeding_rules(consumer::carnivore, resource::T) where {T<:feeding} = 1
feeding_rules(consumer::scavenger, resource::T) where {T<:feeding} = 1
feeding_rules(consumer::parasitic, resource::T) where {T<:feeding} = 1
feeding_rules(consumer::parasitic, resource::parasitic) = 1

# motility rules
motility_rules(consumer::T1, resource::T2) where {T1<:motility,T2<:motility} = 0
motility_rules(consumer::fast_moving, resource::T) where {T<:motility} = 1
motility_rules(consumer::fast_moving, resource::fast_moving) = 1
motility_rules(consumer::fast_moving, resource::non_motile) = 0
motility_rules(consumer::slow_moving, resource::T) where {T<:motility} = 1
motility_rules(consumer::slow_moving, resource::fast_moving) = 0
motility_rules(consumer::slow_moving, resource::slow_moving) = 1
motility_rules(consumer::slow_moving, resource::non_motile) = 0
motility_rules(consumer::facultative, resource::facultative) = 1
motility_rules(consumer::facultative, resource::non_motile) = 1
motility_rules(consumer::grazer_carnivore, resource::suspension) = 1

# tiering rules
tiering_rules(consumer::T1, resource::T2) where {T1<:tier,T2<:tier} = 0
tiering_rules(consumer::T, resource::T) where {T<:tier} = 1
tiering_rules(consumer::nektonic, resource::epifaunal_erect) = 1
tiering_rules(consumer::nektonic, resource::epifaunal_surficial) = 1
tiering_rules(consumer::epifaunal_erect, resource::nektonic) = 1
tiering_rules(consumer::epifaunal_erect, resource::epifaunal_surficial) = 1
tiering_rules(consumer::epifaunal_surficial, resource::epifaunal_erect) = 1
tiering_rules(consumer::epifaunal_surficial, resource::semi_infaunal) = 1
tiering_rules(consumer::semi_infaunal, resource::epifaunal_surficial) = 1
tiering_rules(consumer::semi_infaunal, resource::shallow_infaunal) = 1
tiering_rules(consumer::semi_infaunal, resource::deep_infaunal) = 1
tiering_rules(consumer::shallow_infaunal, resource::semi_infaunal) = 1
tiering_rules(consumer::shallow_infaunal, resource::deep_infaunal) = 1
tiering_rules(consumer::deep_infaunal, resource::semi_infaunal) = 1
tiering_rules(consumer::deep_infaunal, resource::shallow_infaunal) = 1

# size rules
size_rules(consumer::T1, resource::T2) where {T1<:sizes,T2<:sizes} = 0
size_rules(consumer::T, resource::T) where {T<:sizes} = 1
size_rules(consumer::very_large, resource::T) where {T<:sizes} = 1
size_rules(consumer::large, resource::medium) = 1
size_rules(consumer::large, resource::small) = 1
size_rules(consumer::large, resource::tiny) = 1
size_rules(consumer::medium, resource::small) = 1
size_rules(consumer::medium, resource::tiny) = 1
size_rules(consumer::small, resource::tiny) = 1
