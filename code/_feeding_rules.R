#library
library(here)
library(igraph)

# set path to code sub dir
setwd(here())

df <- read.csv(here("data/raw/feeding_rules.csv"))

traits_classes <- df %>%
  select(trait_type_resource) %>%
  distinct() %>%
  pull()

#### Maximal traits ####

trait_networks <- vector(mode = "list", length = 4)

for (i in seq_along(traits_classes)) {
  
  df_sub <- df %>%
    filter(trait_type_resource == traits_classes[i]) %>%
    select(trait_resource, trait_consumer)
  
  # create data:
  links <- data.frame(
    source = df_sub$trait_resource,
    target = df_sub$trait_consumer
  )
  
  # create the network object
  trait_networks[[i]] <- graph_from_data_frame(d=links, directed=T) 
  
  
}

png(file = paste0(here('figures/feeding_rules'),"/maximal.png"),
    width = 1100, height = 900, units = "px", res = 100)
# plot it
par(mfrow=c(2,2), mar=c(1,1,1,1))
plot(trait_networks[[1]], layout = layout.circle, main = "Tiering", edge.arrow.size =.3)
plot(trait_networks[[2]], layout = layout.circle, main = "Motility", edge.arrow.size =.3)
plot(trait_networks[[3]], layout = layout.circle, main = "Feeding", edge.arrow.size =.3)
plot(trait_networks[[4]], layout = layout.circle, main = "Size", edge.arrow.size =.3)
dev.off()

#### Minimal traits ####

df_min <- df %>%
  mutate(trait_resource = case_when(trait_resource == "epifaunal_erect" ~ "epifaunal",
                                    trait_resource == "epifaunal_surficial" ~ "epifaunal",
                                    trait_resource == "shallow-infaunal" ~ "infaunal",
                                    trait_resource == "deep-infaunal" ~ "infaunal",
                                    trait_resource == "non-motile_attached" ~ "attached",
                                    trait_resource == "non-motile_byssate" ~ "attached",
                                    trait_resource == "deposit_surficial" ~ "herbivore",
                                    trait_resource == "deposit_mining" ~ "herbivore",
                                    trait_resource == "deposit_mining_chemosymbiotic" ~ "herbivore",
                                    trait_resource == "grazer_herbivore" ~ "herbivore",
                                    trait_resource == "suspension" ~ "herbivore",
                                    trait_resource == "suspension_chemosymbiotic" ~ "herbivore",
                                    .default = trait_resource),
         trait_consumer = case_when(trait_consumer == "epifaunal_erect" ~ "epifaunal",
                                    trait_consumer == "epifaunal_surficial" ~ "epifaunal",
                                    trait_consumer == "shallow-infaunal" ~ "infaunal",
                                    trait_consumer == "deep-infaunal" ~ "infaunal",
                                    trait_consumer == "non-motile_attached" ~ "attached",
                                    trait_consumer == "non-motile_byssate" ~ "attached",
                                    trait_consumer == "deposit_surficial" ~ "herbivore",
                                    trait_consumer == "deposit_mining" ~ "herbivore",
                                    trait_consumer == "deposit_mining_chemosymbiotic" ~ "herbivore",
                                    trait_consumer == "grazer_herbivore" ~ "herbivore",
                                    trait_consumer == "suspension" ~ "herbivore",
                                    trait_consumer == "suspension_chemosymbiotic" ~ "herbivore",
                                    .default = trait_consumer)) %>%
  distinct()

trait_networks <- vector(mode = "list", length = 4)

for (i in seq_along(traits_classes)) {
  
  df_sub <- df_min %>%
    filter(trait_type_resource == traits_classes[i]) %>%
    select(trait_resource, trait_consumer)
  
  # create data:
  links <- data.frame(
    source = df_sub$trait_resource,
    target = df_sub$trait_consumer
  )
  
  # create the network object
  trait_networks[[i]] <- graph_from_data_frame(d=links, directed=T) 
  
  
}

png(file = paste0(here('figures/feeding_rules'),"/minimal.png"),
    width = 1100, height = 900, units = "px", res = 100)
# plot it
par(mfrow=c(2,2), mar=c(1,1,1,1))
plot(trait_networks[[1]], layout = layout.circle, main = "Tiering", edge.arrow.size =.3)
plot(trait_networks[[2]], layout = layout.circle, main = "Motility", edge.arrow.size =.3)
plot(trait_networks[[3]], layout = layout.circle, main = "Feeding", edge.arrow.size =.3)
plot(trait_networks[[4]], layout = layout.circle, main = "Size", edge.arrow.size =.3)
dev.off()

#### No scavengers/parasites ####

df_scav <- df_min %>%
  filter(!str_detect(trait_consumer, "scavenger")) %>%
  filter(!str_detect(trait_consumer, "parasitic")) %>%
  filter(!str_detect(trait_resource, "scavenger")) %>%
  filter(!str_detect(trait_resource, "parasitic")) %>%
  distinct()

trait_networks <- vector(mode = "list", length = 4)

for (i in seq_along(traits_classes)) {
  
  df_sub <- df_scav %>%
    filter(trait_type_resource == traits_classes[i]) %>%
    select(trait_resource, trait_consumer)
  
  # create data:
  links <- data.frame(
    source = df_sub$trait_resource,
    target = df_sub$trait_consumer
  )
  
  # create the network object
  trait_networks[[i]] <- graph_from_data_frame(d=links, directed=T) 
  
  
}

png(file = paste0(here('figures/feeding_rules'),"/no_scav.png"),
    width = 1100, height = 900, units = "px", res = 100)
# plot it
par(mfrow=c(2,2), mar=c(1,1,1,1))
plot(trait_networks[[1]], layout = layout.circle, main = "Tiering", edge.arrow.size =.3)
plot(trait_networks[[2]], layout = layout.circle, main = "Motility", edge.arrow.size =.3)
plot(trait_networks[[3]], layout = layout.circle, main = "Feeding", edge.arrow.size =.3)
plot(trait_networks[[4]], layout = layout.circle, main = "Size", edge.arrow.size =.3)
dev.off()

