#library
library(here)
library(igraph)
library(tidyverse)

# set path to code sub dir
setwd(here())

df <- read.csv(here("data/raw/feeding_rules.csv")) %>%
  mutate_all(~str_replace_all(., "nonmotile", "non_motile")) %>%
  mutate_all(~str_replace_all(., "-", "_")) 

traits_classes <- df %>%
  select(trait_type_resource) %>%
  distinct() %>%
  pull()

#### Maximal traits ####

# write rules .csv
write.csv(df, here("data/clean/feeding_rules/feeding_rules_maximal.csv"),
          row.names = FALSE)

trait_networks <- vector(mode = "list", length = 4)

for (i in seq_along(traits_classes)) {
  
  df_sub <- df %>%
    filter(trait_type_resource == traits_classes[i]) %>%
    select(trait_resource, trait_consumer) %>%
    # remove primary node
    filter(!str_detect(trait_consumer, "primary")) %>%
    filter(!str_detect(trait_resource, "primary"))
  
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
  mutate(trait_resource = case_when(trait_resource == "shallow-infaunal" ~ "infaunal",
                                    trait_resource == "deep-infaunal" ~ "infaunal",
                                    str_detect(trait_resource, "^.*epifaunal.*$") ~ "epifaunal",
                                    trait_resource == "non-motile_attached" ~ "attached",
                                    trait_resource == "non-motile_byssate" ~ "attached",
                                    str_detect(trait_resource, "^.*deposit.*$") ~ "herbivore",
                                    str_detect(trait_resource, "^.*suspension.*$") ~ "herbivore",
                                    str_detect(trait_resource, "grazer_herbivore") ~ "herbivore",
                                    TRUE ~ as.character(trait_resource)),
         trait_consumer = case_when(trait_consumer == "shallow-infaunal" ~ "infaunal",
                                    trait_consumer == "deep-infaunal" ~ "infaunal",
                                    str_detect(trait_consumer, "^.*epifaunal.*$") ~ "epifaunal",
                                    trait_consumer == "non-motile_attached" ~ "attached",
                                    trait_consumer == "non-motile_byssate" ~ "attached",
                                    str_detect(trait_consumer, "^.*deposit.*$") ~ "herbivore",
                                    str_detect(trait_consumer, "^.*suspension.*$") ~ "herbivore",
                                    str_detect(trait_consumer, "grazer_herbivore") ~ "herbivore",
                                    TRUE ~ as.character(trait_consumer))) %>%
  distinct()

# write rules .csv
write.csv(df_min, here("data/clean/feeding_rules/feeding_rules_minimum.csv"),
          row.names = FALSE)

trait_networks <- vector(mode = "list", length = 4)

for (i in seq_along(traits_classes)) {
  
  df_sub <- df_min %>%
    filter(trait_type_resource == traits_classes[i]) %>%
    select(trait_resource, trait_consumer) %>%
    # remove primary node
    filter(!str_detect(trait_consumer, "primary")) %>%
    filter(!str_detect(trait_resource, "primary"))
  
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

# write rules .csv
write.csv(df_scav, here("data/clean/feeding_rules/feeding_rules_noscavs.csv"),
          row.names = FALSE)

trait_networks <- vector(mode = "list", length = 4)

for (i in seq_along(traits_classes)) {
  
  df_sub <- df_scav %>%
    filter(trait_type_resource == traits_classes[i]) %>%
    select(trait_resource, trait_consumer) %>%
    # remove primary node
    filter(!str_detect(trait_consumer, "primary")) %>%
    filter(!str_detect(trait_resource, "primary"))
  
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

