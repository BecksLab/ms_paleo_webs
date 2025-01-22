library(here)
library(tidyverse)

# set path to code sub dir
setwd(here("code"))

# clean traits data
trait_files <- list.files(path = "../data/raw", pattern = "traits.csv", full.names = TRUE)

for (i in seq_along(trait_files)) {
  
  df <- read.csv(trait_files[i])   %>%
    # standardise strings
    # mutate(species = str_replace_all(species, "[:punct:]", "")) %>% 
    # then all whitespace
    mutate(species = str_replace_all(species, " ", "_")) %>%
    # standardise strings
    mutate_all(~str_replace_all(., "-", "_"))  %>%
    # change motiliy class - we can remove this later
    mutate(motility = case_when( 
      motility == "nonmotile" ~ "non_motile",
      TRUE ~ as.character(motility))) %>% 
    # removes species with no size data
    filter(species != "Hemigordius_baoqingensis") %>%
    # assign size class to species
    mutate(size = case_when(
      species == "Hemigordius baoqingensis" ~ "medium",
      TRUE ~ as.character(size)))
  
  # write as clean data
  write.csv(df, str_replace(trait_files[i], "raw.", "clean/trait_maximal/"),
            row.names = FALSE)
  
  df  %>%
    # combine some trait classes
    mutate(tiering = case_when(tiering == "shallow_infaunal" ~ "infaunal",
                               tiering == "deep_infaunal" ~ "infaunal",
                               str_detect(tiering, "^.*epifaunal.*$") ~ "epifaunal",
                               TRUE ~ as.character(tiering)),
           feeding = case_when(str_detect(feeding, "^.*deposit.*$") ~ "herbivore",
                               str_detect(feeding, "^.*suspension.*$") ~ "herbivore",
                               str_detect(feeding, "grazer_herbivore") ~ "herbivore",
                               TRUE ~ as.character(feeding)),
           motility = case_when(motility == "non-motile_attached" ~ "attached",
                                motility == "non-motile_byssate" ~ "attached",
                                TRUE ~ as.character(motility)),
           size = case_when(str_detect(size, "^.*large.*$") ~ "large",
                            str_detect(size, "^.*medium.*$") ~ "medium",
                            str_detect(size, "^.*small.*$") ~ "small",
                            str_detect(size, "^.*tiny.*$") ~ "tiny")) %>%
    # write as clean data
    write.csv(., str_replace(trait_files[i], "raw.", "clean/trait_minimum/"),
              row.names = FALSE)
}

# clean size data

size_files <- list.files(path = "../data/raw", pattern = "size.csv", full.names = TRUE)

for (i in seq_along(size_files)) {
  
  df <- read.csv(size_files[i])  %>%
    # standardise strings
    # mutate(species = str_replace_all(species, "[:punct:]", "")) %>% 
    # then all whitespace
    mutate(species = str_replace_all(species, " ", "_")) %>%
    # standardise strings
    mutate(species = str_replace_all(species, "-", "_")) %>%
    select(species, Height..interpreted., Length..interpreted.)  %>%
    group_by(species)  %>%
    summarise(height = mean(Height..interpreted.),
              length = mean(Length..interpreted.))  %>% 
    reframe(
      species = species,
      size = length * height)
  
  # write as clean data
  write.csv(df, str_replace(size_files[i], "raw.", "clean/size/"),
            row.names = FALSE)
}

# create simple traits df for extinction sims

list.files(path = "../data/clean/trait_maximal", pattern = ".csv", full.names = TRUE) %>% 
  lapply(read_csv) %>% 
  bind_rows %>%
  select(-feeding) %>%
  mutate(motility = case_when(str_detect(motility, "^.*non_motile.*$") ~ "non_motile",
                              motility == "slow_moving" ~ "slow",
                              motility == "fast_moving" ~ "fast",
                              TRUE ~ as.character(motility)),
         tiering = case_when(str_detect(tiering, "^.*infaunal.*$") ~ "infaunal",
                             str_detect(tiering, "^.*epifaunal.*$") ~ "epifaunal",
                             TRUE ~ as.character(tiering))
  ) %>%
  group_by(species, motility, tiering, size) %>% 
  distinct() %>%
  write.csv(., "../data/clean/extinction/traits.csv",
            row.names = FALSE)

# create simple traits df for extinction sims - but for the trophic species

list.files(path = "../data/clean/trophic_maximal", pattern = ".csv", full.names = TRUE) %>% 
  lapply(read_csv) %>% 
  bind_rows %>%
  select(-feeding) %>%
  mutate(motility = case_when(str_detect(motility, "^.*non_motile.*$") ~ "non_motile",
                              motility == "slow_moving" ~ "slow",
                              motility == "fast_moving" ~ "fast",
                              TRUE ~ as.character(motility)),
         tiering = case_when(str_detect(tiering, "^.*infaunal.*$") ~ "infaunal",
                             str_detect(tiering, "^.*epifaunal.*$") ~ "epifaunal",
                             TRUE ~ as.character(tiering))
  ) %>%
  group_by(species, motility, tiering, size) %>% 
  distinct() %>%
  write.csv(., "../data/clean/extinction/traits_trophic.csv",
            row.names = FALSE)

# create a minimum trophic community (collapse all species with the same traits)

trait_files <- list.files(path = "../data/clean/trait_minimum", pattern = ".csv", full.names = TRUE)

all_spp <- trait_files %>%
  lapply(read_csv) %>%
  bind_rows %>%
  select(-species) %>%
  group_by(feeding, motility, tiering, size) %>% 
  distinct() %>%
  ungroup() %>%
  mutate(species = paste0("sp_", row_number()))

for (i in seq_along(trait_files)) {
  
  df <- read.csv(trait_files[i]) %>%
    select(-species) %>%
    group_by(feeding, motility, tiering, size) %>% 
    distinct() %>% 
    ungroup() %>%
    left_join(all_spp)
  
  # write as clean data
  write.csv(df, str_replace(trait_files[i], "trait_minimum", "trophic_minimum"),
            row.names = FALSE)
}

# create a maximal trophic community (collapse all species with the same traits)

trait_files <- list.files(path = "../data/clean/trait_maximal", pattern = ".csv", full.names = TRUE)

all_spp <- trait_files %>%
  lapply(read_csv) %>%
  bind_rows %>%
  select(-species) %>%
  group_by(feeding, motility, tiering, size) %>% 
  distinct() %>%
  ungroup() %>%
  mutate(species = paste0("sp_", row_number()))

for (i in seq_along(trait_files)) {
  
  df <- read.csv(trait_files[i]) %>%
    select(-species) %>%
    group_by(feeding, motility, tiering, size) %>% 
    distinct() %>% 
    ungroup() %>%
    left_join(all_spp)
  
  # write as clean data
  write.csv(df, str_replace(trait_files[i], "trait_maximal", "trophic_maximal"),
            row.names = FALSE)
}
