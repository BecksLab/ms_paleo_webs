library(tidyverse)

# clean traits data
trait_files <- list.files(path = "data/raw", pattern = "traits.csv", full.names = TRUE)

for (i in seq_along(trait_files)) {

    df <- read.csv(trait_files[i])  %>%
        # standardise strings
        mutate_all(~str_replace_all(., "-", "_"))  %>%
        # remove for now
        filter(!(feeding %in% c("primary", "zooplankton")))  %>%
        # change motiliy class - we can remove this later
        mutate(motility = case_when( 
               motility == "nonmotile" ~ "non_motile",
               TRUE ~ as.character(motility)))  %>%
        # divorce feeding from size
        mutate(size = case_when(
            str_detect(size, "^.*small.*$") ~ "small",
            str_detect(size, "^.*large.*$") ~ "large",
            str_detect(size, "^.*medium.*$") ~ "medium",
            str_detect(size, "^.*tiny.*$") ~ "tiny",
            species == "Hemigordius baoqingensis" ~ "medium",
            TRUE ~ as.character(size)),
            # same for feeding
            feeding = case_when(
                        feeding == "microcarnivore" ~ "carnivore",
                        TRUE ~ as.character(feeding)))

    # write as clean data
    write.csv(df, str_replace(trait_files[i], "raw.", "clean/trait/"),
            row.names = FALSE)
}

# clean size data

size_files <- list.files(path = "data/raw", pattern = "size.csv", full.names = TRUE)

for (i in seq_along(size_files)) {
   
   df <- read.csv(size_files[i])  %>%
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