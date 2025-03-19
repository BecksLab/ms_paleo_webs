library(here)
library(tidyverse)

# set path to code sub dir
setwd(here("code"))

# create simple traits df for extinction sims

file_list <- list.files(path = "../data/raw", pattern = "species.csv", full.names = TRUE)

for (i in seq_along(file_list)) {
  
  read.csv(file_list[i])%>%
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
  distinct()  %>%
    write.csv(., str_replace(file_list[i], "raw", "extinction"),
              row.names = FALSE)

}
