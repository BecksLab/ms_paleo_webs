library(here)
library(readr)
library(scales)
library(tidyverse)

# set path to code sub dir
setwd(here("notebooks/code"))

#### Structure ####

df <- list.files(path = "../../data/processed/", pattern = ".csv", full.names = TRUE) %>% 
  lapply(read_csv) %>% 
  bind_rows %>% 
  mutate(across(matches("S[[:digit:]]"), log)) %>% 
  select(-c(network, richness, distance)) %>% 
  filter(str_detect(model, "^.*pfim.*$")) %>%
  # remove
  mutate(top = NULL,
         basal = NULL) %>%
  pivot_longer(
    cols = -c(id, model), 
    names_to = "stat",
    values_to = "stat_val")  %>% 
  # standardise names
  mutate(id = case_when(
    str_detect(id, "^.*pre.*$") ~ "pre",
    str_detect(id, "^.*during.*$") ~ "post",
    str_detect(id, "^.*post.*$") ~ "during",
    TRUE ~ as.character(id)),
       stat = case_when(stat == "S1" ~ "No. of linear chains",
                        stat == "S2" ~ "No. of omnivory motifs",
                        stat == "S4" ~ "No. of apparent competition motifs",
                        stat == "S5" ~ "No. of direct competition motifs",
                        .default = as.character(stat)))

df$id <- ordered(df$id, levels=c("pre", "during", "post"))

ggplot(df,
       aes(x = factor(`id`), 
           y = stat_val, 
           colour = model,
           group = model)) +
  geom_line() +
  geom_point() +
  facet_wrap(vars(stat),
             scales = 'free') +
  scale_size(guide = 'none') +
  theme_classic() +
  xlab("time") +
  ylab("value") +
  coord_cartesian(clip = "off") +
  theme(panel.border = element_rect(colour = 'black',
                                    fill = "#ffffff00"),
        axis.ticks.x = element_blank())

ggsave("../figures/summary_pfim.png",
       width = 11000,
       height = 6000,
       units = "px",
       dpi = 600)
