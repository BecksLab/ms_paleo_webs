library(RColorBrewer)
library(here)
library(patchwork)
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
  select(-c(network, distance)) %>% 
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
                        .default = as.character(stat))) %>%
  mutate(level = case_when(
    stat %in% c("richness", "deficiency", "complexity", "connectance") ~ "Macro",
    stat %in% c("generality", "vulnerability") ~ "Micro",
    .default = "Meso"
  ))

df$id <- ordered(df$id, levels=c("pre", "during", "post"))

plot_list <- vector(mode = "list", length = 3)
levs = c("Macro", "Meso", "Micro")

for (i in seq_along(plot_list)) {
  
  plot_list[[i]] <- ggplot(df %>% 
              filter(level == levs[i]),
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
  scale_colour_brewer(palette = "Dark2") +
  labs(title = levs[i]) +
  theme(panel.border = element_rect(colour = 'black',
                                    fill = "#ffffff00"),
        axis.ticks.x = element_blank())
}

plot_list[[1]] / plot_list[[2]] / plot_list[[3]] +
  plot_layout(guides = 'collect') +
  plot_layout(height = c(2, 2, 1))

ggsave("../figures/summary_pfim.png",
       width = 4500,
       height = 7000,
       units = "px",
       dpi = 600)
