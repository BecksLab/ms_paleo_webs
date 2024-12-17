library(here)
library(patchwork)
library(RColorBrewer)
library(readr)
library(scales)
library(tidyverse)

# set path to code sub dir
setwd(here("notebooks/code"))

#### Structure ####

df <- read_csv("../data/processed/nz_summary.csv") %>% 
  mutate(model = "real") %>% 
  rbind(.,
        read_csv("../data/processed/topology_models.csv")) %>% 
  select(-richness) %>% 
  mutate(across(matches("S[[:digit:]]"), log)) %>% 
  mutate(distance = NULL,
         top = NULL,
         basal = NULL) %>%
  pivot_longer(
    cols = -c(id, model), 
    names_to = "stat",
    values_to = "stat_val") %>% 
  group_by(id, model, stat) %>% 
  mutate(stat_val = mean(stat_val, na.rm = TRUE),
         sd = sd(stat_val, na.rm = TRUE),
         stat = case_when(stat == "S1" ~ "No. of linear chains",
                          stat == "S2" ~ "No. of omnivory motifs",
                          stat == "S4" ~ "No. of apparent competition motifs",
                          stat == "S5" ~ "No. of direct competition motifs",
                          .default = as.character(stat))) %>%
  mutate(level = case_when(
    stat %in% c("richness", "deficiency", "complexity", "connectance") ~ "Macro",
    stat %in% c("generality", "vulnerability") ~ "Micro",
    .default = "Meso"))

plot_list <- vector(mode = "list", length = 3)
levs = c("Macro", "Meso", "Micro")

for (i in seq_along(plot_list)) {
  
  plot_list[[i]] <- ggplot() +
  geom_point(data = df %>% filter(model == "real") %>% filter(level == levs[i]),
                aes(x = stat_val, 
                    y = id),
                shape = 4,
                size = 3) +
  geom_point(data = df %>% filter(model != "real") %>% filter(level == levs[i]),
                aes(x = stat_val, 
                    y = id,
                    colour = model)) +
  facet_wrap(vars(stat),
             scales = 'free',
             ncol = 2) +
  theme_classic() +
  labs(title = levs[i]) +
  xlab("value") +
  ylab("site") +
  scale_colour_brewer(palette = "Dark2") +
  coord_cartesian(clip = "off") +
  theme(panel.border = element_rect(colour = 'black',
                                    fill = "#ffffff00"))
}

plot_list[[1]] / plot_list[[2]] / plot_list[[3]] +
  plot_layout(guides = 'collect') +
  plot_layout(height = c(2, 2, 1))

ggsave("../figures/summary_contemporary.png",
       width = 5000,
       height = 8000,
       units = "px",
       dpi = 600)

real_nets <- read_csv("../data/processed/nz_summary.csv") %>% 
  select(-richness) %>% 
  mutate(across(matches("S[[:digit:]]"), log)) %>% 
  # to get the ratio
  mutate(distance = NULL,
         top = NULL,
         basal = NULL) %>%
  pivot_longer(
    cols = -c(id), 
    names_to = "stat",
    values_to = "stat_val") %>% 
    group_by(stat, id) %>% 
    reframe(real_mu = mean(stat_val, na.rm = TRUE))

mod_nets <- read_csv("../data/processed/topology_models.csv") %>% 
  select(-richness) %>% 
  mutate(across(matches("S[[:digit:]]"), log)) %>% 
  # to get the ratio
  mutate(distance = NULL,
         top = NULL,
         basal = NULL) %>%
  pivot_longer(
    cols = -c(id, model), 
    names_to = "stat",
    values_to = "stat_val") %>% 
    group_by(model, stat, id) %>% 
    reframe(model_mu = mean(stat_val, na.rm = TRUE),
            model_sd = sd(stat_val, na.rm = TRUE)) %>% 
   left_join(., real_nets) %>% 
   mutate(z_score = (real_mu - model_mu)/model_sd,
       stat = case_when(stat == "S1" ~ "No. of linear chains",
                        stat == "S2" ~ "No. of omnivory motifs",
                        stat == "S4" ~ "No. of apparent competition motifs",
                        stat == "S5" ~ "No. of direct competition motifs",
                        .default = as.character(stat))) %>%
   mutate(level = case_when(
    stat %in% c("richness", "deficiency", "complexity", "connectance") ~ "Macro",
    stat %in% c("generality", "vulnerability") ~ "Micro",
    .default = "Meso"))

for (i in seq_along(plot_list)) {

  plot_list[[i]] <-  ggplot(mod_nets %>% filter(level == levs[i])) +
    geom_vline(aes(xintercept = 0)) +
    geom_histogram(aes(x = z_score,
                    fill = model),
                colour = "#ffffff00",
                alpha = 0.7) +
    facet_grid(rows = vars(model),
                cols = vars(stat),
                scales = "free") +
    scale_fill_brewer(palette = "Dark2") +
    coord_cartesian(expand = FALSE, clip = "off") +
    theme_classic() +
    labs(title = levs[i]) +
    theme(panel.border = element_rect(colour = 'black',
                                      fill = "#ffffff00"),
          legend.position = 'none')
  
}

plot_list[[1]] / plot_list[[2]] / plot_list[[3]] +
  plot_layout(guides = 'collect') +
  plot_layout(height = c(1, 1, 1))

ggsave("../figures/zscore_contemporary.png",
       width = 11000,
       height = 17000,
       units = "px",
       dpi = 600)
