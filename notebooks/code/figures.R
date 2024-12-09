library(here)
library(readr)
library(scales)
library(tidyverse)

# set path to code sub dir
setwd(here("notebooks/code"))

#### Structure ####

df <- read_csv("../data/processed/nz_summary.csv") %>% 
  mutate(model = "real") %>% 
  mutate(across(matches("S[[:digit:]]"), log)) %>% 
  rbind(.,
        read_csv("../data/processed/topology_models.csv")) %>% 
  select(-richness) %>% 
  # to get the ratio
  mutate(ratio = top/basal,
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
                          .default = as.character(stat)))

ggplot() +
  geom_point(data = df %>% filter(model == "real"),
                aes(x = stat_val, 
                    y = id),
                shape = 4,
                size = 3) +
  geom_point(data = df %>% filter(model != "real"),
                aes(x = stat_val, 
                    y = id,
                    colour = model)) +
  facet_wrap(vars(stat),
             scales = 'free') +
  theme_classic() +
  xlab("value") +
  ylab("site") +
  coord_cartesian(clip = "off") +
  theme(panel.border = element_rect(colour = 'black',
                                    fill = "#ffffff00"))

ggsave("../figures/summary.png",
       width = 11000,
       height = 6000,
       units = "px",
       dpi = 600)

real_nets <- read_csv("../data/processed/nz_summary.csv") %>% 
  select(-richness) %>% 
  # to get the ratio
  mutate(ratio = top/basal,
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
  # to get the ratio
  mutate(ratio = top/basal,
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
                        .default = as.character(stat)))


ggplot(mod_nets) +
    geom_vline(aes(xintercept = 0)) +
    geom_histogram(aes(x = z_score,
                    fill = model),
                colour = "#ffffff00",
                alpha = 0.7) +
    facet_grid(rows = vars(model),
                cols = vars(stat),
                scales = "free") +
    scale_fill_discrete(guide = 'none') +
    coord_cartesian(expand = FALSE, clip = "off") +
    theme_classic() +
    theme(panel.border = element_rect(colour = 'black',
                                      fill = "#ffffff00"))

ggsave("../figures/zscore.png",
       width = 11000,
       height = 6000,
       units = "px",
       dpi = 600)
