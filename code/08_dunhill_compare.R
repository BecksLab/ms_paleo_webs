# libraries
library(here)
library(patchwork)
library(tidyverse)

# set path to code sub dir
setwd(here("code"))

#load script that determines plotting aesthetics
source("lib/plotting_theme.R")

# import simulated data

df <- read_csv("../data/processed/topology.csv") %>%
  # remove metaweb pfims
  filter(model != "pfim_metaweb") %>%
  # rename the remianing pfim col
  mutate(model = case_when(model == "pfim_downsample" ~ "pfim",
                           .default = as.character(model)))  %>%
  mutate(model = str_replace(model, "bodymassratio", "log ratio")) %>%
  #mutate(across(matches("S[[:digit:]]"), log)) %>%
  select(-c(distance, redundancy, complexity, diameter, n_rep)) %>%
  pivot_longer(
    cols = -c(model, time),
    names_to = "stat",
    values_to = "stat_val") %>%
  # get mean values
  group_by(model, time, stat) %>%
  summarise(mean = mean(stat_val)) %>%
  ungroup() %>%
  # standardise names
  mutate(stat = case_when(stat == "S1" ~ "No. of linear chains",
                          stat == "S2" ~ "No. of omnivory motifs",
                          stat == "S5" ~ "No. of apparent competition motifs",
                          stat == "S4" ~ "No. of direct competition motifs",
                          .default = as.character(stat))) %>%
  mutate(level = case_when(stat %in% c("richness", "connectance", "generality", "vulnerability", "trophic_level") ~ "struct",
                           .default = "motif"),
         model = factor(model, ordered = TRUE, 
                        levels = c("niche", "random", "adbm", "lmatrix", "log ratio", "pfim")),
         time = str_extract(time, "\\d+"))

struct <- ggplot(df %>% 
                   filter(level == "struct"),
                 aes(x = time,
                     y = mean,
                     colour = model,
                     group = model)) +
  geom_line() +
  facet_wrap(vars(stat),
             scales = 'free_y',
             ncol = 1) +
  scale_size(guide = 'none') +
  xlab(NULL) +
  ylab(NULL)  +
  coord_cartesian(clip = "off")  +
  scale_colour_manual(values = pal_df$c,
                      breaks = pal_df$l) +
  scale_x_discrete(labels = c("pre", "during", "early", "late")) +
  figure_theme +
  theme(legend.position = 'none')

motif <- ggplot(df %>% 
                  filter(level == "motif"),
                aes(x = time,
                    y = mean,
                    colour = model,
                    group = model)) +
  geom_line() +
  facet_wrap(vars(stat),
             scales = 'free_y',
             ncol = 1) +
  scale_size(guide = 'none') +
  xlab(NULL) +
  ylab(NULL) +
  coord_cartesian(clip = "off") +
  scale_colour_manual(values = pal_df$c,
                      breaks = pal_df$l) +
  scale_x_discrete(labels = c("pre", "during", "early", "late")) +
  figure_theme +
  theme(legend.position = 'none')

# means abs differences between real and extinction simulations

df <- read_csv("../data/processed/extinction_topology.csv") %>%
  # remove metaweb pfims
  filter(model != "pfim_metaweb") %>%
  # rename the remianing pfim col
  mutate(model = case_when(model == "pfim_downsample" ~ "pfim",
                           .default = as.character(model))) %>%
  select(-c(rep, distance, redundancy, complexity, diameter, S1, S2, S4, S5, resilience)) %>%
  pivot_longer(
    cols = -c(model, extinction_mechanism, n_rep),
    names_to = "stat",
    values_to = "sim_val") %>%
  full_join(.,
            read_csv("../data/processed/topology.csv") %>%
              # remove metaweb pfims
              filter(model != "pfim_metaweb") %>%
              # rename the remianing pfim col
              mutate(model = case_when(model == "pfim_downsample" ~ "pfim",
                                       .default = as.character(model))) %>%
              filter(time == "G2") %>%
              select(-c(distance, redundancy, complexity, diameter, time, S1, S2, S4, S5, trophic_level)) %>%
              pivot_longer(
                cols = -c(model, n_rep),
                names_to = "stat",
                values_to = "real_val")) %>%
  mutate(model = str_replace(model, "bodymassratio", "log ratio")) %>%
  mutate(diff = real_val - sim_val) %>%
  group_by(model, extinction_mechanism) %>%
  summarise(mean_diff = abs(mean(diff, na.rm = TRUE))) %>%
  mutate(model = factor(model, ordered = TRUE, 
                        levels = c("niche", "random", "adbm", "lmatrix", "log ratio", "pfim")))

mean_diff <- ggplot(df,
                    aes(y = extinction_mechanism,
                        x = mean_diff,
                        fill = model)) +
  geom_bar(stat="identity", position=position_dodge()) +
  coord_cartesian(clip = "off") +
  ylab(NULL)  +
  scale_fill_manual(values = pal_df$c,
                      breaks = pal_df$l) +
  figure_theme

# tss scores

df <- read_csv("../data/processed/extinction_tss.csv") %>%
  # remove metaweb pfims
  filter(model != "pfim_metaweb") %>%
  # rename the remianing pfim col
  mutate(model = case_when(model == "pfim_downsample" ~ "pfim",
                           .default = as.character(model)))%>%
  mutate(model = str_replace(model, "bodymassratio", "log ratio")) %>%
  group_by(model, extinction_mechanism) %>%
  summarise(mean_tss = mean(tss, na.rm = TRUE))

tss <-  ggplot(df,
               aes(y = extinction_mechanism,
                   x = mean_tss,
                   fill = model)) +
  geom_bar(stat="identity", position=position_dodge()) +
  coord_cartesian(clip = "off") +
  ylab(NULL)  +
  scale_fill_manual(values = pal_df$c,
                      breaks = pal_df$l) +
  figure_theme


(tss + struct) / (mean_diff + motif)  +
  plot_layout(guides = 'collect') +
  plot_annotation(tag_levels = 'A')

ggsave("../figures/dunhill_comp.png",
       width = 7000,
       height = 7500,
       units = "px",
       dpi = 600)

