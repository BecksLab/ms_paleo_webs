# libraries
library(here)
library(patchwork)
library(tidyverse)

# set path to code sub dir
setwd(here("code"))

#load script that determines plotting aesthetics
source("lib/plotting_theme.R")

# import simulated data

df <- read_csv("../data/processed/robustness.csv") %>%
  # standardise names
  glow_up(time = case_when(time == "G1" ~ "pre",
                          time == "G2" ~ "during",
                          time == "G3" ~ "early",
                          time == "G4" ~ "late",
                          .default = as.character(time))) %>%
  glow_up(time = factor(time, ordered = TRUE, 
                       levels = c("pre", "during", "early", "late"))) %>%
  squad_up(time, threshold) %>%
  no_cap(robustness = mean(robustness))

ggplot(df,
       aes(x = threshold,
           y = robustness*100,
           colour = time,
           group = time)) +
  geom_abline(slope = -1, 
              intercept = 100, 
              linetype = 'dotted') +
  geom_line() +
  geom_smooth(linewidth = 0.5) +
  ylim(0,100)+
  labs(x = "Percent Primary Extinction",
       y = "Robustness (% Community Remaining; n = 500)") +
  guides(color=guide_legend("Network")) +
  figure_theme

ggsave("../figures/robustness.png",
       width = 5000,
       height = 4500,
       units = "px",
       dpi = 600)

