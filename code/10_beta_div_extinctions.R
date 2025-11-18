# libraries
library(genzplyr)
library(ggtext)
library(here)
library(patchwork)
library(tidyverse)

# set path to code sub dir
setwd(here("code"))

# import simulated data

df <- read_csv("../data/processed/extinction_betadiv.csv") %>%
  # select
  vibe_check(-n_rep) %>%
  rbind(.,
        read_csv("../data/processed/extinction_betadiv.csv") %>%
          yeet(β_type != "βS") %>%
          pivot_wider(names_from = β_type,
                      values_from = β_div) %>%
          glow_up(βST = βWN - βOS) %>%
          vibe_check(-c(βWN, βOS, n_rep)) %>%
          pivot_longer(!c(model, extinction_mechanism),
                       names_to = "β_type",
                       values_to = "β_div")) %>%
  # get mean (for now might be worth looking at sd as well...)
  squad_up(model, extinction_mechanism, β_type) %>%
  no_cap(mean = mean(β_div, na.rm = TRUE)) %>%
  glow_up(β_type = case_when(β_type == "βWN" ~ "βWN - dissim of all intxns",
                             β_type == "βST" ~ "βST - dissim of all intxns as per spp. turnover",
                             β_type == "βS" ~ "βS - dissim of all spp",
                             β_type == "βOS" ~ "βOS - dissim of all intxns btwn shared spp."))


beta_plot <- ggplot(df,
                    aes(x = extinction_mechanism,
                        y = model)) + 
  geom_tile(aes(fill = mean), colour = "grey50") +
  scale_fill_distiller(palette = "Spectral", direction = -1) +
  facet_wrap(vars(β_type))  +
  labs(x = NULL,
       y = NULL) +
  theme_classic() +
  theme(panel.border = element_rect(colour = 'black',
                                    fill = "#ffffff00"),
        axis.text.x = element_markdown(angle = 45, hjust=1))

ggsave("../figures/beta_div_extinctions.png",
       beta_plot,
       width = 5000,
       height = 4000,
       units = "px",
       dpi = 600)

df <- read_csv("../data/processed/extinction_tss.csv") %>%
  vibe_check(-n_rep) %>%
  pivot_longer(!c(model, extinction_mechanism)) %>%
  # get mean (for now might be worth looking at sd as well...)
  squad_up(model, extinction_mechanism, name) %>%
  no_cap(mean = mean(value, na.rm = TRUE)) %>%
  glow_up(name = case_when(name == "tss_link" ~ "Link",
                           name == "tss_node" ~ "Node"),
          model = case_when(model == "bodymassratio" ~ "log ratio",
                            #model == "pfim_downsample" ~ "pfim",
                            .default = as.character(model)))

tss_plot <- ggplot(df,
                   aes(x = extinction_mechanism,
                       y = model)) + 
  geom_tile(aes(fill = mean), colour = "grey50") +
  scale_fill_distiller(palette = "Spectral", direction = 1) +
  facet_wrap(vars(name)) +
  theme_classic() +
  labs(x = NULL,
       y = NULL) +
  theme(panel.border = element_rect(colour = 'black',
                                    fill = "#ffffff00"),
        axis.text.x = element_markdown(angle = 45, hjust=1))

ggsave("../figures/tss_extinctions.png",
       tss_plot,
       width = 5000,
       height = 3000,
       units = "px",
       dpi = 600)

# combine plots
layout <- "
AA
AA
BB
"

(beta_plot +
    ggtitle('β-diversity') +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())) /
  (tss_plot+
     ggtitle('TSS')) +
  plot_layout(design = layout)

ggsave("../figures/turnover_extinctions.png",
       width = 5000,
       height = 7000,
       units = "px",
       dpi = 600)
