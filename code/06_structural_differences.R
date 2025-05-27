
# libraries
library(ggpubr)
library(here)
library(patchwork)
library(rstatix)
library(tidyverse)

# set path to code sub dir
setwd(here("code"))

# import data

df <- read_csv("../data/processed/topology.csv") %>%
  #mutate(across(matches("S[[:digit:]]"), log)) %>%
  select(-c(richness, redundancy)) %>%
  # to get the ratio
  pivot_longer(
    cols = -c(model, time),
    names_to = "stat",
    values_to = "stat_val") %>%
  # standardise names
  mutate(stat = case_when(stat == "S1" ~ "No. of linear chains",
                          stat == "S2" ~ "No. of omnivory motifs",
                          stat == "S5" ~ "No. of apparent competition motifs",
                          stat == "S4" ~ "No. of direct competition motifs",
                          .default = as.character(stat))) %>%
  mutate(level = case_when(
    stat %in% c("complexity", "connectance", "diameter", "distance") ~ "Macro",
    stat %in% c("generality", "vulnerability") ~ "Micro",
    .default = "Meso"
  )) %>%
  mutate(model = factor(model, ordered = TRUE, 
                        levels = c("niche", "random", "adbm", "lmatrix", "pfim", "bodymassratio")))

lm1 <- lm(stat_val ~ model + time, df %>% filter(stat == "connectance"))
summary(lm1)

df %>% 
  filter(stat == "generality") %>%
  filter(time == "G1") %>%
  wilcox_test(stat_val ~ model,
              paired = FALSE)

df %>% 
  filter(stat == "generality") %>%
  filter(time == "G1") %>%
  t_test(stat_val ~ model,
         var.equal = TRUE,
         detailed = TRUE) %>% view()

plot_list <- vector(mode = "list", length = 3)
levs = c("Macro", "Meso", "Micro")

for (i in seq_along(plot_list)) {
  
  plot_list[[i]] <- ggplot(df %>% 
                             filter(level == levs[i]),
                           aes(x = time,
                               y = stat_val,
                               colour = model)) +
    geom_boxplot(position=position_dodge(1)) +
    facet_wrap(vars(stat),
               scales = 'free',
               ncol = 2) +
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

ggsave("../figures/summary.png",
       width = 5000,
       height = 7000,
       units = "px",
       dpi = 600)
