library(RColorBrewer)
library(here)
library(patchwork)
library(readr)
library(scales)
library(stats)
library(tidyverse)

# set path to code sub dir
setwd(here("code"))

#### Structure ####

df <- read_csv("../data/processed/topology.csv") %>%
  mutate(across(matches("S[[:digit:]]"), log),
         top = top / richness,
         basal = basal / richness) %>%
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
    stat %in% c("richness", "complexity", "connectance", "links", "diameter", "distance", "redundancy") ~ "Macro",
    stat %in% c("generality", "vulnerability", "top", "basal") ~ "Micro",
    .default = "Meso"
  ))

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
  plot_layout(height = c(2, 1, 1))

ggsave("../figures/summary.png",
       width = 4500,
       height = 7000,
       units = "px",
       dpi = 600)

#### Extinctions ####

df_ext <- read_csv("../data/processed/extinction_topology.csv") %>%
  mutate(across(matches("S[[:digit:]]"), log),
         top = top / richness,
         basal = basal / richness) %>%
  # to get the ratio
  pivot_longer(
    cols = -c(model, extinction_mechanism, rep),
    names_to = "stat",
    values_to = "stat_val") %>%
  # standardise names
  mutate(stat = case_when(stat == "S1" ~ "No. of linear chains",
                          stat == "S2" ~ "No. of omnivory motifs",
                          stat == "S5" ~ "No. of apparent competition motifs",
                          stat == "S4" ~ "No. of direct competition motifs",
                          .default = as.character(stat))) %>% 
  mutate(level = case_when(
    stat %in% c("richness", "complexity", "connectance", "links", "diameter", "distance", "redundancy") ~ "Macro",
    stat %in% c("generality", "vulnerability", "top", "basal") ~ "Micro",
    .default = "Meso"
  ))

df_ext_summ <-
  df_ext %>%
  group_by(model, stat, level, extinction_mechanism) %>% 
  summarise(y = mean(stat_val, na.rm = TRUE))

ext_plot_list <- vector(mode = "list", length = 3)

for (i in seq_along(ext_plot_list)) {
  
  ext_plot_list[[i]] <- 
    ggplot(data = df_ext %>% 
             #filter(extinction_mechanism == "random") %>%
             filter(level == levs[i])) +
    geom_boxplot(aes(x = model,
                     y = stat_val, 
                     colour = model)) +
    facet_wrap(vars(stat),
               scales = 'free',
               ncol = 2) +
    scale_size(guide = 'none') +
    theme_classic() +
    xlab("extinction mechanism") +
    ylab("value") +
    coord_cartesian(clip = "off") +
    scale_colour_brewer(palette = "Dark2") +
    labs(title = levs[i]) +
    theme(panel.border = element_rect(colour = 'black',
                                      fill = "#ffffff00"),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank())
}

ext_plot_list[[1]] / ext_plot_list[[2]] / ext_plot_list[[3]] +
  plot_layout(guides = 'collect') +
  plot_layout(height = c(2, 2, 1))

ggsave("../figures/extinction.png",
       width = 4500,
       height = 9000,
       units = "px",
       dpi = 600)

# groups plot
read_csv("../data/processed/extinction_topology.csv") %>%
  group_by(model, extinction_mechanism) %>% 
  count() %>%
  ggplot(., aes(x=extinction_mechanism, y=n, fill=model))  +
  scale_fill_brewer(palette = "Dark2")+ 
  geom_bar(position="dodge", stat="identity")


ext_lins =
  df %>%
  filter(time %in% c("G1", "G2")) %>%
  mutate(known = "known") %>%
  rbind(df_ext %>%
          mutate(time = "G2",
                 known = "simulated",
                 rep = NULL,
                 extinction_mechanism = NULL)) %>%
  rbind(df %>%
          filter(time == "G1") %>%
          mutate(known = "simulated")) %>%
  group_by(model, time, stat, level, known) %>%
  summarise(stat_val = mean(stat_val, na.rm = TRUE)) %>%
  filter(stat != "robustness") %>%
  filter(stat != "r50")

for (i in seq_along(ext_plot_list)) {
  
  ext_plot_list[[i]] <- ggplot(ext_lins %>% 
                                 filter(level == levs[i])) +
    geom_line(aes(x = time,
                  y = stat_val,
                  colour = model,
                  group = paste(model, known),
                  linetype = known)) +
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

ext_plot_list[[1]] / ext_plot_list[[2]] / ext_plot_list[[3]] +
  plot_layout(guides = 'collect') +
  plot_layout(height = c(2, 1, 1))

ggsave("../figures/extinction_all_results.png",
       width = 4500,
       height = 7000,
       units = "px",
       dpi = 600)

#### PCA ####

df_pca <- list.files(path = "../data/processed/topology/", pattern = ".csv", full.names = TRUE) %>% 
  lapply(read_csv) %>% 
  bind_rows %>% 
  mutate(across(matches("S[[:digit:]]"), log)) %>% 
  select(-c(network, distance)) %>% 
  filter(model %in% c("adbm", "bodymassratio", "niche", "pfim_minimum", "random", "lmatrix")) %>%
  # remove
  mutate(top = NULL,
         basal = NULL,
         id = case_when(
           str_detect(id, "^.*pre.*$") ~ "pre",
           str_detect(id, "^.*during.*$") ~ "post",
           str_detect(id, "^.*post.*$") ~ "during",
           TRUE ~ as.character(id))) %>%
  drop_na()

pca_res <- prcomp(df_pca[3:8], scale. = TRUE)

autoplot(pca_res, data = df_pca, colour = 'model', shape = 'id', size = 4, alpha = 0.7) +
  theme_classic() +
  #scale_colour_brewer(palette = "Dark2") +
  theme(panel.border = element_rect(colour = 'black',
                                    fill = "#ffffff00"),
        axis.ticks.x = element_blank())

ggsave("../figures/pca.png",
       width = 4500,
       height = 4000,
       units = "px",
       dpi = 600)
