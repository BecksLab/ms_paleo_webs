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
    cols = -c(location, model, time),
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
                             filter(level == levs[i]) %>%
                             filter(time <= 3),
                           aes(x = time,
                               y = stat_val,
                               colour = model)) +
    geom_smooth(alpha = 0.3,
                weight = 0.5, 
                se = FALSE) +
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

df_ext <- read_csv("../data/processed/extinctions/extinctions.csv") %>% 
  bind_rows %>% 
  mutate(across(matches("S[[:digit:]]"), log)) %>% 
  select(-c(id, links, richness)) %>% 
  # to get the ratio
  mutate(top = NULL,
         basal = NULL,
         id = time,
         time = NULL,
         distance = NULL) %>%
  pivot_longer(
    cols = -c(id, model, extinction_mechanism), 
    names_to = "stat",
    values_to = "end_val") %>% 
  filter(id == "pre") %>% 
  mutate(id = "post",
         stat = case_when(stat == "S1" ~ "No. of linear chains",
                          stat == "S2" ~ "No. of omnivory motifs",
                          stat == "S4" ~ "No. of apparent competition motifs",
                          stat == "S5" ~ "No. of direct competition motifs",
                          .default = as.character(stat))) %>% 
  left_join(.,
            df %>%
              filter(id == "pre") %>% 
              select(-id)) %>% 
  mutate(xstart = "pre",
         xend = id,
         start_val = stat_val,
         stat_val = NULL,
         id = NULL) %>%
  mutate(level = case_when(
    stat %in% c("richness", "deficiency", "complexity", "connectance") ~ "Macro",
    stat %in% c("generality", "vulnerability") ~ "Micro",
    .default = "Meso"
  )) %>%
  filter(model %in% c("adbm", "bodymassratio", "niche", "pfim_minimum", "random", "lmatrix"))

df_ext$xstart <- ordered(df_ext$xstart, levels = c("pre", "during", "post"))
df_ext$xend <- ordered(df_ext$xend, levels = c("pre", "during", "post"))

df_ext_summ <-
  df_ext %>%
  select(model, stat, end_val, start_val, level) %>% 
  pivot_longer(cols = c(end_val, start_val)) %>% 
  group_by(model, stat, name, level) %>% 
  summarise(y = mean(value, na.rm = TRUE)) %>%
  mutate(name = ifelse(name == "end_val",
                       "post",
                       "pre"))

ext_plot_list <- vector(mode = "list", length = 3)

for (i in seq_along(ext_plot_list)) {
  
  ext_plot_list[[i]] <- plot_list[[i]] +
    geom_line(data = df_ext_summ %>% 
                filter(level == levs[i]),
              aes(x = factor(name),
                  y = y, 
                  group = model),
              linetype = "dashed")
}

ext_plot_list[[1]] / ext_plot_list[[2]] / ext_plot_list[[3]] +
  plot_layout(guides = 'collect') +
  plot_layout(height = c(2, 2, 1))

ggsave("../figures/extinction.png",
       width = 4500,
       height = 7000,
       units = "px",
       dpi = 600)

for (i in seq_along(ext_plot_list)) {
  
  ext_plot_list[[i]] <- ggplot() +
    geom_segment(data = df_ext %>% 
                   filter(level == levs[i]),
                 aes(x = xstart,
                     y = start_val, 
                     xend = xend,
                     yend = end_val,
                     colour = model,
                     group = model,
                     linetype = extinction_mechanism),
                 alpha = 0.3) +
    geom_line(data = df %>% filter(id != "during")%>% 
                filter(level == levs[i]),
              aes(x = factor(id),
                  y = stat_val,
                  colour = model,
                  group = model)) +
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
  plot_layout(height = c(2, 2, 1))

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
