library(ggfortify)
library(here)
library(patchwork)
library(RColorBrewer)
library(readr)
library(scales)
library(stats)
library(tidyverse)
library(vegan)

# set path to code sub dir
setwd(here("notebooks/code"))

#### Structure ####

df <- list.files(path = "../../data/processed/topology/", pattern = ".csv", full.names = TRUE) %>% 
  lapply(read_csv) %>% 
  bind_rows %>% 
  mutate(across(matches("S[[:digit:]]"), log)) %>% 
  select(-c(network, deficiency)) %>% 
  filter(str_detect(model, "^.*pfim.*$")) %>%
  # remove
  mutate(top = NULL,
         basal = NULL) %>%
  pivot_longer(
    cols = -c(id, model), 
    names_to = "stat",
    values_to = "stat_val")  %>% 
  # standardise names
  mutate(id = case_when(str_detect(id, "^.*pre.*$") ~ "pre",
                        str_detect(id, "^.*during.*$") ~ "during",
                        str_detect(id, "^.*post.*$") ~ "post",
                        TRUE ~ as.character(id)),
         stat = case_when(stat == "S1" ~ "No. of linear chains",
                          stat == "S2" ~ "No. of omnivory motifs",
                          stat == "S4" ~ "No. of apparent competition motifs",
                          stat == "S5" ~ "No. of direct competition motifs",
                          stat == "distance" ~ "chain length",
                          .default = as.character(stat)),
         model_broad = case_when(str_detect(model, "maximal") ~ "maximal",
                                 TRUE ~ "minmum"),
         node = case_when(str_detect(model, "trophic") ~ "trophic",
                          TRUE ~ "taxonomic"),
         downsample = case_when(str_detect(model, "downsample") ~ "downsample",
                                TRUE ~ "metaweb"),
         model_simple = str_replace_all(model, "_(minimum|maximal)", "")) %>%
  mutate(level = case_when(
    stat %in% c("richness", "chain length", "complexity", "connectance") ~ "Macro",
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
                               colour = model_simple,
                               linetype = model_broad,
                               group = model)) +
    geom_line(alpha = 0.7) +
    geom_point(alpha = 0.7) +
    facet_wrap(vars(stat),
               scales = 'free') +
    scale_size(guide = 'none') +
    theme_classic() +
    xlab("time") +
    ylab("value") +
    coord_cartesian(clip = "off") +
    #scale_colour_brewer(palette = "Dark2") +
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

#### Redundancy ####

df_redun <- list.files(path = "../../data/processed/topology/", pattern = ".csv", full.names = TRUE) %>% 
  lapply(read_csv) %>% 
  bind_rows %>% 
  select(-network) %>%
  filter(str_detect(model, "^.*pfim.*$")) %>%
  mutate(links = round((connectance * richness^2), 0),
         redundancy = links - (richness - 1),
         id = case_when(str_detect(id, "^.*pre.*$") ~ "pre",
                        str_detect(id, "^.*during.*$") ~ "during",
                        str_detect(id, "^.*post.*$") ~ "post",
                        TRUE ~ as.character(id)),
         model_broad = case_when(str_detect(model, "maximal") ~ "maximal",
                                 TRUE ~ "minmum"),
         species = case_when(str_detect(model, "trophic") ~ "trophic",
                             TRUE ~ "taxonomic"),
         downsample = case_when(str_detect(model, "downsample") ~ "downsample",
                                TRUE ~ "metaweb"),
         model_simple = str_replace_all(model, "_(minimum|maximal|trophic)", ""))

df_redun$id <- ordered(df_redun$id, levels=c("pre", "during", "post"))

ggplot(df_redun,
       aes(x = factor(`id`), 
           y = log(redundancy), 
           colour = model_simple,
           group = model,
           linetype = species)) +
  geom_line(alpha = 0.7) +
  geom_point(alpha = 0.7) +
  facet_wrap(vars(model_broad)) +
  scale_size(guide = 'none') +
  theme_classic() +
  xlab("time") +
  ylab("log(redundancy)") +
  coord_cartesian(clip = "off") +
  theme(panel.border = element_rect(colour = 'black',
                                    fill = "#ffffff00"),
        axis.ticks.x = element_blank())

ggsave("../figures/redundancy_pfim.png",
       width = 4500,
       height = 3500,
       units = "px",
       dpi = 600)

#### PCA ####

df_pca <- list.files(path = "../../data/processed/topology/", pattern = ".csv", full.names = TRUE) %>% 
  lapply(read_csv) %>% 
  bind_rows %>% 
  select(-c(network)) %>% 
  filter(str_detect(model, "^.*pfim.*$")) %>%
  # remove
  mutate(top = NULL,
         basal = NULL,
         links = round((connectance * richness^2), 0),
         redundancy = links - (richness - 1),,
         id = case_when(str_detect(id, "^.*pre.*$") ~ 1,
                        str_detect(id, "^.*during.*$") ~ 2,
                        str_detect(id, "^.*post.*$") ~ 3),
         model_broad = case_when(str_detect(model, "maximal") ~ "maximal",
                                 TRUE ~ "minmum"),
         node = case_when(str_detect(model, "trophic") ~ "trophic",
                          TRUE ~ "taxonomic"),
         downsample = case_when(str_detect(model, "downsample") ~ "downsample",
                                TRUE ~ "metaweb"),
         model_simple = str_replace_all(model, "_(minimum|maximal|trophic)", "")) %>%
  drop_na()

ord <- metaMDS(log1p(df_pca[3:15]))
fit <- envfit(ord, df_pca[c(1,16:18)], perm = 9999)

df_pca <-
  cbind(df_pca, ord[["points"]])

ggplot() +
  geom_point(data = df_pca,
             aes(x = MDS1,
                 y = MDS2,
                 colour = model_simple,
                 shape = node,
                 size = as.factor(id)),
             alpha = 0.5) +
  geom_segment(aes(x = 0, y = 0, 
                   xend = c(fit[["vectors"]][["arrows"]][1]), 
                   yend = fit[["vectors"]][["arrows"]][2]),
               colour = "black",
               arrow = arrow(length = unit(0.03, "npc"))) +
  geom_text(aes(label = "time", 
                x = fit[["vectors"]][["arrows"]][1], 
                y = fit[["vectors"]][["arrows"]][2]), 
            nudge_x = 0.15,
            colour = "black") +
  theme_classic()

ggsave("../figures/pca_pfim.png",
       width = 3500,
       height = 2500,
       units = "px",
       dpi = 600)

#### Extinctions ####

df_ext <- read_csv("../../data/processed/extinctions/extinctions.csv") %>% 
  filter(str_detect(model, "^.*pfim.*$")) %>% 
  mutate(across(matches("S[[:digit:]]"), log)) %>% 
  mutate(top = NULL,
         basal = NULL,
         id = time,
         time = NULL) %>%
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
                          stat == "distance" ~ "chain length",
                          .default = as.character(stat))) %>% 
  inner_join(.,
             list.files(path = "../../data/processed/topology/", pattern = ".csv", full.names = TRUE) %>% 
               lapply(read_csv) %>% 
               bind_rows %>% 
               select(-network) %>%
               filter(str_detect(model, "^.*pfim.*$")) %>%
               mutate(links = round((connectance * richness^2), 0),
                      redundancy = links - (richness - 1),
                      id = case_when(str_detect(id, "^.*pre.*$") ~ 1,
                                     str_detect(id, "^.*during.*$") ~ 2,
                                     str_detect(id, "^.*post.*$") ~ 3))  %>%
               filter(id == 1) %>%
               filter(redundancy == max(redundancy) | redundancy == min(redundancy)) %>% 
               select(-id) %>%
               pivot_longer(
                 cols = -c(model), 
                 names_to = "stat",
                 values_to = "stat_val"),
             mutate(stat = case_when(stat == "S1" ~ "No. of linear chains",
                                     stat == "S2" ~ "No. of omnivory motifs",
                                     stat == "S4" ~ "No. of apparent competition motifs",
                                     stat == "S5" ~ "No. of direct competition motifs",
                                     stat == "distance" ~ "chain length",
                                     .default = as.character(stat))))%>%
  mutate(level = case_when(
    stat %in% c("richness", "chain length", "complexity", "connectance") ~ "Macro",
    stat %in% c("generality", "vulnerability") ~ "Micro",
    .default = "Meso"
  )) %>% 
  mutate(xstart = "pre",
         xend = id,
         start_val = stat_val,
         stat_val = NULL,
         id = NULL)

df_ext$xstart <- ordered(df_ext$xstart, levels = c("pre", "during", "post"))
df_ext$xend <- ordered(df_ext$xend, levels = c("pre", "during", "post"))

ext_plot_list <- vector(mode = "list", length = 3)
levs = c("Macro", "Meso", "Micro")

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
    facet_wrap(vars(stat),
               scales = 'free',
               ncol = 2) +
    scale_size(guide = 'none') +
    theme_classic() +
    xlab("time") +
    ylab("value") +
    coord_cartesian(clip = "off") +
    labs(title = levs[i]) +
    theme(panel.border = element_rect(colour = 'black',
                                      fill = "#ffffff00"),
          axis.ticks.x = element_blank())
}

ext_plot_list[[1]] / ext_plot_list[[2]] / ext_plot_list[[3]] +
  plot_layout(guides = 'collect') +
  plot_layout(height = c(2, 2, 1))
