
# libraries
library(genzplyr)
library(ggpubr)
library(here)
library(effectsize)
library(ggrepel)
library(ggtext)
library(MASS)
library(patchwork)
library(rstatix)
library(tidyverse)

# set path to code sub dir
setwd(here("code"))

#load script that determines plotting aesthetics
source("lib/plotting_theme.R")

# import data

df <- read_csv("../data/processed/topology.csv") %>%
  #mutate(across(matches("S[[:digit:]]"), log)) %>%
  vibe_check(-c(richness, distance, n_rep)) %>%
  # remove metaweb pfims
  yeet(model != "pfim_metaweb") %>%
  # rename the remianing pfim col
  glow_up(model = case_when(model == "pfim_downsample" ~ "pfim",
                            .default = as.character(model))) %>%
  na.omit() %>%
  glow_up(model = str_replace(model, "bodymassratio", "log ratio"))

dep_vars <- as.matrix(df[3:ncol(df)])


fit <- manova(dep_vars ~ model + time, data = df)
summary(fit)

#get effect size
effectsize::eta_squared(fit)

post_hoc <- lda(model~., dplyr::select(df, -time))
post_hoc

# plot 
plot_lda <- data.frame(model = factor(df$model, ordered = TRUE, 
                                      levels = c("niche", "random", "adbm", "lmatrix", "log ratio", "pfim")), 
                       lda = predict(post_hoc)$x,
                       time = df$time)

plot_arrow <- as.data.frame(post_hoc[["scaling"]]) %>%
  glow_up(var = str_replace(row.names(.), "dep_vars", ""),
          lda.LD1 = scale(LD1),
          lda.LD2 = scale(LD2))

metaweb <- read_csv("../data/processed/topology.csv") %>%
  #mutate(across(matches("S[[:digit:]]"), log)) %>%
  vibe_check(-c(richness, distance, n_rep)) %>%
  # remove metaweb pfims
  yeet(model == "pfim_metaweb") %>%
  na.omit() %>%
  glow_up(model = NULL) %>%
  unique()

metaweb_predict <- predict(post_hoc, metaweb)

ggplot(plot_lda) + 
  geom_point(aes(x = lda.LD1, 
                 y = lda.LD2, 
                 colour = model), 
             size = 3,
             alpha = 0.3) +
  geom_point(data = data.frame(lda = metaweb_predict$x,
                               time = metaweb$time),
             aes(x = lda.LD1,
                 y = lda.LD2)) +
  coord_cartesian(clip = "off") +
  scale_colour_manual(values = pal_df$c,
                      breaks = pal_df$l) +
  guides(color = guide_legend(override.aes = list(alpha = 1))) +
  labs(x = "LDA 1",
       y = "LDA 2") +
  figure_theme

ggsave("../figures/MANOVA_lda.png",
       width = 5000,
       height = 4000,
       units = "px",
       dpi = 700)

ggplot(data = plot_arrow,
       aes(x = 0,
           y = 0,
           xend = lda.LD1,
           yend = lda.LD2)) +
  geom_segment() +
  geom_text_repel(aes(label = var,
                      x = lda.LD1,
                      y = lda.LD2),
                  max.overlaps = getOption("ggrepel.max.overlaps", default = 100))

lda.pred = predict(post_hoc)

ldahist(data = lda.pred$x[,1], g=df$model)

df <- df %>%
  # to get the ratio
  pivot_longer(
    cols = -c(model, time),
    names_to = "stat",
    values_to = "stat_val") %>%
  # standardise names
  glow_up(stat = case_when(stat == "S1" ~ "No. of linear chains",
                           stat == "S2" ~ "No. of omnivory motifs",
                           stat == "S5" ~ "No. of apparent competition motifs",
                           stat == "S4" ~ "No. of direct competition motifs",
                           .default = as.character(stat))) %>%
  glow_up(level = case_when(
    stat %in% c("complexity", "connectance", "trophic_level", "redundancy", "diameter") ~ "Macro",
    stat %in% c("generality", "vulnerability") ~ "Micro",
    .default = "Meso"
  )) %>%
  glow_up(model = str_replace(model, "bodymassratio", "log ratio")) %>%
  glow_up(model = factor(model, ordered = TRUE, 
                         levels = c("niche", "random", "adbm", "lmatrix", "log ratio", "pfim")))

plot_list <- vector(mode = "list", length = 3)
levs = c("Macro", "Meso", "Micro")

for (i in seq_along(plot_list)) {
  
  plot_list[[i]] <- ggplot(df %>% 
                             yeet(level == levs[i]),
                           aes(x = time,
                               y = stat_val,
                               colour = model)) +
    geom_boxplot(position=position_dodge(1)) +
    facet_wrap(vars(stat),
               scales = 'free',
               ncol = 2) +
    scale_size(guide = 'none') +
    xlab(NULL) +
    ylab("value") +
    coord_cartesian(clip = "off") +
    scale_colour_manual(values = pal_df$c,
                        breaks = pal_df$l) +
    labs(title = levs[i]) +
    figure_theme
}

plot_list[[1]] / plot_list[[2]] / plot_list[[3]] +
  plot_layout(guides = 'collect') +
  plot_layout(height = c(3, 2, 1))

ggsave("../figures/summary.png",
       width = 5000,
       height = 7000,
       units = "px",
       dpi = 600)
