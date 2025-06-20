library(here)
library(effectsize)
library(MASS)
library(patchwork)
library(rstatix)
library(tidyverse)

# set path to code sub dir
setwd(here("code"))

#load script that determines plotting aesthetics
source("lib/plotting_theme.R")

# import data
metaweb <- read_csv("../data/processed/topology.csv") %>%
  filter(model == "pfim_metaweb") %>%
  select(-c(distance, n_rep)) %>%
  # remove metaweb pfims
  pivot_longer(
    cols = -c(time, model),
    names_to = "stat",
    values_to = "stat_val") %>%
  group_by(time, stat) %>%
  summarise(mean = mean(stat_val, na.rm = TRUE))

read_csv("../data/dunhill/topology.csv") %>%
  select(-c(distance, n_rep)) %>%
  # remove metaweb pfims
  pivot_longer(
    cols = -c(y_val, time),
    names_to = "stat",
    values_to = "stat_val") %>%
  group_by(time, y_val, stat) %>% 
  summarise(mean = mean(stat_val, na.rm = TRUE)) %>%
  ggplot(.,
         aes(y = mean,
             x = y_val,
             group = time,
             colour = time)) +
  geom_vline(xintercept = 2.5,
             alpha = 0.2) +
  geom_line() +
  geom_hline(data = metaweb,
             aes(yintercept = mean,
                 colour = time),
             linetype = 3) +
  facet_wrap(vars(stat),
             scales = 'free',
             ncol = 2) +
  scale_size(guide = 'none') +
  xlab("y value") +
  ylab("value") +
  coord_cartesian(clip = "off") +
  figure_theme

ggsave("../figures/pfim_downsample.png",
       width = 5000,
       height = 6000,
       units = "px",
       dpi = 600)


df <- read_csv("../data/dunhill/topology.csv") %>%
  #mutate(across(matches("S[[:digit:]]"), log)) %>%
  select(-c(richness, distance, n_rep, trophic_level)) %>%
  na.omit()

dep_vars <- as.matrix(df[3:ncol(df)])


fit <- manova(dep_vars ~ y_val + time, data = df)
summary(fit)

#get effect size
effectsize::eta_squared(fit)

post_hoc <- lda(y_val~., df[2:ncol(df)])
post_hoc

# plot 
plot_lda <- data.frame(y_val = as.numeric(as.character(predict(post_hoc)$class)),
                       lda = predict(post_hoc)$x,
                       time = df$time)

# predict real metaweb location
metaweb <- read_csv("../data/processed/topology.csv") %>%
  filter(model == "pfim_metaweb") %>%
  select(-c(richness, distance, n_rep)) %>%
  mutate(model = NULL) %>%
  unique()

metaweb_predict <- predict(post_hoc, metaweb)


ggplot(plot_lda) + 
  geom_point(aes(x = lda.LD1, 
                 y = lda.LD2, 
                 colour = y_val), 
             size = 3,
             alpha = 0.1) +
  scale_colour_distiller(palette = "YlGnBu") +
  geom_point(data = data.frame(lda = metaweb_predict$x,
                               time = metaweb$time),
             aes(x = lda.LD1,
                 y = lda.LD2)) +
  facet_wrap(vars(time)) +
  figure_theme

ggsave("../figures/pfim_downsample_lda.png",
       width = 10000,
       height = 6000,
       units = "px",
       dpi = 600)
