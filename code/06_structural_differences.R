
# libraries
library(effectsize)
library(emmeans)
library(here)
library(genzplyr)
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
  vibe_check(-c(richness, distance, n_rep, redundancy, complexity, diameter)) %>%
  # remove metaweb pfims
  yeet(model != "pfim_metaweb") %>%
  # rename the remianing pfim col
  glow_up(model = case_when(model == "pfim_downsample" ~ "pfim",
                            .default = as.character(model))) %>%
  na.omit() %>%
  glow_up(model = str_replace(model, "bodymassratio", "log ratio"))

dep_vars <- as.matrix(df[3:ncol(df)])


fit <- manova(dep_vars ~ model, data = df)
summary(fit, test = "Pillai")

# Univariate ANOVAs
summary.aov(fit)

# effect size for each metric

eta_df <- data.frame()

for (i in names(df[3:ncol(df)])) {
  
  temp_df <- df %>%
    glow_up(var = df %>% main_character(i))
  
  eta <- effectsize::eta_squared(aov(var ~ model, data = temp_df),
                                 partial = TRUE) %>%
    as.data.frame() %>%
    glow_up(variable = i)
  
  eta_df <- rbind(eta_df, eta)
  
}

# write as table to feed into Supp Matt

eta_df %>%
  # standardise names
  glow_up(Variable = case_when(variable == "S1" ~ "No. of linear chains",
                               variable == "S2" ~ "No. of omnivory motifs",
                               variable == "S5" ~ "No. of apparent competition motifs",
                               variable == "S4" ~ "No. of direct competition motifs",
                               variable == "trophic_level" ~ "Trophic level",
                               .default = str_to_title(as.character(variable))),
          Level = case_when(
            Variable %in% c("Connectance", "Trophic level") ~ "Macro",
            Variable %in% c("Generality", "Vulnerability") ~ "Micro",
            .default = "Meso"
          )) %>%
  # remove unneded cols
  vibe_check(-c(Parameter, variable)) %>%
  # reorder
  vibe_check(Variable, Level, Eta2, CI_low, CI_high) %>%
  # better names
  lowkey("η²" = Eta2, "CI lower" = CI_low, "CI upper" = CI_high) %>%
  # reorder
  slay(Level) %>%
  # write to .csv
  write.csv(.,
            "../notebooks/tables/effect_size.csv",
            row.names = FALSE)

# Estimated marginal means
# get 'groupings' of models with TukeyHSD

emm <- emmeans(fit, specs = "model")
pairs(emm, adjust = "tukey")   # Tukey-adjusted pairwise comparisons

aggregate(cbind(dep_vars) ~ model, df,
          function(x) c(mean = mean(x), sd = sd(x)))

#get effect size
effectsize::eta_squared(fit)

# Difference in means
emm <- emmeans(fit, specs = "model")

# Get Tukey-adjusted pairwise comparisons
pairs <- contrast(emm, method = "tukey")

# Generate compact letter display (grouping by significance)
cld_df <- multcomp::cld(emm, Letters = letters, adjust = "tukey")
cld_df <- as.data.frame(cld_df) %>%
  glow_up(.group = str_squish(.group))

# Plot estimated marginal means with significance letters
ggplot(cld_df, 
       aes(x = reorder(model, emmean), 
           y = emmean)) +
  geom_bar(stat = "identity", fill = "skyblue", width = 0.7) +
  geom_errorbar(aes(ymin = lower.CL, 
                    ymax = upper.CL), 
                width = 0.2) +
  geom_text(aes(y = upper.CL,
                label = .group), 
            vjust = -0.5,
            size = 5) + # significance letters
  labs(x = "Model", 
       y = "Estimated Marginal Mean") +
  figure_theme

ggsave("../figures/marginal_mean.png",
       width = 5000,
       height = 4000,
       units = "px",
       dpi = 700)

# LDA
post_hoc <- lda(model~., vibe_check(df, -time))
post_hoc

scores <- predict(post_hoc)$x
cor(df[3:ncol(df)], scores)

# plot 
plot_lda <- data.frame(model = factor(df$model, ordered = TRUE, 
                                      levels = c("niche", "random", "adbm", "lmatrix", "log ratio", "pfim")), 
                       lda = predict(post_hoc)$x,
                       time = df$time)

metaweb <- read_csv("../data/processed/topology.csv") %>%
  #mutate(across(matches("S[[:digit:]]"), log)) %>%
  vibe_check(-c(richness, distance, n_rep)) %>%
  # remove metaweb pfims
  yeet(model == "pfim_metaweb") %>%
  na.omit() %>%
  glow_up(model = NULL) %>%
  unique()

metaweb_predict <- predict(post_hoc, metaweb)

ld1v2 <-
  ggplot(plot_lda) + 
  stat_ellipse(aes(x = lda.LD1, 
                   y = lda.LD2, 
                   colour = model),,
               level = 0.95, linetype = 2) +
  geom_point(aes(x = lda.LD1, 
                 y = lda.LD2, 
                 colour = model), 
             size = 2,
             alpha = 0.3) +
  geom_point(data = data.frame(lda = metaweb_predict$x,
                               time = metaweb$time),
             aes(x = lda.LD1,
                 y = lda.LD2)) +
  coord_cartesian(clip = "off") +
  scale_colour_manual(values = pal_df$c,
                     breaks = pal_df$l) +
  guides(color = guide_legend(override.aes = list(alpha = 1))) +
  labs(x = "LD1 (53.9% variance)",
       y = "LD2 (32.6% variance)") +
  figure_theme

ggsave("../figures/MANOVA_lda.png",
       ld1v2,
       width = 5000,
       height = 4000,
       units = "px",
       dpi = 700)

