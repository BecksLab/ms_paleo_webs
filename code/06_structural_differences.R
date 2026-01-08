
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
  vibe_check(-c(richness, distance, n_rep, redundancy)) %>%
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
            Variable %in% c("Complexity", "Connectance", "Trophic level", "Diameter", "Redundancy") ~ "Macro",
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
  labs(x = "LD1 (53% variance)",
       y = "LD2 (32% variance)") +
  figure_theme

ggsave("../figures/MANOVA_lda.png",
       ld1v2,
       width = 5000,
       height = 4000,
       units = "px",
       dpi = 700)

# other LDA combos
ld1v3 <-
  ggplot(plot_lda) + 
  stat_ellipse(aes(x = lda.LD1, 
                   y = lda.LD3, 
                   colour = model),,
               level = 0.95, linetype = 2) +
  geom_point(aes(x = lda.LD1, 
                 y = lda.LD3, 
                 colour = model), 
             size = 2,
             alpha = 0.3) +
  geom_point(data = data.frame(lda = metaweb_predict$x,
                               time = metaweb$time),
             aes(x = lda.LD1, 
                 y = lda.LD3)) +
  coord_cartesian(clip = "off") +
  scale_colour_manual(values = pal_df$c,
                      breaks = pal_df$l) +
  guides(color = guide_legend(override.aes = list(alpha = 1))) +
  labs(x = "LD1 (52% variance)",
       y = "LD3 (10% variance)") +
  figure_theme

ld2v3 <-
  ggplot(plot_lda) + 
  stat_ellipse(aes(x = lda.LD3, 
                   y = lda.LD2, 
                   colour = model),,
               level = 0.95, linetype = 2) +
  geom_point(aes(x = lda.LD3, 
                 y = lda.LD2, 
                 colour = model), 
             size = 2,
             alpha = 0.3) +
  geom_point(data = data.frame(lda = metaweb_predict$x,
                               time = metaweb$time),
             aes(x = lda.LD3, 
                 y = lda.LD2)) +
  coord_cartesian(clip = "off") +
  scale_colour_manual(values = pal_df$c,
                      breaks = pal_df$l) +
  guides(color = guide_legend(override.aes = list(alpha = 1))) +
  labs(y = "LD2 (32% variance)",
       x = "LD3 (10% variance)") +
  figure_theme

# going to do the same with the corr plots

# Get the linear discriminant scores
lda_values <- predict(post_hoc)

# Correlation of original variables with the linear discriminants
correlations <- cor(dep_vars, lda_values$x)

# 4. Convert to a data frame for plotting
corr_df <- as.data.frame(correlations)
corr_df$Variable <- rownames(corr_df)

# 5. Create correlation circle plot
ld1v2_corr <- ggplot(corr_df) +
  geom_hline(yintercept = 0, 
             linetype = "dashed", 
             color = "grey70") +
  geom_vline(xintercept = 0, 
             linetype = "dashed", 
             color = "grey70") +
  # Add a unit circle
  annotate("path",
           x = cos(seq(0, 2 * pi, length.out = 200)),
           y = sin(seq(0, 2 * pi, length.out = 200)),
           color = "grey50")  +
  geom_segment( 
    aes(x = 0,
        y = 0,
        xend = LD1, 
        yend = LD2),
    arrow = arrow(length = unit(0.1,"cm")),
    color = "skyblue") +
  geom_text_repel(
    aes(x = LD1, 
        y = LD2, 
        label = Variable),
    max.overlaps = getOption("ggrepel.max.overlaps", default = 100), 
    size = 3) +
  coord_equal()+
  labs(
    x = "LD1",
    y = "LD2"
  ) +
  theme_classic()

ld1v3_corr <- ggplot(corr_df) +
  geom_hline(yintercept = 0, 
             linetype = "dashed", 
             color = "grey70") +
  geom_vline(xintercept = 0, 
             linetype = "dashed", 
             color = "grey70") +
  # Add a unit circle
  annotate("path",
           x = cos(seq(0, 2 * pi, length.out = 200)),
           y = sin(seq(0, 2 * pi, length.out = 200)),
           color = "grey50")  +
  geom_segment( 
    aes(x = 0,
        y = 0,
        xend = LD1, 
        yend = LD3),
    arrow = arrow(length = unit(0.1,"cm")),
    color = "skyblue") +
  geom_text_repel(
    aes(x = LD1, 
        y = LD3, 
        label = Variable),
    max.overlaps = getOption("ggrepel.max.overlaps", default = 100), 
    size = 3) +
  coord_equal()+
  labs(
    x = "LD1",
    y = "LD3"
  ) +
  theme_classic()

ld2v3_corr <- ggplot(corr_df) +
  geom_hline(yintercept = 0, 
             linetype = "dashed", 
             color = "grey70") +
  geom_vline(xintercept = 0, 
             linetype = "dashed", 
             color = "grey70") +
  # Add a unit circle
  annotate("path",
           x = cos(seq(0, 2 * pi, length.out = 200)),
           y = sin(seq(0, 2 * pi, length.out = 200)),
           color = "grey50")  +
  geom_segment( 
    aes(x = 0,
        y = 0,
        xend = LD3, 
        yend = LD2),
    arrow = arrow(length = unit(0.1,"cm")),
    color = "skyblue") +
  geom_text_repel(
    aes(x = LD3, 
        y = LD2, 
        label = Variable),
    max.overlaps = getOption("ggrepel.max.overlaps", default = 100), 
    size = 3) +
  coord_equal()+
  labs(
    x = "LD3",
    y = "LD2"
  ) +
  theme_classic()


# patchwork them up

layout <- "
AD
BE
CF
"

(ld1v2 + labs(tag = "A")) + (ld1v3 + labs(tag = "B")) + (ld2v3 + labs(tag = "C")) + 
  ld1v2_corr + ld1v3_corr + ld2v3_corr +
  plot_layout(guides = 'collect',
              design = layout) &
  theme(legend.position='bottom')

ggsave("../figures/MANOVA_3_panel.png",
       width = 6000,
       height = 6300,
       units = "px",
       dpi = 700)


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

