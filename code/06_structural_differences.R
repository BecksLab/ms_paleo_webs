# libraries: Statistical and Plotting ecosystem
library(candisc)
library(effectsize) # For eta-squared calculations
library(emmeans)    # For estimated marginal means
library(here)       # For robust file path management
library(genzplyr)   # User-specific utility package (wrappers for tidyverse)
library(ggrepel)
library(ggtext)
library(MASS)       # For Linear Discriminant Analysis (LDA)
library(patchwork)
library(rstatix)
library(tidyverse)

# Set path to code sub directory for consistent relative paths
setwd(here("code"))

# Load script that determines plotting aesthetics (colors, fonts, etc.)
source("lib/plotting_theme.R")

# 1. Data Cleaning & Preparation
df <- read_csv("../data/processed/topology.csv") %>%
  # Remove non-numeric or irrelevant metadata columns for MANOVA
  vibe_check(-c(richness, distance, n_rep, redundancy, complexity, diameter)) %>%
  # Focus only on specific model variants (excludes the pfim_metaweb)
  yeet(model != "pfim_metaweb") %>%
  # Standardize model names for publication-quality legends
  glow_up(model = case_when(model == "pfim_downsample" ~ "PFIM",
                            model == "bodymassratio" ~ "Body-size ratio",
                            model == "adbm" ~ "ADBM",
                            model == "lmatrix" ~ "ATN",
                            .default = as.character(model)),
          trophic_level = round(trophic_level, 0)) %>%
  na.omit()

# Separate the dependent variables (topological metrics) into a matrix
dep_vars <- as.matrix(df[3:ncol(df)])

# 2. Multivariate Statistical Testing
# Does the choice of 'model' significantly affect the collective suite of network metrics?
fit <- manova(dep_vars ~ model, data = df)
summary(fit, test = "Pillai")


# 3: Canonical discriminant analysis
cda <- candisc(fit)

# Eigenvalues and variance explained
eig_df <- data.frame(
  CV = paste0("CV", seq_along(cda$eig)),
  Eigenvalue = cda$eig,
  Variance = 100 * cda$eig / sum(cda$eig)
)

write.csv(eig_df,
          "../notebooks/tables/cda_eigenvalues.csv",
          row.names = FALSE)

loadings_df <- as.data.frame(cda$structure[, 1:2]) %>%
  rownames_to_column("metric") %>%
  glow_up(CV1 = Can1,
          CV2 = Can2) %>%
  mutate(
    abs_loading = pmax(abs(CV1), abs(CV2)),
    keep = abs_loading >= 0.4
  ) %>%
  mutate(Level = case_when(
    metric %in% c("connectance", "trophic_level") ~ "Macro",
    metric %in% c("generality", "vulnerability") ~ "Micro",
    .default = "Meso"
  ))

loading_plot <- ggplot(loadings_df, aes(x = CV1, y = CV2)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey70") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey70") +
  geom_segment(
    aes(x = 0, y = 0, xend = CV1, yend = CV2),
    arrow = arrow(length = unit(0.02, "npc")),
    linewidth = 0.8,
    alpha = 0.8
  ) +
  geom_segment(
    aes(x = 0, y = 0, xend = CV1, yend = CV2, colour = Level),
    arrow = arrow(length = unit(0.02, "npc")),
    linewidth = 0.9
  ) +
  scale_colour_manual(values = c(
    Macro = "#1b9e77",
    Meso  = "#d95f02",
    Micro = "#7570b3"
  )) +
  geom_text_repel(
    data = subset(loadings_df, keep),
    aes(label = metric),
    size = 4,
    box.padding = 0.4,
    point.padding = 0.3,
    max.overlaps = Inf
  ) +
  coord_equal(xlim = c(-1, 1), ylim = c(-1, 1)) +
  labs(
    x = "CV1",
    y = "CV2"
  ) +
  figure_theme

ggsave("../figures/loading_plot.png", 
       loading_plot, 
       width = 5000, height = 4000, 
       units = "px", dpi = 700)

# 4: Canonical structure coefficients

structure_df <- as.data.frame(cda$structure) %>%
  rownames_to_column("Metric")

write.csv(structure_df,
          "../notebooks/tables/cda_structure_coefficients.csv",
          row.names = FALSE)

# 5: LDA for visualization only

post_hoc <- lda(model ~ ., data = vibe_check(df, -time))
scores <- as.data.frame(predict(post_hoc)$x)

scores <- scores %>%
  glow_up(
    CV1 = LD1,
    CV2 = LD2,
    model = df$model,
    time = df$time
  )

plot_lda <- data.frame(model = factor(df$model, ordered = TRUE, 
                                      levels = c("niche", "random", "ADBM", "ATN", "Body-size ratio", "PFIM")), 
                       lda = predict(post_hoc)$x,
                       time = df$time)

# Extract proportion of trace for axis labels
lda_variance <- round((post_hoc$svd^2 / sum(post_hoc$svd^2)) * 100, 1)
ld1_lab <- paste0("LD1 (", lda_variance[1], "% variance)")
ld2_lab <- paste0("LD2 (", lda_variance[2], "% variance)")

# Overlay the "metaweb" versions (non-downsampled) as single points in the same space
metaweb <- read_csv("../data/processed/topology.csv") %>%
  vibe_check(-c(richness, distance, n_rep)) %>%
  yeet(model == "pfim_metaweb") %>%
  na.omit() %>%
  glow_up(model = NULL) %>%
  unique()

metaweb_predict <- predict(post_hoc, metaweb)

# Final Plot: LD1 vs LD2 with 95% confidence ellipses
ld1v2 <- ggplot(plot_lda) + 
  stat_ellipse(aes(x = lda.LD1, 
                   y = lda.LD2, 
                   colour = model), 
               level = 0.95, linetype = 2) +
  geom_point(aes(x = lda.LD1, 
                 y = lda.LD2, 
                 colour = model), 
             size = 2, alpha = 0.3) +
  geom_point(data = data.frame(lda = metaweb_predict$x, time = metaweb$time),
             aes(x = lda.LD1, 
                 y = lda.LD2)) +
  stat_ellipse(aes(x = lda.LD1, 
                   y = lda.LD2, 
                   colour = model), 
               level = 0.95, linetype = 2) +
  coord_cartesian(clip = "off") +
  scale_colour_manual(values = pal_df$c, 
                      breaks = pal_df$l) +
  labs(x = ld1_lab,
       y = ld2_lab) +
  figure_theme


ggsave("../figures/MANOVA_lda.png", 
       ld1v2, 
       width = 5000, height = 4000, 
       units = "px", dpi = 700)

# 6: Convex hull canonical plot

hull_df <- cda_scores %>%
  group_by(model) %>%
  slice(chull(CV1, CV2)) %>%
  ungroup()

centroids <- cda_scores %>%
  group_by(model) %>%
  summarise(
    CV1 = mean(CV1),
    CV2 = mean(CV2)
  )

plot_cda_hull <- ggplot(cda_scores, aes(x = CV1, y = CV2, colour = model)) +
  
  # convex hulls
  geom_polygon(
    data = hull_df,
    aes(x = CV1, y = CV2, fill = model),
    alpha = 0.15,
    colour = NA
  ) +
  geom_polygon(
    data = hull_df,
    aes(x = CV1, y = CV2, colour = model),
    fill = NA,
    linewidth = 0.8
  ) +
  geom_point(
    size = 2,
    alpha = 0.35
  )  +
  geom_point(
    data = centroids,
    aes(x = CV1, y = CV2),
    shape = 21,
    fill = "white",
    size = 3,
    stroke = 1
  ) +
  coord_cartesian(clip = "off") +
  scale_colour_manual(values = pal_df$c, breaks = pal_df$l) +
  scale_fill_manual(values = pal_df$c, breaks = pal_df$l) +
  labs(
    x = cv1_lab,
    y = cv2_lab
  ) +
  figure_theme

ggsave("../figures/cda_hull.png", 
       plot_cda_hull, 
       width = 5000, height = 4000, 
       units = "px", dpi = 700)
