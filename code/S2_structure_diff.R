# call in the analysis code

source("06_structural_differences.R")

# ---- Setup ----
library(tidyverse)
library(effectsize)
library(emmeans)
library(rstatix)
library(MASS)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(multcomp)
library(knitr)
library(kableExtra)

df <-  read_csv("../data/processed/topology.csv") %>%
  #mutate(across(matches("S[[:digit:]]"), log)) %>%
  vibe_check(-c(richness, distance, n_rep, redundancy, complexity, diameter)) %>%
  # remove metaweb pfims
  yeet(model != "pfim_metaweb") %>%
  # rename the remianing pfim col
  glow_up(model = case_when(model == "pfim_downsample" ~ "pfim",
                            .default = as.character(model))) %>%
  na.omit() %>%
  glow_up(model = str_replace(model, "bodymassratio", "log ratio"))

## Table S1 – Descriptive statistics

desc_stats <- df %>%
  group_by(model) %>%
  summarise(across(connectance:S5, list(mean=mean, sd=sd))) %>%
  ungroup()

write_csv(desc_stats, "../notebooks/tables/Table_S1_descriptive_stats.csv")

# Table S2 – Univariate ANOVA with partial η²

# Univariate ANOVAs
anova_list <- summary.aov(fit)

eta_df <- data.frame()
for (i in names(df[3:ncol(df)])) {
  temp_df <- df %>% glow_up(var = df %>% main_character(i))
  eta <- effectsize::eta_squared(aov(var ~ model, data=temp_df), partial=TRUE) %>%
    as.data.frame() %>%
    glow_up(variable=i)
  eta_df <- rbind(eta_df, eta)
}

anova_table <- map_df(anova_list, ~as.data.frame(.x[1,])) %>%
  mutate(Variable = names(anova_list)) %>%
  left_join(eta_df %>% select(variable, Eta2, CI_low, CI_high), by=c("Variable"="variable")) %>%
  select(Variable, `F value`, Df, `Pr(>F)`) %>%
  rename(`F (df1, df2)`=`F value`, `p-value`=`Pr(>F)`)

write_csv(anova_table, "../notebooks/tables/Table_S2_ANOVA_eta.csv")

# Table S3 – Tukey pairwise comparisons

emm <- emmeans(fit, specs = "model")
pairs_df <- as.data.frame(pairs(emm, adjust="tukey"))

write_csv(pairs_df, "../notebooks/tables/Table_S3_Tukey_pairs.csv")

# Table S4 – Compact letter display

cld_df <- multcomp::cld(emm, Letters=letters, adjust="sidak")
cld_df <- as.data.frame(cld_df) %>% glow_up(.group = str_squish(.group))

write_csv(cld_df %>% 
            select(model, emmean, SE, .group), "../notebooks/tables/Table_S4_CLD.csv")

# Table S5 – LDA coefficients

post_hoc <- lda(model ~ ., vibe_check(df, -time))
coef_df <- as.data.frame(post_hoc$scaling) %>% rownames_to_column("Variable")

write_csv(coef_df, "../notebooks/tables/Table_S5_LDA_coeffs.csv")

# Table S6 – Correlation of variables with LD axes

scores <- predict(post_hoc)$x
cor_df <- cor(df[3:ncol(df)], scores) %>% as.data.frame() %>% rownames_to_column("Variable")

write_csv(cor_df, "../notebooks/tables/Table_S6_LD_correlations.csv")

# Figure S1 – Estimated marginal means

png("../notebooks/figures/Figure_S1_EMMs.png", width=5000, height=4000, res=700)

ggplot(cld_df, 
       aes(x=reorder(model, emmean), 
           y=emmean)) +
  geom_bar(stat="identity", 
           fill="skyblue", width=0.7) +
  geom_errorbar(aes(ymin=emmean-SE, 
                    ymax=emmean+SE), 
                width=0.2) +
  geom_text(aes(y=emmean+SE, 
                label=.group), 
            vjust=-0.5, size=6) +
  labs(x="Model", 
       y="Estimated Marginal Mean") +
  theme_minimal()
dev.off()

# Figure S2 – LDA scatterplots

library(patchwork)
library(ggrepel)
library(grid)

# ---- Prepare correlation data ----
# Use absolute correlation magnitude to scale arrows
cor_df_scaled <- cor(dep_vars, predict(post_hoc)$x[,1:3]) %>% as.data.frame()
cor_df_scaled$Variable <- rownames(cor_df_scaled)

# Function to rescale arrows to fit nicely in plot
scale_arrow <- function(x, factor=1.5) x * factor

cor_df_scaled <- cor_df_scaled %>%
  mutate(LD1 = scale_arrow(LD1),
         LD2 = scale_arrow(LD2),
         LD3 = scale_arrow(LD3))

# Calculate proportion of variance explained by each LD axis
# MASS::lda stores singular values in post_hoc$svd
var_explained <- post_hoc$svd^2 / sum(post_hoc$svd^2) * 100  # %
names(var_explained) <- paste0("LD", 1:length(var_explained))
var_explained

lda_scores <- as.data.frame(predict(post_hoc)$x)  # LD1, LD2, LD3...
plot_lda <- cbind(
  model = factor(df$model, ordered = TRUE,
                 levels = c("niche", "random", "adbm", "lmatrix", "log ratio", "pfim")),
  lda_scores,
  time = df$time
)

# ---- Function to make correlation biplot for any combination of axes ----
make_corr_plot <- function(x_axis, y_axis) {
  ggplot(cor_df_scaled) +
    geom_hline(yintercept=0, linetype="dashed", color="grey70") +
    geom_vline(xintercept=0, linetype="dashed", color="grey70") +
    annotate("path", x=cos(seq(0,2*pi,length.out=200)), 
             y=sin(seq(0,2*pi,length.out=200)), color="grey50") +
    geom_segment(aes_string(x=0, y=0, xend=x_axis, yend=y_axis),
                 arrow=arrow(length=unit(0.15,"cm")), color="skyblue") +
    geom_text_repel(aes_string(x=x_axis, y=y_axis, label="Variable"), size=3) +
    coord_equal() + theme_classic() +
    labs(x=x_axis, y=y_axis)
}

# ---- LDA scatterplots ----
make_lda_plot <- function(x_axis, y_axis) {
  ggplot(plot_lda, aes(x=!!sym(x_axis), y=!!sym(y_axis), colour=model)) +
    stat_ellipse(level=0.95, linetype=2) +
    geom_point(size=2, alpha=0.3) +
    scale_colour_manual(values=pal_df$c, breaks=pal_df$l) +
    labs(x=paste0(x_axis, " (", round(var_explained[x_axis],1),"%)"),
         y=paste0(y_axis, " (", round(var_explained[y_axis],1),"%)")) +
    figure_theme
}

# ---- Generate all panels ----
lda1v2 <- make_lda_plot("LD1", "LD2")
lda1v3 <- make_lda_plot("LD1", "LD3")
lda2v3 <- make_lda_plot("LD3", "LD2")

corr1v2 <- make_corr_plot("LD1", "LD2")
corr1v3 <- make_corr_plot("LD1", "LD3")
corr2v3 <- make_corr_plot("LD3", "LD2")

# ---- Combine panels using patchwork ----
combined_panel <- (lda1v2 + corr1v2) /
  (lda1v3 + corr1v3) /
  (lda2v3 + corr2v3) +
  plot_layout(guides="collect") & theme(legend.position="bottom")

# ---- Save final multi-panel figure ----
ggsave("../notebooks/figures/MANOVA_LDA_biplot_panel.png",
       combined_panel,
       width=6000, height=6300, units="px", dpi=700)
