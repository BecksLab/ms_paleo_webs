library(tidyverse)
library(viridis)
library(reshape2)
library(genzplyr)
library(here)

# set working directory
setwd(here("code"))

# 1. Load beta_divergence data
beta_df <- read_csv("../data/processed/beta_div.csv") %>%
  # remove metaweb PFIMs and βOS
  filter(!str_detect(left, "pfim_metaweb"),
         !str_detect(right, "pfim_metaweb"),
         β_type != "βOS") %>%
  mutate(
    left = str_replace(left, "bodymassratio", "log ratio"),
    right = str_replace(right, "bodymassratio", "log ratio")
  ) %>%
  select(Model1 = left, Model2 = right, Beta_turnover = β_div)

# 2. Compute mean and SD of beta_turnover for each model pair
summary_df <- beta_df %>%
  group_by(Model1, Model2) %>%
  summarise(
    Beta_mean = mean(Beta_turnover, na.rm = TRUE),
    Beta_sd   = sd(Beta_turnover, na.rm = TRUE),
    .groups = "drop"
  )

# 3. Convert to wide matrices for plotting
mean_matrix <- summary_df %>%
  pivot_wider(names_from = Model2, values_from = Beta_mean)

sd_matrix <- summary_df %>%
  pivot_wider(names_from = Model2, values_from = Beta_sd)

# 4. Prepare for ggplot (long format)
mean_long <- summary_df %>% select(Model1, Model2, Beta_mean)
sd_long   <- summary_df %>% select(Model1, Model2, Beta_sd)

# 5. Plot mean beta turnover
p_mean <- ggplot(mean_long, aes(x = Model1, y = Model2, fill = Beta_mean)) +
  geom_tile(color = "white") +
  scale_fill_viridis_c(option = "plasma", name = "Mean β-turnover") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = NULL, y = NULL) +
  ggtitle("Mean Pairwise β-turnover")

# 6. Plot SD of beta turnover
p_sd <- ggplot(sd_long, aes(x = Model1, y = Model2, fill = Beta_sd)) +
  geom_tile(color = "white") +
  scale_fill_viridis_c(option = "magma", name = "SD β-turnover") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = NULL, y = NULL) +
  ggtitle("SD of Pairwise β-turnover")

# 7. Combine into a single panel using patchwork
library(patchwork)
fig_beta <- p_mean + p_sd + plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "bottom")

# 8. Save figure
ggsave("../notebooks/figures/beta_div_mean_sd.png",
       fig_beta,
       width = 8000, height = 4000,
       units = "px", dpi = 600)
