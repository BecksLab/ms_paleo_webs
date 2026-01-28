# libraries
library(genzplyr)
library(ggtext)
library(here)
library(reshape2)
library(tidyverse)

# set path to code sub dir
setwd(here("code"))


# import and clean
beta_df <- read_csv("../data/processed/beta_div.csv") %>%
  # remove metaweb pfims
  yeet(!str_detect(left, "pfim_metaweb")) %>%
  yeet(!str_detect(right, "pfim_metaweb")) %>%
  yeet(β_type == "βOS") %>%
  glow_up(left = str_extract(left, "[a-z]+"),
          right = str_extract(right, "[a-z]+")) %>%
  glow_up(left = str_replace(left, "bodymassratio", "log ratio"),
          right = str_replace(right, "bodymassratio", "log ratio")) %>%
  vibe_check(left, right, β_div)

colnames(beta_df) <- c("Model1", "Model2", "Beta_turnover")

# write raw beta divergence values to CSV for supp mat
write.csv(beta_df, "../notebooks/tables/beta_divergence.csv", row.names = FALSE)

# compare differences

lm_df <- 
  beta_df %>%
  glow_up(combo = paste0(str_extract(Model1, "[a-z]+"), "_", str_extract(Model2, "[a-z]+"))) %>%
  glow_up(combo = case_when(combo == "bodymassratio_adbm" ~ "adbm_bodymassratio",
                            combo == "lmatrix_adbm" ~ "adbm_lmatrix",
                            combo == "pfim_adbm" ~ "adbm_pfim",
                            combo == "pfim_bodymassratio" ~ "bodymassratio_pfim",
                            combo == "random_bodymassratio" ~ "bodymassratio_random",
                            combo == "pfim_lmatrix" ~ "lmatrix_pfim",
                            combo == "bodymassratio_lmatrix" ~ "lmatrix_bodymassratio",
                            .default = as.character(combo))) 

anova(lm(Beta_turnover ~ combo, lm_df))

# write ANOVA results to CSV
anova_res <- anova(lm(Beta_turnover ~ combo, lm_df))
write.csv(as.data.frame(anova_res), "../notebooks/tables/beta_divergence_anova.csv")

# plotting

tile_df <-
  beta_df %>%
  squad_up(Model1, Model2) %>%
  no_cap(Beta_turnover = mean(Beta_turnover, na.rm = TRUE))

ggplot(tile_df, 
       aes(x = Model1, 
           y = Model2, 
           fill = Beta_turnover)) +
  geom_tile(color = "white") +
  scale_fill_viridis_c(option = "plasma", name = "Beta Turnover") +
  theme_minimal() +
  labs(x = NULL,
       y = NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.background = element_rect(fill = "white"))

ggsave("../figures/beta_div.png",
       width = 4000,
       height = 3500,
       units = "px",
       dpi = 600)
