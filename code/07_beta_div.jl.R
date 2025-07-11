# libraries
library(ggtext)
library(here)
library(tidyverse)

# set path to code sub dir
setwd(here("code"))

# import simulated data

df <- read_csv("../data/processed/beta_div.csv") %>%
  # remove metaweb pfims
  filter(!str_detect(left, "pfim_metaweb")) %>%
  filter(!str_detect(right, "pfim_metaweb")) %>%
  # rename the remianing pfim col
  mutate(left = str_replace(left, "_downsample", ""),
         right = str_replace(right, "_downsample", "")) %>%
  # remove all instances where niche model is present (not a fair comparison?)
  filter(!str_detect(left, "niche")) %>%
  filter(!str_detect(right, "niche")) %>%
  filter(!str_detect(right, "random")) %>%
  filter(!str_detect(left, "random")) %>%
  mutate(combo = paste0(str_extract(left, "[a-z]+"), "_", str_extract(right, "[a-z]+"))) %>%
  mutate(combo = case_when(combo == "bodymassratio_adbm" ~ "adbm_bodymassratio",
                           combo == "lmatrix_adbm" ~ "adbm_lmatrix",
                           combo == "pfim_adbm" ~ "adbm_pfim",
                           combo == "pfim_bodymassratio" ~ "bodymassratio_pfim",
                           combo == "random_bodymassratio" ~ "bodymassratio_random",
                           combo == "pfim_lmatrix" ~ "lmatrix_pfim",
                           combo == "bodymassratio_lmatrix" ~ "lmatrix_bodymassratio",
                           .default = as.character(combo))) %>%
  mutate(shared_pct = β_int_shared/((links_left-β_int_shared) + (links_right-β_int_shared) + β_int_shared),
         model_combo = case_when(combo %in% c("pfim_pfim", "lmatrix_lmatrix", "bodymassratio_bodymassratio",
                                              "adbm_adbm") ~ "Within model",
                                 .default = "Between model")) %>%
  mutate(combo = case_when(model_combo == "Within model" ~ str_extract(combo, "^.+?(?=_)"),
                           .default = combo)) %>%
  filter(model_combo == "Between model") %>%
  mutate(combo = str_replace(combo, "_", "-"),
         combo = str_replace(combo, "bodymassratio", "log ratio"))

ggplot(df,
       aes(x = combo,
           y = shared_pct*100)) + 
  geom_boxplot(position = position_dodge(1),
               outliers = FALSE,
               fill = NA) +
  geom_jitter(data = df %>% group_by(combo) %>% slice_sample(n=100),
              alpha = 0.2) +
  xlab(NULL) +
  ylab("% of interactions shared") +
  coord_cartesian(clip = "off") +
  theme_classic() +
  theme(panel.border = element_rect(colour = 'black',
                                    fill = "#ffffff00"),
        axis.text.x = element_markdown(angle = 45, hjust=1))

ggsave("../figures/beta_div.png",
       width = 4000,
       height = 4000,
       units = "px",
       dpi = 600)
