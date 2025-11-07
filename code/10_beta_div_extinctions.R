# libraries
library(ggtext)
library(here)
library(tidyverse)

# set path to code sub dir
setwd(here("code"))

# import simulated data

df <- read_csv("../data/processed/extinction_betadiv.csv") %>%
  select(-n_rep) %>%
  rbind(.,
        read_csv("../data/processed/extinction_betadiv.csv") %>%
          filter(β_type != "βS") %>%
          pivot_wider(names_from = β_type,
                      values_from = β_div) %>%
          mutate(βST = βWN - βOS) %>%
          select(-c(βWN, βOS, n_rep)) %>%
          pivot_longer(!c(model, extinction_mechanism),
                       names_to = "β_type",
                       values_to = "β_div")) %>%
  # get mean (for now might be worth looking at sd as well...)
  group_by(model, extinction_mechanism, β_type) %>%
  summarise(mean = mean(β_div, na.rm = TRUE)) 
  

ggplot(df,
       aes(x = extinction_mechanism,
           y = model)) + 
  geom_tile(aes(fill = mean), colour = "grey50") +
  scale_fill_distiller(palette = "Spectral", direction = 1) +
  facet_wrap(vars(β_type)) +
  theme_classic() +
  theme(panel.border = element_rect(colour = 'black',
                                    fill = "#ffffff00"),
        axis.text.x = element_markdown(angle = 45, hjust=1))

ggsave("../figures/beta_div_extinctions.png",
       width = 5000,
       height = 4000,
       units = "px",
       dpi = 600)
