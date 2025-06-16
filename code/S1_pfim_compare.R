library(tidyverse)

# set path to code sub dir
setwd(here("code"))

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
  xlab("time") +
  ylab("value") +
  coord_cartesian(clip = "off") +
  theme_classic() +
  theme(panel.border = element_rect(colour = 'black',
                                    fill = "#ffffff00"),
        axis.ticks.x = element_blank())

ggsave("../figures/pfim_downsample.png",
       width = 5000,
       height = 6000,
       units = "px",
       dpi = 600)
