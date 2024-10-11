library(readr)
library(scales)
library(tidyverse)

#### Structure ####

df = list.files(path = "data/processed/", pattern = ".csv", full.names = TRUE) %>% 
  lapply(read_csv) %>% 
  bind_rows %>% 
        select(-c(network, richness)) %>% 
        # to get the ratio
        mutate(ratio = top/basal,
                top = NULL,
                basal = NULL) %>%
        pivot_longer(
            cols = -c(id, model), 
            names_to = "stat",
            values_to = "stat_val")  %>% 
    # standardise names
        mutate(id = case_when(
            str_detect(id, "^.*pre.*$") ~ "pre",
            str_detect(id, "^.*during.*$") ~ "post",
            str_detect(id, "^.*post.*$") ~ "during",
            TRUE ~ as.character(id)))

df$id <- ordered(df$id, levels=c("pre", "during", "post"))

ggplot(df,
        aes(x = factor(`id`), 
            y = stat_val, 
            colour = model,
            group = model)) +
    geom_line() +
    geom_point() +
    facet_wrap(vars(stat),
                scales = 'free') +
    scale_size(guide = 'none') +
    theme_classic() +
    xlab("time") +
    ylab("value") +
    theme(panel.border = element_rect(colour = 'black',
                                      fill = "#ffffff00"),
            axis.ticks.x = element_blank())

ggsave("figures/summary.png",
       width = 9000,
       height = 5000,
       units = "px",
       dpi = 600)

#### Extinctions ####

df_ext <- read_csv("data/processed/extinctions/extinctions.csv") %>% 
        select(-c(id, links, richness, extinction_mechanism)) %>% 
        # to get the ratio
        mutate(ratio = top/basal,
                top = NULL,
                basal = NULL,
                id = time,
                time = NULL) %>%
        pivot_longer(
            cols = -c(id, model), 
            names_to = "stat",
            values_to = "stat_val")

df_ext$id <- ordered(df_ext$id, levels=c("pre", "during"))

ggplot(df_ext,
        aes(x = factor(`id`), 
            y = stat_val, 
            colour = model,
            group = model)) +
    geom_line() +
    geom_point() +
    facet_wrap(vars(stat),
                scales = 'free') +
    scale_size(guide = 'none') +
    theme_classic() +
    xlab("time") +
    ylab("value") +
    theme(panel.border = element_rect(colour = 'black',
                                      fill = "#ffffff00"),
            axis.ticks.x = element_blank())

ggsave("figures/extinction.png",
       width = 9000,
       height = 5000,
       units = "px",
       dpi = 600)
