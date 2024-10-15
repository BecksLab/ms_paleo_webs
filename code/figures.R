library(here)
library(readr)
library(scales)
library(tidyverse)

# set path to code sub dir
setwd(here("code"))

#### Structure ####

df <- list.files(path = "../data/processed/", pattern = ".csv", full.names = TRUE) %>% 
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

summary <-
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
  coord_cartesian(clip = "off") +
  theme(panel.border = element_rect(colour = 'black',
                                    fill = "#ffffff00"),
        axis.ticks.x = element_blank())

ggsave("../figures/summary.png",
       summary,
       width = 9000,
       height = 5000,
       units = "px",
       dpi = 600)

#### Extinctions ####

df_ext <- read_csv("../data/processed/extinctions/extinctions.csv") %>% 
  select(-c(id, links, richness)) %>% 
  # to get the ratio
  mutate(ratio = top/basal,
         top = NULL,
         basal = NULL,
         id = time,
         time = NULL) %>%
  pivot_longer(
    cols = -c(id, model, extinction_mechanism), 
    names_to = "stat",
    values_to = "end_val") %>% 
  filter(id == "pre") %>% 
  mutate(id = "post") %>% 
  left_join(.,
            df %>%
              filter(id == "pre") %>% 
              select(-id)) %>% 
  mutate(xstart = "pre",
         xend = id,
         start_val = stat_val,
         stat_val = NULL,
         id = NULL)

df_ext$xstart <- ordered(df_ext$xstart, levels = c("pre", "during", "post"))
df_ext$xend <- ordered(df_ext$xend, levels = c("pre", "during", "post"))

df_ext_summ <-
  df_ext %>%
  select(model, stat, end_val, start_val) %>% 
  pivot_longer(cols = c(end_val, start_val)) %>% 
  group_by(model, stat, name) %>% 
  summarise(y = mean(value)) %>%
  mutate(name = ifelse(name == "end_val",
                       "post",
                       "pre"))

summary +
  geom_line(data = df_ext_summ,
            aes(x = factor(name),
                y = y, 
                group = model),
            linetype = "dashed")

ggsave("../figures/extinction.png",
       width = 9000,
       height = 6000,
       units = "px",
       dpi = 600)

ggplot() +
  geom_segment(data = df_ext,
               aes(x = xstart,
                   y = start_val, 
                   xend = xend,
                   yend = end_val,
                   colour = model,
                   group = model,
                   linetype = extinction_mechanism),
               alpha = 0.3) +
  geom_line(data = df_ext_summ,
            aes(x = factor(name),
                y = y,
                colour = model,
                group = model)) +
  facet_wrap(vars(stat),
             scales = 'free') +
  scale_size(guide = 'none') +
  theme_classic() +
  xlab("time") +
  ylab("value") +
  coord_cartesian(clip = "off") +
  theme(panel.border = element_rect(colour = 'black',
                                    fill = "#ffffff00"),
        axis.ticks.x = element_blank())

ggsave("../figures/extinction_all_results.png",
       width = 9000,
       height = 5000,
       units = "px",
       dpi = 600)
