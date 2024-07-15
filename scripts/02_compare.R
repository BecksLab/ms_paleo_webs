##Load other packages that are needed
library(tidyverse)

# functions needed to infer food webs
source("lib/internals.R")

## combine summaries for different models
networks <- 
  tribble(
    ~model, ~links, ~connectance,
    "pfim", sum(pfim_web), connectance(pfim_web),
    "niche", sum(niche_web), connectance(niche_web)) %>%
  pivot_longer(
    cols = !model,
    names_to = "measure",
    values_to = "value"
  )

## plot 

ggplot(networks,
       aes(x = model,
           y = value,
           colour = model)) +
  geom_point() +
  facet_wrap(vars(measure),
             scales = 'free') +
  theme_classic()