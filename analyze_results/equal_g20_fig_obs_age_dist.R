### Script for generating figures with observed age distributions
### This is for the simulations with equal initial cohort variance

library(ggplot2) # 3.4.3
library(dplyr)   # 1.1.3
library(tidyr)   # 1.3.0

rm(list = ls())

v.age = read.csv('run_sims/out/sim_results_m1_ages_equalg20.csv')

nrow(v.age) # 0.9 million rows
head(v.age) # that's a lot of info

# Number of trials per parameter combination here is 200 (small batch)
trys.per = 200

# Get only extant time steps for each treatment
v.age.extant = v.age %>%
  group_by(s.max, h2, t) %>%
  mutate(n.trials = length(unique(trial))) %>%
  group_by(s.max, h2) %>%
  mutate(max.t = max(t[n.trials == trys.per])) %>%
  ungroup() %>%
  filter(t <= max.t)

# Summarise
age.dist = v.age.extant %>%
  group_by(trial, t) %>%
  mutate(n.total = sum(n)) %>%
  group_by(s.max, h2, t, age) %>%
  summarise(p.age = sum(n / n.total) / trys.per) %>%
  group_by(s.max, h2, t) %>%
  mutate(
    p.age1 = cumsum(p.age),
    p.age0 = c(0, p.age1[-n()])
  ) %>%
  ungroup() %>%
  mutate(
    age2 = ifelse(age > 8, '9+', as.character(age)),
    age3 = ifelse(age > 3, '4+', as.character(age))
  )

# Nine panel plot (all data, all longevities)

age.dist %>%
  mutate(
    # long = longevity
    long = factor(s.max, labels = c('low', 'medium', 'high'), levels = c(0.9, 0.5, 0.1)),
    # hert = heritability
    hert = factor(paste("h^2 ==", h2))
  ) %>%
  ggplot(aes(x = t)) +
  geom_ribbon(
    aes(
      ymin = p.age0,
      ymax = p.age1,
      group = age,
      fill = age2
    ),
    colour = 'gray22',
    linewidth = 0.25
  ) +
  labs(x = 'Time step', y = 'Proportion of population in age group') +
  scale_fill_viridis_d(option = 'B', 'age') +
  facet_grid(rows = vars(long), cols = vars(hert), labeller = label_parsed) +
  theme(
    panel.background = element_rect(fill = 'gray99'),
    panel.grid = element_blank()
  )


### Plot for *only* high longevity

age.dist %>%
  filter(s.max > 0.5) %>%
  # hert = heritability (for facet labeling)
  mutate(hert = factor(paste("h^2 ==", h2))) %>%
  ggplot(aes(x = t)) +
  geom_ribbon(
    aes(
      ymin = p.age0,
      ymax = p.age1,
      group = age,
      fill = age2
    ),
    colour = 'gray22',
    linewidth = 0.25
  ) +
  labs(x = 'Time step', y = 'Proportion of population in age group') +
  scale_fill_viridis_d(option = 'B', 'age') +
  facet_grid(cols = vars(hert), labeller = label_parsed) +
  theme(
    panel.background = element_rect(fill = 'gray99'),
    panel.grid = element_blank()
  )

ggsave('analyze_results/figs_out/fig_supp_agedist_highlong_eg20.png', width = 8, height = 3)

# ### Plot for medium and low longevity groups
# 
# age.dist %>%
#   filter(s.max < 0.9) %>%
#   mutate(
#     # long = longevity
#     long = factor(s.max, labels = c('low longevity', 'medium longevity')),
#     # hert = heritability
#     hert = factor(paste("h^2 ==", h2))
#   ) %>%
#   mutate(hert = factor(paste("h^2 ==", h2))) %>%
#   ggplot(aes(x = t)) +
#   geom_ribbon(
#     aes(
#       ymin = p.age0,
#       ymax = p.age1,
#       group = age,
#       fill = age3
#     ),
#     colour = 'gray22',
#     linewidth = 0.25
#   ) +
#   labs(x = 'Time step', y = 'Proportion of population in age group') +
#   scale_fill_viridis_d(option = 'B', 'age') +
#   facet_grid(
#     rows = vars(long), cols = vars(hert), 
#     labeller = labeller(hert = label_parsed, long = label_value)
#   ) +
#   theme(
#     panel.background = element_blank()
#   )
# 
# ggsave('analyze_results/figs_out/figS20_agedist_medlowlong.png', width = 8, height = 5)

