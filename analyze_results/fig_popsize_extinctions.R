### Script for reproducing figs. 2 and 3 (population size and extinction rates
### over time) from the large batch of simulations

library(ggplot2)
library(dplyr)
library(tidyr)

rm(list = ls())

# Read in source model functions just in case
source('model_source/sim_model1_functions.R')

# Number of trials per group
trys.per = 1000

# Read in data (all population sizes)
n.all = read.csv('run_sims/out/sim_results_m1_allsizes_n.csv') %>%
  # Factor for plotting aesthetics
  mutate(
    long = factor(s.max, labels = c('low', 'medium', 'high'), levels = c(0.9, 0.5, 0.1)),
    hert = factor(paste0('h^2 == ', h2)),
    varn = factor(paste0('gamma^2 == ', var.z))
  )

nrow(n.all)
head(n.all)

# Read in data for individual trials
n.ind = read.csv('run_sims/out/sim_results_m1_disaggregated_n.csv') %>% 
  # Factor for plotting aesthetics
  mutate(
    long = factor(s.max, labels = c('low', 'medium', 'high'), levels = c(0.9, 0.5, 0.1)),
    hert = factor(paste0('h^2 == ', h2)),
    varn = factor(paste0('gamma^2 == ', var.z))
  )

nrow(n.ind)
head(n.ind)

# Make plot of mean and example population sizes
n.all %>%
  ggplot(aes(x = t, group = long)) +
  geom_line(
    aes(
      y = nbar,
      colour = long
    ),
    linewidth = 1.2
  ) +
  geom_ribbon(
    aes(
      ymin = nbar - 2 * sqrt(nvar / nn),
      ymax = nbar + 2 * sqrt(nvar / nn),
      fill = long
    ),
    alpha = 0.2
  ) +
  geom_line(
    data = n.ind,
    aes(
      y = n,
      colour = long,
      group = trial
    ),
    linewidth = 0.125
  ) +
  labs(x = 'Time step', y = 'Population size') +
  scale_colour_manual(values = c("#999999", "#56B4E9", "#E69F00"), 'longevity') +
  scale_fill_manual(values = c("#999999", "#56B4E9", "#E69F00"), 'longevity') +
  # scale_colour_manual(values = c("#E69F00", "#56B4E9", "#999999"), 'longevity') +
  # scale_fill_manual(values = c("#E69F00", "#56B4E9", "#999999"), 'longevity') +
  # scale_colour_brewer(palette = 'Dark2', 'longevity') +
  # scale_fill_brewer(palette = 'Dark2', 'longevity') +
  facet_grid(rows = vars(varn), cols = vars(hert), labeller = label_parsed) +
  theme(
    panel.background = element_blank(),
    legend.position = 'top'
  )

ggsave('analyze_results/figs_out/fig2_popsize.png',
       width = 5, height = 5)

# Plot of extinctions over time
n.all %>%
  ggplot(aes(x = t, y = psrv)) +
  geom_line(aes(colour = long), linewidth = 1.2) +
  labs(x = 'Time step', y = 'Proportion surviving') +
  scale_colour_manual(values = c("#999999", "#56B4E9", "#E69F00"), 'longevity') +
  # scale_colour_manual(values = c("#E69F00", "#56B4E9", "#999999"), 'longevity')
  # scale_colour_brewer(palette = 'Dark2', 'longevity') +
  facet_grid(rows = vars(varn), cols = vars(hert), labeller = label_parsed) +
  theme(
    panel.background = element_blank(),
    legend.position = 'top',
  )

ggsave('analyze_results/figs_out/fig3_extinctions.png',
       width = 5, height = 5)
