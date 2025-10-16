# Equal-lambda plot that shows first four generations (all pre-extinction):
# - Population size (top row)
# - Coefficient of variation in 
# - Phenotypic components
# Includes heritability levels (across which rates of adaptation change, due to
# decoupling), all with the same set of phenotypic variance (high)
# This produces Figure S22

##### Setup

library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)
library(ggborderline) # for the line borders

# Clear namespace
rm(list = ls())


##### Read in data

# As it happens, I 

# Read in mean population size/growth rate data (all population sizes)
n.all = read.csv('run_sims/out/equal_lambda/sim_results_m1_ages.csv') %>%
  # Get n (file here has size of each age group/cohort)
  group_by(t, trial, h2, lstar) %>%
  summarise(n = sum(n)) %>%
  ungroup() %>%
  # Factor for plotting aesthetics
  mutate(
    # (mistake somehow in original code... s.max not recorded; we'll extract this from lstar)
    long = cut(lstar, breaks = c(0, 1.185, 1.19, Inf), labels = c('low', 'medium', 'high')),
    # (relevel to make sure the high longevity is reference)
    long = factor(long, levels = c('high', 'medium', 'low')),
    hert = factor(paste0('h^2 == ', h2))
  )

# Read in individual trial-level data on population growth
n.ind = n.all %>% filter(((trial - 1) %% 200) <= 20)
  
# Coefficient of variation for lambda
# n.cvl is n (population size) - coefficient variation for lambda
n.cvl = n.all %>%
  # Give us only the first five time steps (we want four lambdas)
  filter(t < 6) %>%
  # Estimate lambda
  group_by(trial, long, hert) %>%
  mutate(lambda = c(exp(diff(log(n))), NA)) %>%
  # get rid of terminating NA
  filter(!is.na(lambda)) %>%
  group_by(t, long, hert) %>%
  summarise(cv.lambda = sd(lambda) / mean(lambda)) %>%
  ungroup()

# Summarize (mean and standard error) population size in each time step
n.sum = n.all %>%
  # Adding in zeros (not recorded when simulation ended but useful for visualizing results in figure)
  complete(t, nesting(trial, long, hert)) %>%
  mutate(n = ifelse(is.na(n), 0, n)) %>% 
  group_by(t, long, hert) %>%
  summarise(
    nbar = mean(n), 
    nvar = var(n),
    nn = n()
  ) %>%
  ungroup()

# Read in and get means for phenotypic component data
# Get means per group
z.bar = read.csv('run_sims/out/equal_lambda/sim_results_m1_phenotypic_components.csv') %>%
  group_by(t, s.max, h2) %>%
  summarise(
    across(
      contains('bar'),
      list(bar = function(x) mean(x, na.rm = TRUE), var = function(x) var(x, na.rm = TRUE)),
      .names = "{.col}.{.fn}"
    ),
    nn = n()
  ) %>%
  pivot_longer(
    cols = contains('ar'),
    names_to = 'vartype',
    values_to = 'varval'
  ) %>%
  separate(col = 'vartype', into = c('comp', 'stat'), sep = '\\.') %>%
  pivot_wider(
    names_from = stat,
    values_from = varval
  )


##### Make two-panel plot using cowplot

pan.a = n.sum %>%
  filter(t < 5) %>%
  ggplot(aes(x = t, group = long)) +
  geom_line(
    data = n.ind %>% filter(t < 5),
    aes(
      y = n,
      colour = long,
      group = trial
    ),
    linewidth = 0.125
  ) +
  geom_borderline(
    aes(
      y = nbar,
      colour = long
    ),
    bordercolour = 'black', linewidth = 1.2, borderwidth = 0.125
  ) +
  geom_ribbon(
    aes(
      ymin = nbar - 2 * sqrt(nvar / nn),
      ymax = nbar + 2 * sqrt(nvar / nn),
      fill = long
    ),
    alpha = 0.2
  ) +
  labs(x = '', y = 'Population size') +
  scale_colour_manual(values = c("#999999", "#56B4E9", "#E69F00"), 'longevity') +
  scale_fill_manual(values = c("#999999", "#56B4E9", "#E69F00"), 'longevity') +
  # scale_colour_brewer(palette = 'Dark2', 'longevity') +
  # scale_fill_brewer(palette = 'Dark2', 'longevity') +
  scale_y_log10(breaks = c(1000, 10000)) +
  facet_grid(cols = vars(hert), labeller = label_parsed) +
  theme(
    panel.background = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    # panel.grid.major = element_line(colour = 'gray88', linewidth = 0.25),
    legend.position = 'top',
    text = element_text(family = 'ArialMT')
    # used for visuals when fidgeting with cowplot margins only
    # panel.border = element_rect(fill = NA),
  )

# Panel with coefficient of variation for lambda
pan.b = n.cvl %>%
  ggplot(aes(x = t, y = cv.lambda)) + 
  geom_line(aes(group = long, colour = long), linewidth = 1.1) +
  scale_colour_manual(values = c("#999999", "#56B4E9", "#E69F00"), 'longevity') +
  facet_wrap(~ hert, labeller = label_parsed) +
  labs(x = '', y = expression('CV' ~ lambda)) +
  theme(
    panel.background = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y.left = element_text(vjust = 5),
    strip.text = element_blank(),
    legend.position = 'none',
    # used for visuals when fidgeting with cowplot margins only
    # panel.border = element_rect(fill = NA),
    plot.margin = margin(l = 17),
    text = element_text(family = 'ArialMT')
  )

pan.c = z.bar %>% 
  # For only extant populations
  filter(t < 5) %>%
  ungroup() %>%
  # Factor for plotting aesthetics
  mutate(
    long = factor(s.max, labels = c('high', 'medium', 'low'), levels = c(0.9, 0.5, 0.1)),
    # long = relevel(long, ref = 'high'),
    # hert = factor(paste0('heritability ', h2)),
    hert = factor(paste0("h^2 == ", h2)),
    comp = relevel(factor(comp), ref = 'zbar')
  ) %>%
  # Re-arrange to get the low longevity on top
  arrange(desc(long)) %>%
  ggplot(aes(x = t)) +
  geom_ribbon(
    aes(
      ymin = bar - 2 * sqrt(var / nn),
      ymax = bar + 2 * sqrt(var / nn),
      group = interaction(comp, long),
      fill = long
    ),
    alpha = 0.2
  ) +
  geom_line(
    aes(
      y = bar,
      group = interaction(comp, long),
      colour = long,
      linetype = comp
    ),
    linewidth = 1.1
  ) +
  scale_colour_manual(values = c("#999999", "#56B4E9", "#E69F00"), 'longevity') +
  scale_fill_manual(values = c("#999999", "#56B4E9", "#E69F00"), 'longevity') +
  # scale_colour_manual(values = c("#E69F00", "#56B4E9", "#999999"), 'longevity') +
  # scale_fill_manual(values = c("#E69F00", "#56B4E9", "#999999"), 'longevity') +
  # scale_colour_brewer(palette = 'Dark2', 'longevity') +
  # scale_fill_brewer(palette = 'Dark2', 'longevity') +
  scale_linetype(
    'component', 
    labels = c(expression(bar(z)), expression(bar(b)), expression(bar(e)))
  ) +
  scale_y_continuous(breaks = 0:2) +
  labs(x = 'Time step', y = 'Value') +
  guides(
    colour = 'none',
    fill   = 'none',
    linetype = guide_legend(order = 1, override.aes = list(linewidth = 1))
  ) +
  facet_wrap(~ hert, labeller = labeller(hert = label_parsed)) +
  theme(
    panel.background = element_blank(), 
    strip.text = element_blank(), 
    legend.position = 'top',
    # used for visuals when fidgeting with cowplot margins only
    # panel.border = element_rect(fill = NA),
    axis.title.y.left = element_text(vjust = 8),
    plot.margin = margin(l = 25),
    text = element_text(family = 'ArialMT')
  )

# Combine into single figure and export (Fig. S22)
plot_grid(
  pan.a, pan.b, pan.c, nrow = 3,
  rel_heights = c(1.1, .65, 1),
  labels = c('a)', 'b)', 'c)')
) %>%
  save_plot(
    filename = 'analyze_results/figs_out/figS22_combined_popsize_cvlambda_components.pdf',
    base_width = 5, base_height = 5
  )


