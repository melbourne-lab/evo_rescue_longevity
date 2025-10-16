### Code for re-creating figures of phenotypic component means and variances
### from the shorter batch (200 sims only) over time.
### (In the MS these are Figs. 4 and 6)

# Setup
library(ggplot2)
library(dplyr)
library(tidyr)

rm(list = ls())

# Load in data
z.com = read.csv('run_sims/out/sim_results_m1_phtype.csv')

# Assess
nrow(z.com)
head(z.com)

# Trials per combination here is 200
trys.per = 200

# Get means per group
z.bar = z.com %>%
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

head(z.bar)

# Get variances per group
z.var = z.com %>%
  group_by(t, s.max, h2) %>%
  summarise(
    across(
      contains('var'),
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

head(z.var)

# Plot of component means over time (Fig. 4)
z.bar %>% 
  # For only extant populations
  filter(nn == trys.per) %>%
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
  labs(x = 'Time step', y = 'Value') +
  guides(
    colour = guide_legend(order = 2),
    fill   = guide_legend(order = 2),
    linetype = guide_legend(order = 1, override.aes = list(linewidth = 1))
  ) +
  facet_wrap(~ hert, labeller = labeller(hert = label_parsed)) +
  theme(
    panel.background = element_blank(), 
    legend.position = 'top',
    text = element_text(family = 'ArialMT')
  )

ggsave('analyze_results/figs_out/fig4_pheno_comp_mean.pdf',
       width = 8, height = 3)


# Phenotypic component variance plot (Fig. 6)
z.var %>% 
  # Get only extant populations
  # get rid of 'w' vars - i.e., variance components after selection
  filter(nn == trys.per, !grepl('w', comp)) %>%
  ungroup() %>%
  # Factor for plotting aesthetics
  mutate(
    long = factor(s.max, labels = c('high', 'medium', 'low'), levels = c(0.9, 0.5, 0.1)),
    # long = relevel(long, ref = 'high'),
    # hert = factor(paste0('heritability ', h2)),
    hert = factor(paste0("h^2 == ", h2)),
    comp = relevel(factor(comp), ref = 'zvar')
  ) %>%
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
  scale_linetype(
    'component', 
    labels = c(expression(gamma^2), expression(gamma[a]^2), expression(gamma[e]^2))
  ) +
  scale_colour_manual(values = c("#999999", "#56B4E9", "#E69F00"), 'longevity') +
  scale_fill_manual(values = c("#999999", "#56B4E9", "#E69F00"), 'longevity') +
  # scale_colour_manual(values = c("#E69F00", "#56B4E9", "#999999"), 'longevity') +
  # scale_fill_manual(values = c("#E69F00", "#56B4E9", "#999999"), 'longevity') +
  # scale_colour_brewer(palette = 'Dark2', 'longevity') +
  # scale_fill_brewer(palette = 'Dark2', 'longevity') +  
  labs(x = 'Time step', y = 'Value') +
  guides(linetype = guide_legend(order = 2, override.aes = list(linewidth = 1))) +
  facet_wrap(~ hert, labeller = labeller(hert = label_parsed)) +
  theme(
    panel.background = element_blank(), 
    legend.position = 'top',
    text = element_text(family = 'ArialMT')
  )

ggsave('analyze_results/figs_out/fig6_pheno_comp_vars.pdf',
       width = 8, height = 3)
