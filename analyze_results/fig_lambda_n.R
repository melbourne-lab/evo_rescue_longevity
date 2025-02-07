### Script for re-making Figure 1 in manuscript (two panels)
### - a) lambda as a function of z for for LH + variance groups
### - b) population size over time for LH groups

# Load in libraries

library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)
library(cowplot)

rm(list = ls())

# Load in aux functions
source('model_source/sim_model1_functions.R')

lambdas = data.frame(s.max = c(.1, .5, .9)) %>%
  mutate(
    w.max = 3,
    r = w.max * (1 - s.max) / s.max,
    gamma = sqrt(0.1)
  )

lambdas = cbind(
  lambdas,
  lambdas %>%
    split(~ s.max) %>%
    map(\(x) newt.method.2d(1, .25, x$s.max, x$r, x$gamma^2, 1e-6)) %>%
    do.call(rbind, .)
) %>%
  rename(
    lstar = lam,
    sig.z = g2_0
  ) %>%
  mutate(sig.z = sqrt(sig.z))

# Plot for panel (a),

lambda.z = lambdas %>%
  merge(y = expand.grid(k = 0:500, z = (0:300)/100)) %>%
  mutate(
    pk = (r / (1+r)) * (s.max / lstar)^k * 1 / sqrt(1 + k*sig.z^2),
    sbark = s.max * sqrt((1 + k*sig.z^2) / (1 + (k+1)*sig.z^2)) * exp(-z^2 / 2)
  ) %>%
  group_by(s.max, r, z, sig.z, lstar) %>%
  summarise(lambz = sum(pk * sbark * (1 + r))) %>%
  ungroup()

lambda.plot = lambda.z %>%
  # mutate(s.max = factor(s.max, levels = c(.9, .5, .1))) %>%
  mutate(long = factor(s.max, labels = c('high', 'medium', 'low'), levels = c(0.9, 0.5, 0.1))) %>%
  ggplot(aes(x = z, y = lambz, group = long)) +
  geom_segment(
    aes(x = 0, xend = 3, y = 1, yend = 1),
    linetype = 5, colour = 'gray88'
  ) +
  # geom_line(aes(colour = s.max)) +
  geom_line(aes(colour = long)) +
  # x-axis ticks are z-tilde
  # geom_rug(
  #   data = lambda.z %>%
  #     distinct(gamma, s.max, w.max, var.z) %>%
  #     mutate(ztilde = sqrt(2 * log(w.max * (1-s.max) + s.max) - log(1 + gamma^2))) %>%
  #     mutate(s.max = factor(s.max, levels = c(.9, .5, .1))),
  #   inherit.aes = FALSE,
  #   aes(x = ztilde, colour = s.max, linetype = var.z),
  #   sides = 'b', length = grid::unit(0.075, 'npc')
  # ) +
  labs(
    x = expression(paste('Mean population phenotype ', bar(z))),
    y = expression(lambda)
  ) +
  scale_colour_manual(values = c("#999999", "#56B4E9", "#E69F00"), 'longevity') +
  # scale_colour_brewer(palette = 'Dark2', 'longevity') +
  # scale_linetype_manual(values = c(1, 5, 2)) +
  theme(
    panel.background = element_blank(),
    legend.position = 'none'
  )

lambda.plot

#####

# Plot for panel (b)

g.and.h = expand.grid(
  # Life history groups
  s.max = c(.1, .5, .9),
  # Time steps (x-axis)
  tstep = 0:30
) %>%
  # Get lambda* and gamma^2
  merge(lambdas) %>%
  mutate(
    n0 = log(20000),
    k  = (1 + .5*gamma^2) / (1 + gamma^2),
    nt = n0 + 
      tstep * (log(lstar)) -
      (2^2 / (2 * (1 + gamma^2))) * (1 - k^(2*tstep)) / (1 - k^2)
  )

n.plot = g.and.h %>%
  mutate(s.max = factor(s.max, levels = c(.9, .5, .1))) %>%
  filter(nt > 0, nt <= log(20000)) %>%
  ggplot(aes(x = tstep, y = exp(nt), group = s.max)) + 
  geom_line(aes(colour = s.max)) +
  scale_y_log10() +
  scale_colour_manual(values = c("#999999", "#56B4E9", "#E69F00"), 'longevity') +
  # scale_colour_brewer(palette = 'Dark2', 'longevity') +
  labs(x = expression(Time ~ step ~ t), y = expression(N[t])) +
  theme(
    panel.background = element_blank(),
    legend.position = 'none'
  )

n.plot
# (ignore the warnings - they are just about points getting chopped from the plot)

fig.hypo.legend = get_plot_component(
  lambda.plot +
    guides(colour = guide_legend(linetype = 'none')) +
    theme(legend.position = 'top'),
  'guide-box', return_all = TRUE
)[[4]]

plot_grid(
  fig.hypo.legend,
  plot_grid(
    lambda.plot,
    n.plot,
    nrow = 1,
    labels = c('a)', 'b)')
  ),
  rel_heights = c(.1, 1),
  nrow = 2
)

ggsave(filename = 'analyze_results/figs_out/fig1_lambda_n.png',
       width = 8, height = 3)
