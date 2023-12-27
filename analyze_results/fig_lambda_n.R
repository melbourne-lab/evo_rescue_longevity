### Script for re-making Figure 1 in manuscript (two panels)
### - a) lambda as a function of z for for LH + variance groups
### - b) population size over time for LH groups

# Load in libraries

library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)

rm(list = ls())

# Plot for panel (a),

lambda.z = expand.grid(
  # Three different phenotypic standard deviation values
  gamma = sqrt(c(.1, .25, .4)),
  # Three different maximum survivals (defining life histories)
  s.max = c(.1, .5, .9)
) %>%
  mutate(
    # Get equilibrium population growth rates, lambda^*
    w.max = 3,
    l.max = s.max + (1 - s.max) * w.max,
    lstar = l.max / sqrt(1 + gamma^2)
  ) %>%
  # Get the lambda values as a function of z
  merge(
    y = expand.grid(
      lstar = .$lstar, #  %>% pull(lstar),
      z     = (0:300)/100
    ) %>%
      mutate(lamdz = exp(-z^2 / 2) * lstar)
  ) %>%
  # Some labels for plotting and faceting
  mutate(var.z = factor(round(gamma^2, 2)))

lambda.plot = lambda.z %>%
  mutate(s.max = factor(s.max, levels = c(.9, .5, .1))) %>%
  ggplot(aes(x = z, y = lamdz, group = interaction(s.max, gamma))) +
  geom_segment(
    aes(x = 0, xend = 3, y = 1, yend = 1),
    linetype = 5, colour = 'gray88'
  ) +
  geom_line(aes(colour = s.max, linetype = var.z)) +
  # x-axis ticks are z-tilde
  geom_rug(
    data = lambda.z %>%
      distinct(gamma, s.max, w.max, var.z) %>%
      mutate(ztilde = sqrt(2 * log(w.max * (1-s.max) + s.max) - log(1 + gamma^2))) %>%
      mutate(s.max = factor(s.max, levels = c(.9, .5, .1))),
    inherit.aes = FALSE,
    aes(x = ztilde, colour = s.max, linetype = var.z),
    sides = 'b', length = grid::unit(0.075, 'npc')
  ) +
  labs(
    x = expression(paste('Mean population phenotype ', bar(z))),
    y = expression(lambda)
  ) +
  scale_colour_manual(values = c("#999999", "#56B4E9", "#E69F00"), 'longevity') +
  # scale_colour_brewer(palette = 'Dark2', 'longevity') +
  scale_linetype_manual(values = c(1, 5, 2)) +
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
  mutate(
    n0 = log(20000),
    w.max = 3,
    l.max = s.max + w.max * (1 - s.max),
    gamma = sqrt(.1),
    k  = (1 + .5*gamma^2) / (1 + gamma^2),
    nt = n0 + 
      tstep * (log(l.max) - .5 * log(1 + gamma^2)) -
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

fig.hypo.legend = get_legend(
  lambda.plot +
    guides(
      linetype = guide_legend(expression(gamma^2)),
      colour   = guide_legend(expression(hat(s)))
    ) +
    theme(legend.position = 'top', legend.key = element_rect(fill = NA))
)

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
