# Script for plotting example age distributions (supplemental figure S1 in )
# - 25 May 2023

library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)

rm(list = ls())

source('model_source/sim_model1_functions.R')

backbone = expand.grid(
  # (Equilibrium) size of initial cohort
  s.max = c(0.7, 0.8, 0.9),
  # Gamma squared (pheno variance / sel pressure)
  sig.z = sqrt((1:20)/50)
) %>%
  # Demographic rates
  mutate(
    # Maximum expected lifetime fitness
    w.max = 2,
    # Equilibrium lifetime fitness
    wstar = w.max * (1 - s.max) / (sqrt(1 + sig.z^2) - s.max),
    # Mean fecundity
    r     = w.max * (1 - s.max) / s.max,
    # Equilibrium population growth rate
    lstar = (s.max + w.max * (1 - s.max)) / (s.max + (w.max/wstar) * (1 - s.max)),
    # p0 in case it is useful
    p0    = (w.max * (1 - s.max)) / (w.max * (1 - s.max) + s.max)
  ) %>%
  # Genetic info
  filter(s.max < lstar) %>%
  group_by(lstar, s.max, p0) %>%
  mutate(
    # Gamma-parameterization
    # wfitn = 1 in gamma parameterization
    wfitn = 1,
    # Phenotypic standard deviation in new cohorts
    sig.0 = sqrt(newt.method.g1(ifelse(sig.z^2 < .1, .001, .1), 1e-8, s.max / lstar, r)),
    # Breeding value standard deviation in new cohorts
  ) 

for.plot = merge(
  backbone,
  expand.grid(
    # (Equilibrium) size of initial cohort
    s.max = c(0.7, 0.8, 0.9),
    # Gamma squared (pheno variance / sel pressure)
    sig.z = sqrt((1:40)/100),
    # Age
    k = 0:15
  )
) %>%
  mutate(p.k = p0 * (s.max/lstar)^k * sqrt(1 / (1 + k*sig.0^2)))

# Plot with probability of age on the natural scale
plot.naturale = for.plot %>%
  ggplot(aes(x = k, y = p.k)) +
  geom_line(aes(group = sig.z, colour = sig.z^2)) +
  scale_colour_gradient(low = 'darkviolet', high = 'darkorange', breaks = c(.1, .25, .4)) +
  facet_wrap(~ s.max, labeller = label_bquote(cols = hat(s) == .(s.max))) +
  labs(x = '', y = expression(p[age])) +
  theme(
    panel.background = element_blank(),
    legend.position = 'none'
  )

# Plot with probability of age on the log scale
plot.logscale = for.plot %>%
  ggplot(aes(x = k, y = p.k)) +
  geom_line(aes(group = sig.z, colour = sig.z^2)) +
  scale_colour_gradient(low = 'darkviolet', high = 'darkorange', breaks = c(.1, .25, .4)) +
  scale_y_log10() +
  facet_wrap(~ s.max, labeller = label_bquote(cols = hat(s) == .(s.max))) +
  labs(x = 'age', y = expression(p[age])) +
  theme(
    panel.background = element_blank(),
    strip.background = element_blank(),
    strip.text = element_blank(),
    legend.position = 'none'
  )

plot.legend = get_legend(
  plot.naturale +
    guides(
      colour = guide_legend(expression(phenotypic ~ variance ~ gamma^2)),
      override.aes = list(linewidth = 1.5)
    ) +
    theme(legend.position = 'top')
)

plot_grid(
  plot.legend,
  plot.naturale,
  plot.logscale,
  nrow = 3, 
  rel_heights = c(.2, 1, 1)
)

ggsave(filename = 'analyze_results/figs_out/figS1_age_distn.png',
       width = 8, height = 6)
