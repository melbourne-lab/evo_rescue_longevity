# Script for plotting example age distributions (supplemental figure S2)

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
  sig.z = sqrt((0:25)/50)
) %>%
  # Demographic rates
  mutate(
    # Maximum expected lifetime fitness
    w.max = 3,
    # Mean fecundity
    r     = w.max * (1 - s.max) / s.max,
  )

backbone = cbind(
    backbone,
    backbone %>%
      split(~ s.max + sig.z) %>%
      map(\(x) newt.method.2d(1.1, .4, x$s.max, x$r, x$sig.z^2, 1e-5)) %>%
      do.call(rbind, .)
  )

for.plot = merge(
  backbone %>% rename(lstar = lam) %>% mutate(sig.0 = sqrt(g2_0)), 
  data.frame(k = 0:15)
  # expand.grid(
  #   # (Equilibrium) size of initial cohort
  #   s.max = c(0.7, 0.8, 0.9),
  #   # Gamma squared (pheno variance / sel pressure)
  #   sig.z = sqrt((1:40)/100),
  #   # Age
  #   k = 0:15
  # )
) %>%
  mutate(p.k = (r / (1 + r)) * (s.max/lstar)^k * sqrt(1 / (1 + k*sig.0^2)))

# Plot with probability of age on the natural scale
plot.naturale = for.plot %>%
  ggplot(aes(x = k, y = p.k)) +
  geom_line(aes(group = sig.z, colour = sig.z^2)) +
  scale_colour_gradient(low = 'darkviolet', high = 'darkorange', breaks = c(0, .25, .5)) +
  facet_wrap(~ s.max, labeller = label_bquote(cols = hat(s) == .(s.max))) +
  labs(x = '', y = expression(p[age])) +
  theme(
    panel.background = element_blank(),
    legend.position = 'none',
    text = element_text(family = 'ArialMT')
  )

# Plot with probability of age on the log scale
plot.logscale = for.plot %>%
  ggplot(aes(x = k, y = p.k)) +
  geom_line(aes(group = sig.z, colour = sig.z^2)) +
  scale_colour_gradient(low = 'darkviolet', high = 'darkorange', breaks = c(0, .25, .5)) +
  scale_y_log10() +
  facet_wrap(~ s.max, labeller = label_bquote(cols = hat(s) == .(s.max))) +
  labs(x = 'age', y = expression(p[age])) +
  theme(
    panel.background = element_blank(),
    strip.background = element_blank(),
    strip.text = element_blank(),
    legend.position = 'none',
    text = element_text(family = 'ArialMT')
  )

plot.legend = get_plot_component(
  plot.naturale +
    guides(
      colour = guide_coloursteps(expression(phenotypic ~ variance ~ gamma^2)),
      override.aes = list(linewidth = 3)
    ) +
    theme(
      legend.position = 'top', 
      legend.text = element_text(hjust = -0.25),
      legend.title = element_text(vjust = 0.75)
    ),
  pattern = "guide-box", return_all = TRUE
)[[4]]

plot_grid(
  plot.legend,
  plot.naturale,
  plot.logscale,
  nrow = 3, 
  rel_heights = c(.2, 1, 1)
)

ggsave(filename = 'analyze_results/figs_out/figS2_age_distn.pdf',
       width = 8, height = 5)
