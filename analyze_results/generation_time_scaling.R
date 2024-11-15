### Script for reproducing figs. 2 and 3 (population size and extinction rates
### over time) from the large batch of simulations

library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)

rm(list = ls())

# Read in source model functions just in case
source('model_source/sim_model1_functions.R')

pars = data.frame(
  # Longevity groups defined by s.max
  s.max = c(0.1, 0.5, 0.9)
) %>%
  mutate(
    # Max lifetime fitness
    w.max = 3,
    # Genetic fitness groups defined by 
    sig.z = sqrt(0.4),
    # Fecundity per time step
    r = w.max * (1 - s.max) / s.max
  ) 

# Estimate lambda* and gamma^2_0 (equilibrium pop growth rate and phenotypic
# variance in newborn cohorts, respectively)
# This is done numerically with a method for Newton's method
pars = cbind(
  pars,
  pars %>%
    split(~ s.max + sig.z) %>%
    map(\(x) newt.method.2d(1.1, .5, x$s.max, x$r, x$sig.z^2, 1e-5)) %>%
    do.call(rbind, .)
)

# Add remaining parameters
pars = pars %>%
  # Rename variables as needed
  rename(
    lstar = lam, 
    sig.0 = g2_0
  ) %>%
  # Wrapper function above gives the squared phenotypic variance; convert to
  # standard dev.
  mutate(sig.0 = sqrt(sig.0)) %>%
  # Repicate data frame above to include three heritability levels
  merge(data.frame(h2 = c(0.25, .5, 1.0))) %>%
  # Estimate genetic variance and environmental component variance in newborn
  # cohorts
  mutate(
    sig.a = sqrt(h2 * sig.0^2), 
    sig.e = sqrt((1 - h2) * sig.0^2)
  ) %>%
  # Add other necessary parameters as needed
  mutate(
    # Length of simulations
    timesteps = 100,
    # Initial pouplation size
    n.pop0 = 20000,
    # Ceiling-type carrying capacity
    kceil = 20000,
    # Initial population distance from phenotypic optimum (i.e., magnitude of
    # environmental shift)
    gbar0 = 2,
    # Width of fitness landscape
    # (setting this to 1 makes the simulation follow the gamma-parameterization
    # in manuscript)
    wfitn = 1,
    # Per-individual mutation rate
    mu    = 1,
    # Population-wide additive genetic variance (needed to get mutational
    # variance)
    sig.p = sqrt(gamma.a.calc(sig.a^2, s.max / lstar, r, sig.e^2)),
    # Mutational variance
    sig.m = sqrt(wfitn^2 * (sig.a^2 - sig.p^2) * (1 + r)),
  )

# Old pars setup below:
# pars = expand.grid(
#   # (Equilibrium) size of initial cohort
#   s.max = c(0.1, 0.5, 0.9),
#   # Heritability of fitness
#   h2    = c(.25, .5, 1)
# ) %>%
#   # Demographic rates
#   mutate(
#     # Maximum expected lifetime fitness
#     w.max = 3,
#     # Gamma squared (pheno variance / sel pressure)
#     sig.z = sqrt(.4),
#     # Equilibrium lifetime fitness
#     wstar = w.max * (1 - s.max) / (sqrt(1 + sig.z^2) - s.max),
#     # Mean fecundity
#     r     = w.max * (1 - s.max) / s.max,
#     # Equilibrium population growth rate
#     lstar = (s.max + w.max * (1 - s.max)) / (s.max + (w.max/wstar) * (1 - s.max)),
#     # Initial population size
#     n.pop0 = 20000,
#     # Strength of density dependence
#     alpha = log(lstar) / n.pop0,
#     # Ceciling-type carrying capacity just in case
#     kceil = 20000,
#     p0    = (w.max * (1 - s.max)) / (w.max * (1 - s.max) + s.max)
#   ) %>%
#   # Genetic info
#   group_by(lstar, s.max, h2, p0) %>%
#   mutate(
#     # Gamma-parameterization
#     # wfitn = 1 in gamma parameterization
#     wfitn = 1,
#     # Phenotypic standard deviation in new cohorts
#     sig.0 = sqrt(newt.method.g1(.1, 1e-8, s.max / lstar, r)),
#     # Breeding value standard deviation in new cohorts
#     sig.a = sqrt(h2 * sig.0^2),
#     # Non-inherited standard dxeviation in new cohorts
#     sig.e = sqrt((1-h2) * sig.0^2),
#     # Population-wide breeding value standard deviation
#     sig.p = sqrt(gamma.a.calc(sig.a^2, s.max / lstar, r, sig.e^2)),
#     mu    = 1,
#     sig.m = sqrt(wfitn^2 * (sig.a^2 - (sig.p^2 - p0*sig.a^2)/(1-p0))),
#     gbar0 = 2
#   ) %>%
#   ungroup() %>%
#   # Other junk
#   mutate(timesteps = 50)

# Generatio time is:
# T = log(R_0) / R
#   where R_0 is net reproductive number, R is intrinsic population growth rate
#   R = log(lambda*)
#   R_0 is just r (per-mating bout fecundity) times mean life expectancy
#   (unfortunately I don't have a neat, closed form expression for mean life
#   expectancy)

# pars %>%
#   distinct(s.max, r, lstar, sig.0) %>%
#   uncount(weight = 501) %>%
#   group_by(s.max, r, lstar, sig.0) %>%
#   mutate(
#     k = 0:500,
#     p_k = (r/(1+r)) * (s.max / lstar)^k * (1 / sqrt(1 + k*sig.0^2)),
#     mm = cumsum(k * p_k)
#   ) %>%
#   mutate(s.max = factor(s.max)) %>%
#   ggplot(aes(x = k, y = mm, group = s.max)) +
#   geom_line(aes(colour = s.max)) +
#   facet_wrap(~ sig.0^2)

gen.ts = pars %>%
  distinct(s.max, r, lstar, sig.0) %>%
  uncount(weight = 501) %>%
  group_by(s.max, r, lstar, sig.0) %>%
  mutate(
    k = 0:500,
    p_k = (r/(1+r)) * (s.max / lstar)^k * (1 / sqrt(1 + k*sig.0^2)),
    l_k = s.max^k * sqrt(1 / (1 + k*sig.0^2))
  ) %>%
  summarise(
    mean.age = sum(k * p_k),
    p.o.gaps = sum((k / (lstar)^k) * l_k)
  ) %>%
  ungroup() %>%
  mutate(
    Tgap = r * p.o.gaps,
    R_0 = r * mean.age,
    Tgen = log(R_0) / log(lstar)
  )
# generation time of .0001 is kinda silly!

gen.ts

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

gen.adjusted.z.bar = merge(
  x = z.bar,
  y = gen.ts %>% select(-c(r, lstar, sig.0)),
) %>%
  mutate(gen = t / Tgap)

# Plot of components over time
gen.adjusted.z.bar %>% 
  # For only extant populations
  filter(nn == trys.per) %>%
  ungroup() %>%
  # Factor for plotting aesthetics
  mutate(
    long = factor(s.max, labels = c('low', 'medium', 'high'), levels = c(0.9, 0.5, 0.1)),
    # hert = factor(paste0('heritability ', h2)),
    hert = factor(paste0("h^2 == ", h2)),
    comp = relevel(factor(comp), ref = 'zbar')
  ) %>%
  filter(gen < 11) %>%
  ggplot(aes(x = gen)) +
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
  # geom_segment(
  #   data = data.frame(x = 1:9, xend = 1:9),
  #   aes(x = x, xend = xend, y = 0, yend = 2),
  #   colour = 'gray55', linewidth = 0.1
  # ) +
  # geom_point(
  #   aes(
  #     y = bar,
  #     colour = long
  #   ),
  #   size = 5, pch = 4
  # ) +
  scale_colour_manual(values = c("#999999", "#56B4E9", "#E69F00"), 'longevity') +
  scale_fill_manual(values = c("#999999", "#56B4E9", "#E69F00"), 'longevity') +
  # scale_colour_manual(values = c("#E69F00", "#56B4E9", "#999999"), 'longevity') +
  # scale_fill_manual(values = c("#E69F00", "#56B4E9", "#999999"), 'longevity') +  # scale_colour_brewer(palette = 'Dark2', 'longevity') +
  # scale_fill_brewer(palette = 'Dark2', 'longevity') +
  scale_linetype(
    'component', 
    labels = c(expression(bar(z)), expression(bar(b)), expression(bar(e)))
  ) +
  scale_x_continuous(breaks = c(0:5)*2) +
  labs(x = 'Generation', y = 'Value') +
  guides(
    colour = guide_legend(order = 2),
    fill   = guide_legend(order = 2),
    linetype = guide_legend(order = 1, override.aes = list(linewidth = 1))
  ) +
  facet_wrap(~ hert, labeller = labeller(hert = label_parsed)) +
  theme(panel.background = element_blank(), legend.position = 'top')
  
ggsave('analyze_results/figs_out/figS19_generation_time_scaling.png',
       width = 8, height = 5)
