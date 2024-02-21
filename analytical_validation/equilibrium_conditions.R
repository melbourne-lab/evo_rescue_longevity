##########
# Script to validate analytical expressions for steady state
# 
# - Simulations were run at steady state (environmental optimum of zero,
#   population matched to optimum)
# - Parameter combinations are same as in simulated experiment in main text
# - Each trial is run for 10 time steps
# - I aggregate these to get means for each time step; in some instances I
#   also aggregate by age to get one estimate of the mean per time step per age
#   (because this is at steady state, the time steps should all be the same on average)
# - Script contains plots for: 
#   - Population size over time
#   - Age structure - probability of being in each age
#   - Mean of each phenotypic component (population-level)
#   - Variance of each phenotypic component (per age)
# SN - init 11 Dec. 2023
##########

# Packages
library(ggplot2)
library(dplyr)
library(tidyr)
library(parallel)

# Clear namespace
rm(list = ls())

# Load source code
source('model_source/sim_model1_functions.R')

### Load in parameters

# Trials per parameter combo
trys.per = 25

# Parameters

pars = expand.grid(
  # (Equilibrium) size of initial cohort
  s.max = c(0.1, 0.5, 0.9),
  # Heritability of fitness
  h2    = c(.25, .5, 1),
  # Gamma squared (pheno variance / sel pressure)
  sig.z = sqrt(c(.1, .25, .4))
) %>%
  # Demographic rates
  mutate(
    # Maximum expected lifetime fitness
    w.max = 3,
    # Equilibrium lifetime fitness
    wstar = w.max * (1 - s.max) / (sqrt(1 + sig.z^2) - s.max),
    # Mean fecundity
    r     = w.max * (1 - s.max) / s.max,
    # Equilibrium population growth rate
    lstar = (s.max + w.max * (1 - s.max)) / (s.max + (w.max/wstar) * (1 - s.max)),
    # Initial population size
    n.pop0 = 10000,
    # Strength of density dependence
    alpha = log(lstar) / n.pop0,
    # Ceciling-type carrying capacity just in case
    kceil = 20000,
    p0    = (w.max * (1 - s.max)) / (w.max * (1 - s.max) + s.max)
  ) %>%
  # Genetic info
  group_by(lstar, s.max, h2, p0) %>%
  mutate(
    # Gamma-parameterization
    # wfitn = 1 in gamma parameterization
    wfitn = 1,
    # Phenotypic standard deviation in new cohorts
    sig.0 = sqrt(newt.method.g1(.1, 1e-8, s.max / lstar, r)),
    # Breeding value standard deviation in new cohorts
    sig.a = sqrt(h2 * sig.0^2),
    # Non-inherited standard deviation in new cohorts
    sig.e = sqrt((1-h2) * sig.0^2),
    # Population-wide breeding value standard deviation
    sig.p = sqrt(gamma.a.calc(sig.a^2, s.max / lstar, r, sig.e^2)),
    mu    = 1,
    sig.m = sqrt(wfitn^2 * (sig.a^2 - (sig.p^2 - p0*sig.a^2)/(1-p0))),
    gbar0 = 0
  ) %>%
  ungroup() %>%
  # Other junk
  mutate(timesteps = 10)

# Run simulations
set.seed(641876)

sim.out = mclapply(
  pars %>% uncount(trys.per) %>% mutate(try.no = 1:(nrow(.))) %>% split(.$try.no),
  function(pars) {
    sim(pars, theta.t = 0, init.rows = 30 * 30000) %>%
      mutate(e_i = z_i - b_i) %>%
      group_by(t, age) %>%
      summarise(
        n = n(),
        across(
          c(b_i, e_i, z_i), 
          .fns = list(bar = mean, var = ~ ifelse(n == 1, 0, var(.x) * (n/(n-1))))
        ),
        rho = ifelse(n == 1, 0, cor(b_i, e_i))
      )  %>%
      mutate(
        trial = pars$try.no, 
        p0    = pars$p0,
        h2    = pars$h2,
        var.z = pars$sig.z^2
      )
  },
  mc.cores = 4
) %>%
  do.call(rbind, .)

# Get means and variances for each component in each time step for each age
sim.sum = sim.out %>%
  group_by(trial, h2, var.z, p0, t) %>%
  mutate(
    nn = sum(n),
    p_k = n / nn
  ) %>%
  group_by(h2, var.z, p0, k = age, t) %>%
  summarise(
    n   = sum(n) / trys.per,
    p_k = sum(p_k) / trys.per,
    across(contains('_bar'), mean),
    across(contains('_var'), mean),
    rho_k = mean(rho, na.rm = TRUE),
    n.obs = n()
  ) %>%
  ungroup() %>%
  # Add longevity column (for plotting)
  mutate(
    long = factor(p0),
    long = factor(long, labels = c('high', 'medium', 'low')),
    hert = factor(paste0('h^2 == ', h2)),
    varn = factor(paste0('gamma^2 == ', var.z))
  )

head(sim.sum)

# Data frame of analytical expectations
analyticals = merge(
  x = pars %>%
    select(h2, sig.z, p0, s.max, r, sig.0, sig.a, sig.e, lstar),
  y = sim.out %>% 
    mutate(sig.z = sqrt(var.z)) %>% 
    group_by(p0, h2, sig.z) %>% 
    summarise(age = max(age)) %>% 
    ungroup()
) %>%
  uncount(weights = age + 1) %>%
  group_by(p0, h2, sig.z) %>%
  # k column will denote age
  mutate(k = (1:n()) - 1) %>%
  ungroup() %>%
  mutate(
    # Age distribution
    p_k     = (r / (1 + r)) * (s.max / lstar)^k / sqrt(1 + k*sig.0^2),
    # Phenotypic variance components
    sig.z_k = sig.0^2 / (1 + k*sig.0^2),
    sig.a_k = sig.a^2 * (1 + k*sig.e^2) / (1 + k*(sig.a^2 + sig.e^2)),
    sig.e_k = sig.e^2 * (1 + k*sig.a^2) / (1 + k*(sig.a^2 + sig.e^2)),
    rho_k   = -k * sqrt(sig.a^2 * sig.e^2 / ((1+k*sig.a^2) * (1 + k*sig.e^2)))
  ) %>%
  # Add longevity column (for plotting)
  mutate(
    long = factor(p0),
    long = factor(long, labels = c('high', 'medium', 'low')),
    hert = factor(paste0('h^2 == ', h2)),
    varn = factor(paste0('gamma^2 == ', sig.z^2))
  )
  
head(analyticals)

# Phenotpyic component means and population size across entire populations for
# each trial
n.pheno.means = sim.out %>% 
  group_by(h2, var.z, p0, t, trial) %>% 
  summarise(
    across(contains('bar'), ~ (sum(. * n) / sum(n))),
    n = sum(n)
  ) %>% 
  ungroup() %>%
  mutate(
    long = factor(p0),
    long = factor(long, labels = c('high', 'medium', 'low')),
    hert = factor(paste0('h^2 == ', h2)),
    varn = factor(paste0('gamma^2 == ', var.z))
  )

### Begin plots

# Population size plot
n.pheno.means %>%
  ggplot(aes(x = t, y = n)) +
  geom_line(aes(group = trial, colour = long), linewidth = 0.5) +
  # scale_colour_brewer(palette = 'Dark2', 'longevity') +
  scale_colour_manual(values = c("#999999", "#56B4E9", "#E69F00"), 'longevity') +
  scale_x_continuous(breaks = (0:2)*5) +
  labs(x = 'Time step', y = 'Population size') +
  facet_grid(cols = vars(hert), rows = vars(varn), labeller = label_parsed) +
  theme(
    panel.background = element_blank(),
    legend.position = 'top'
  )

ggsave('analytical_validation/validation_figs/fig_v1_popsize.png',
       width = 8, height = 5)

# Age structure (all ages)
sim.sum %>%
  filter(t > 0, n.obs > 3) %>%
  ggplot(aes(x = k, y = p_k)) +
  geom_point(
    data = analyticals %>%
      mutate(
        hert = factor(paste0('h^2 == ', h2)),
        varn = factor(paste0('gamma^2 == ', sig.z^2))
      ),
    aes(x = k, y = p_k, colour = long),
    size = 2, shape = 21
  ) +
  geom_line(aes(group = interaction(p0, t), colour = long), linewidth = 0.1) +
  scale_colour_manual(values = c("#999999", "#56B4E9", "#E69F00"), 'longevity') +
  scale_fill_manual(values = c("#999999", "#56B4E9", "#E69F00"), 'longevity') +
  # scale_colour_manual(values = c('gray77', 'gray44', 'gray11')) +
  # scale_fill_manual(values = c('gray77', 'gray44', 'gray11')) +
  scale_y_log10() +
  labs(x = 'Age', y = 'Probability of being at age') +
  facet_grid(cols = vars(hert), rows = vars(varn), labeller = label_parsed) +
  theme(
    panel.background = element_blank(),
    legend.position = 'top'
  )

ggsave('analytical_validation/validation_figs/fig_v2_age_full.png',
       width = 8, height = 5)

# Age structure (only ages 0-9)
sim.sum %>%
  filter(t > 0, k < 10, n.obs > 3) %>%
  ggplot(aes(x = k, y = p_k)) +
  geom_point(
    data = analyticals %>% 
      filter(k < 10) %>% 
      mutate(
        hert = factor(paste0('h^2 == ', h2)),
        varn = factor(paste0('gamma^2 == ', sig.z^2))
      ),
    aes(x = k, y = p_k, colour = long),
    shape = 21, size = 2
  ) +
  geom_line(aes(group = interaction(p0, t), colour = long), linewidth = 0.1) +
  scale_colour_manual(values = c("#999999", "#56B4E9", "#E69F00"), 'longevity') +
  scale_fill_manual(values = c("#999999", "#56B4E9", "#E69F00"), 'longevity') +
  scale_x_continuous(breaks = (0:3)*3) +
  scale_y_log10() +
  labs(x = 'Age', y = 'Probability of being at age') +
  facet_grid(cols = vars(hert), rows = vars(varn), labeller = label_parsed) +
  theme(
    panel.background = element_blank(),
    legend.position = 'top'
  )

ggsave('analytical_validation/validation_figs/fig_v3_age_0-9.png',
       width = 8, height = 5)

# Mean breeding value over time (should average around zero)
n.pheno.means %>%
  ggplot(aes(x = t, y = b_i_bar, group = trial)) + 
  geom_line(aes(colour = long), linewidth = 0.5) + 
  labs(x = 'Time step', y = expression(Mean ~ breeding ~ value ~ bar(b)[t])) +
  scale_x_continuous(breaks = (0:2)*5) +
  scale_colour_manual(values = c("#999999", "#56B4E9", "#E69F00"), 'longevity') +
  facet_grid(
    cols = vars(hert), rows = vars(varn),
    labeller = label_parsed
  ) +
  theme(
    panel.background = element_blank(),
    legend.position = 'top'
  )

ggsave('analytical_validation/validation_figs/fig_v4_bbar_steady.png',
       width = 8, height = 5)

# Mean non-environmental component (should average around zero)
n.pheno.means %>%
  ggplot(aes(x = t, y = e_i_bar, group = trial)) + 
  geom_line(aes(colour = long), linewidth = 0.5) + 
  labs(x = 'Time step', y = expression(Mean ~ environmental ~ component ~ bar(e)[t])) +
  scale_x_continuous(breaks = (0:2)*5) +
  scale_colour_manual(values = c("#999999", "#56B4E9", "#E69F00"), 'longevity') +
  facet_grid(
    cols = vars(hert), rows = vars(varn),
    labeller = label_parsed
  ) +
  theme(
    panel.background = element_blank(),
    legend.position = 'top'
  )

ggsave('analytical_validation/validation_figs/fig_v5_ebar_steady.png',
       width = 8, height = 5)

# Mean phenotype over time (should average around zero)
n.pheno.means %>%
  ggplot(aes(x = t, y = z_i_bar, group = trial)) + 
  geom_line(aes(colour = long), linewidth = 0.5) + 
  labs(x = 'Time step', y = expression(Mean ~ phenotype ~ bar(z)[t])) +
  scale_x_continuous(breaks = (0:2)*5) +
  scale_colour_manual(values = c("#999999", "#56B4E9", "#E69F00"), 'longevity') +
  facet_grid(
    cols = vars(hert), rows = vars(varn),
    labeller = label_parsed
  ) +
  theme(
    panel.background = element_blank(),
    legend.position = 'top'
  )

ggsave('analytical_validation/validation_figs/fig_v6_zbar_steady.png',
       width = 8, height = 5)

# Additive genetic variance
sim.sum %>%
  filter(t > 0, k < 10, n.obs > 3) %>%
  mutate(
    h2 = factor(h2),
    long = factor(paste(long, 'longevity')),
    long = factor(long, levels = levels(long)[c(2, 3, 1)])
  ) %>%
  ggplot(aes(x = k, y = b_i_var)) +
  geom_line(aes(group = interaction(h2, t), colour = h2), linewidth = 0.1) +
  geom_point(
    data = analyticals %>% 
      filter(k < 10) %>% 
      mutate(
        var.z = sig.z^2, 
        h2 = factor(h2),
        long = factor(paste(long, 'longevity')),
        long = factor(long, levels = levels(long)[c(2, 3, 1)])
      ),
    aes(x = k, y = sig.a_k, colour = h2),
    shape = 21, size = 2
  ) +
  # scale_colour_brewer(palette = 'Dark2', 'longevity') +
  # scale_fill_brewer(palette = 'Dark2', 'longevity') +
  scale_colour_manual(values = c('gray77', 'gray44', 'gray11'), expression(h^2)) +
  scale_x_continuous(breaks = (0:3)*3) +
  labs(x = 'Age', y = expression(gamma[a]^2 ~ at ~ age)) +
  facet_grid(
    cols = vars(long), rows = vars(varn), 
    labeller = labeller(long = label_value, varn = label_parsed)
  ) +
  theme(
    panel.background = element_blank(),
    legend.position = 'top'
  )

ggsave('analytical_validation/validation_figs/fig_v7_bvar_age_steady.png',
       width = 8, height = 5)

# Environmental phenotypic variance
sim.sum %>%
  filter(t > 0, k < 10, n.obs > 3) %>%
  mutate(
    h2 = factor(h2),
    long = factor(paste(long, 'longevity')),
    long = factor(long, levels = levels(long)[c(2, 3, 1)])
  ) %>%
  ggplot(aes(x = k, y = e_i_var)) +
  geom_line(aes(group = interaction(h2, t), colour = h2), linewidth = 0.1) +
  geom_point(
    data = analyticals %>%
      filter(k < 10) %>% 
      mutate(
        var.z = sig.z^2, 
        h2 = factor(h2),
        long = factor(paste(long, 'longevity')),
        long = factor(long, levels = levels(long)[c(2, 3, 1)])
      ),
    aes(x = k, y = sig.e_k, colour = h2),
    size = 2, shape = 21
  ) +
  # scale_colour_brewer(palette = 'Dark2', 'longevity') +
  # scale_fill_brewer(palette = 'Dark2', 'longevity') +
  scale_colour_manual(values = c('gray77', 'gray44', 'gray11'), expression(h^2)) + 
  scale_x_continuous(breaks = (0:3)*3) +
  labs(x = 'Age', y = expression(gamma[e]^2 ~ at ~ age)) +
  facet_grid(
    cols = vars(long), rows = vars(varn), 
    labeller = labeller(long = label_value, varn = label_parsed)
  ) +
  theme(
    panel.background = element_blank(),
    legend.position = 'top'
  )

ggsave('analytical_validation/validation_figs/fig_v8_evar_age_steady.png',
       width = 8, height = 5)

# Phenotypic variance
sim.sum %>%
  filter(t > 0, k < 10, n.obs > 3) %>%
  mutate(
    h2 = factor(h2),
    long = factor(paste(long, 'longevity')),
    long = factor(long, levels = levels(long)[c(2, 3, 1)])
  ) %>%
  ggplot(aes(x = k, y = z_i_var)) +
  geom_line(aes(group = interaction(h2, t), colour = h2), linewidth = 0.1) +
  geom_point(
    data = analyticals %>%
      filter(k < 10) %>% 
      mutate(
        var.z = sig.z^2, 
        h2 = factor(h2),
        long = factor(paste(long, 'longevity')),
        long = factor(long, levels = levels(long)[c(2, 3, 1)])
      ),
    aes(x = k, y = sig.z_k, colour = h2),
    size = 2, shape = 21
  ) +
  # scale_colour_brewer(palette = 'Dark2', 'longevity') +
  # scale_fill_brewer(palette = 'Dark2', 'longevity') +
  scale_colour_manual(values = c('gray77', 'gray44', 'gray11'), expression(h^2)) +
  scale_x_continuous(breaks = (0:3)*3) +
  labs(x = 'Age', y = expression(gamma^2 ~ at ~ age)) +
  facet_grid(
    cols = vars(long), rows = vars(varn), 
    labeller = labeller(long = label_value, varn = label_parsed)
  ) +
  theme(
    panel.background = element_blank(),
    legend.position = 'top'
  )

ggsave('analytical_validation/validation_figs/fig_v9_zvar_age_steady.png',
       width = 8, height = 5)


# Breeding value-environmental component correlation

sim.sum %>%
  filter(t > 0, k < 10, h2 < 1, n.obs > 3) %>%
  mutate(
    h2 = factor(h2),
    long = factor(paste(long, 'longevity')),
    long = factor(long, levels = levels(long)[c(2, 3, 1)])
  ) %>%
  ggplot(aes(x = k, y = rho_k)) +
  geom_line(aes(group = interaction(h2, t), colour = h2), linewidth = 0.1) +
  geom_point(
    data = analyticals %>%
      filter(k < 10, h2 < 1) %>% 
      mutate(
        var.z = sig.z^2, 
        h2 = factor(h2),
        long = factor(paste(long, 'longevity')),
        long = factor(long, levels = levels(long)[c(2, 3, 1)])
      ),
    aes(x = k, y = rho_k, colour = h2),
    size = 2, shape = 21
  ) +
  # scale_colour_brewer(palette = 'Dark2', 'longevity') +
  # scale_fill_brewer(palette = 'Dark2', 'longevity') +
  scale_colour_manual(values = c('gray77', 'gray11'), expression(h^2)) +
  scale_x_continuous(breaks = (0:3)*3) +
  labs(x = 'Age', y = expression(gamma^2 ~ at ~ age)) +
  facet_grid(
    cols = vars(long), rows = vars(varn), 
    labeller = labeller(long = label_value, varn = label_parsed)
  ) +
  theme(
    panel.background = element_blank(),
    legend.position = 'top'
  )

ggsave('analytical_validation/validation_figs/fig_v10_rhok_age_steady.png',
       width = 8, height = 5)


### SesssionInfo (15 Dec 2023)

# R version 4.3.0 (2023-04-21)
# Platform: aarch64-apple-darwin20 (64-bit)
# Running under: macOS Big Sur 11.2.3
# 
# Matrix products: default
# BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
# LAPACK: /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.11.0
# 
# locale:
# [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# time zone: America/Vancouver
# tzcode source: internal
# 
# attached base packages:
# [1] parallel  stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
# [1] cowplot_1.1.1 mc2d_0.1-22   mvtnorm_1.1-3 tidyr_1.3.0   dplyr_1.1.3   ggplot2_3.4.3
# 
# loaded via a namespace (and not attached):
# [1] crayon_1.5.2       vctrs_0.6.3        cli_3.6.1          rlang_1.1.1       
# [5] purrr_1.0.2        generics_0.1.3     glue_1.6.2         labeling_0.4.3    
# [9] colorspace_2.1-0   scales_1.2.1       fansi_1.0.4        grid_4.3.0        
# [13] munsell_0.5.0      tibble_3.2.1       lifecycle_1.0.3    compiler_4.3.0    
# [17] RColorBrewer_1.1-3 pkgconfig_2.0.3    rstudioapi_0.15.0  farver_2.1.1      
# [21] R6_2.5.1           tidyselect_1.2.0   utf8_1.2.3         pillar_1.9.0      
# [25] magrittr_2.0.3     tools_4.3.0        withr_2.5.0        gtable_0.3.4  
