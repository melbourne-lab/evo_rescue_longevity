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
# - init 11 Dec. 2023, updated 6 Feb. 2025 for revision
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
trys.per = 20

# Parameters

# Parameters

pars = expand.grid(
  # Longevity groups defined by s.max
  s.max = c(0.1, 0.5, 0.9),
  # Genetic fitness groups defined by 
  sig.z = sqrt(c(0.1, 0.25, 0.4))
) %>%
  mutate(
    # Max lifetime fitness
    w.max = 3,
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
    timesteps = 10,
    # Initial pouplation size
    n.pop0 = 5000,
    # Ceiling-type carrying capacity
    kceil = 5000,
    # Initial population distance from phenotypic optimum (i.e., magnitude of
    # environmental shift)
    gbar0 = 0,
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

# Run simulations
set.seed(641876)

sim.out = mclapply(
  pars %>% uncount(trys.per) %>% mutate(try.no = 1:(nrow(.))) %>% split(.$try.no),
  function(pars) {
    sim(pars, theta.t = 0, init.rows = 10 * 5000) %>%
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
        s.max = pars$s.max,
        h2    = pars$h2,
        var.z = pars$sig.z^2
      )
  },
  mc.cores = 6
) %>%
  do.call(rbind, .)

# Get means and variances for each component in each time step for each age
sim.sum = sim.out %>%
  group_by(trial, h2, var.z, s.max, t) %>%
  mutate(
    nn = sum(n),
    p_k = n / nn
  ) %>%
  group_by(h2, var.z, s.max, k = age, t) %>%
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
    long = factor(s.max, labels = c('low', 'medium', 'high')),
    hert = factor(paste0('h^2 == ', h2)),
    varn = factor(paste0('gamma^2 == ', var.z))
  )

head(sim.sum)

# Data frame of analytical expectations
analyticals = merge(
  x = pars %>%
    select(h2, sig.z, s.max, r, sig.0, sig.a, sig.e, lstar),
  y = sim.out %>% 
    mutate(sig.z = sqrt(var.z)) %>% 
    group_by(s.max, h2, sig.z) %>% 
    summarise(age = max(age)) %>% 
    ungroup()
) %>%
  uncount(weights = age + 1) %>%
  group_by(s.max, h2, sig.z) %>%
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
    long = factor(s.max, labels = c('low', 'medium', 'high')),
    hert = factor(paste0('h^2 == ', h2)),
    varn = factor(paste0('gamma^2 == ', sig.z^2))
  )
  
head(analyticals)

# Phenotpyic component means and population size across entire populations for
# each trial
n.pheno.means = sim.out %>% 
  group_by(h2, var.z, s.max, t, trial) %>% 
  summarise(
    across(contains('bar'), ~ (sum(. * n) / sum(n))),
    n = sum(n)
  ) %>% 
  ungroup() %>%
  mutate(
    long = factor(s.max, labels = c('low', 'medium', 'high')),
    hert = factor(paste0('h^2 == ', h2)),
    varn = factor(paste0('gamma^2 == ', var.z))
  )

### Begin plots

# Population size plot is no longer useful under a ceiling-type density dependence
# # Population size plot
# n.pheno.means %>%
#   ggplot(aes(x = t, y = n)) +
#   geom_line(aes(group = trial, colour = long), linewidth = 0.5) +
#   # scale_colour_brewer(palette = 'Dark2', 'longevity') +
#   scale_colour_manual(values = c("#999999", "#56B4E9", "#E69F00"), 'longevity') +
#   scale_x_continuous(breaks = (0:2)*5) +
#   labs(x = 'Time step', y = 'Population size') +
#   facet_grid(cols = vars(hert), rows = vars(varn), labeller = label_parsed) +
#   theme(
#     panel.background = element_blank(),
#     legend.position = 'top'
#   )
# 
# ggsave('analytical_validation/validation_figs/fig_v1_popsize.png',
#        width = 8, height = 5)

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
    size = 2, shape = 21, alpha = 0.5
  ) +
  geom_line(aes(group = interaction(s.max, t), colour = long), linewidth = 0.2) +
  scale_colour_manual(
    values = c("#E69F00", "#56B4E9", "#999999"), 
    breaks = c('low', 'medium', 'high'),
    'longevity'
  ) +
  scale_fill_manual(
    values = c("#E69F00", "#56B4E9", "#999999"), 
    breaks = c('low', 'medium', 'high'),
    'longevity'
  ) +
  # scale_colour_manual(values = c("#999999", "#56B4E9", "#E69F00"), 'longevity') +
  # scale_fill_manual(values = c("#999999", "#56B4E9", "#E69F00"), 'longevity') +
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
    shape = 21, size = 2, alpha = 0.5
  ) +
  geom_line(aes(group = interaction(s.max, t), colour = long), linewidth = 0.2) +
  scale_colour_manual(
    values = c("#E69F00", "#56B4E9", "#999999"), 
    breaks = c('low', 'medium', 'high'),
    'longevity'
  ) +
  scale_fill_manual(
    values = c("#E69F00", "#56B4E9", "#999999"), 
    breaks = c('low', 'medium', 'high'),
    'longevity'
  ) +
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
  scale_colour_manual(
    values = c("#E69F00", "#56B4E9", "#999999"), 
    breaks = c('low', 'medium', 'high'),
    'longevity'
  ) +
  scale_fill_manual(
    values = c("#E69F00", "#56B4E9", "#999999"), 
    breaks = c('low', 'medium', 'high'),
    'longevity'
  ) +
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
  scale_colour_manual(
    values = c("#E69F00", "#56B4E9", "#999999"), 
    breaks = c('low', 'medium', 'high'),
    'longevity'
  ) +
  scale_fill_manual(
    values = c("#E69F00", "#56B4E9", "#999999"), 
    breaks = c('low', 'medium', 'high'),
    'longevity'
  ) +
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
  scale_colour_manual(
    values = c("#E69F00", "#56B4E9", "#999999"), 
    breaks = c('low', 'medium', 'high'),
    'longevity'
  ) +
  scale_fill_manual(
    values = c("#E69F00", "#56B4E9", "#999999"), 
    breaks = c('low', 'medium', 'high'),
    'longevity'
  ) +
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


### SessionInfo (6 Feb. 2025)

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
# time zone: America/Los_Angeles
# tzcode source: internal
# 
# attached base packages:
# [1] parallel  stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
# [1] cowplot_1.1.1 mc2d_0.1-22   mvtnorm_1.1-3 purrr_1.0.2   tidyr_1.3.0   dplyr_1.1.3  
# [7] ggplot2_3.5.1
# 
# loaded via a namespace (and not attached):
# [1] crayon_1.5.2      vctrs_0.6.3       cli_3.6.1         rlang_1.1.4      
# [5] generics_0.1.3    textshaping_0.3.6 glue_1.6.2        labeling_0.4.3   
# [9] colorspace_2.1-0  ragg_1.2.5        scales_1.3.0      fansi_1.0.4      
# [13] grid_4.3.0        munsell_0.5.0     tibble_3.2.1      lifecycle_1.0.3  
# [17] compiler_4.3.0    pkgconfig_2.0.3   rstudioapi_0.15.0 systemfonts_1.0.4
# [21] farver_2.1.1      viridisLite_0.4.2 R6_2.5.1          tidyselect_1.2.0 
# [25] utf8_1.2.3        pillar_1.9.0      magrittr_2.0.3    tools_4.3.0      
# [29] withr_2.5.0       gtable_0.3.4   
