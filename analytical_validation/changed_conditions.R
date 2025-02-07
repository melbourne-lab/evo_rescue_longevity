##########
# Script to validate analytical expressions following environmental change
# - Simulations were run following environmental change 
#   (population phenotype is 1 unit away from novel optimum)
# - Parameter combinations are mostly same as in main text; 
#   - differences: population phenotypic variances are larger to accentuate 
#     differences among treatments, and sims are run only for heritability of
#     0.5
# - Each trial is run for three time steps
# - I aggregate data in two ways:
#   - mean phenotypic component variance for age class (only for age-classes 
#     that were alive before the environmental shift)
#     (these should remain constant over time - overlapping lines)
#   - mean phenotypic component mean for each cohort as they age
#     these are presented as a time-series as the validated expressions
#     demonstrate cohort dynamics through time
# - Script contains plots for: 
#   - Mean of each phenotypic component over time (per cohort)
#   - Variance of each phenotypic component (per age)
# - init 18 Dec. 2023, revised 6 Feb. 2025 for resubmission
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
trys.per = 100

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
  mutate(
    sig.0 = sqrt(sig.0),
    h2 = 0.5
  ) %>%
  # Estimate genetic variance and environmental component variance in newborn
  # cohorts
  mutate(
    sig.a = sqrt(h2 * sig.0^2), 
    sig.e = sqrt((1 - h2) * sig.0^2)
  ) %>%
  # Add other necessary parameters as needed
  mutate(
    # Length of simulations
    timesteps = 3,
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
set.seed(820991)

sim.out = mclapply(
  pars %>% uncount(trys.per) %>% mutate(try.no = 1:(nrow(.))) %>% split(.$try.no),
  function(pars) {
    sim(pars, theta.t = 0, init.rows = 5 * 5000) %>%
      mutate(e_i = z_i - b_i) %>%
      group_by(t, age) %>%
      summarise(
        n = n(),
        across(
          c(b_i, e_i, z_i), 
          .fns = list(bar = mean, var = ~ ifelse(n == 1, 0, var(.x) * (n/(n-1))))
        ),
        rho = cor(b_i, e_i)
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

# Simulation summary by cohort
# (to get dynamics of cohort after env. shift)
coh.sum = sim.out %>%
  filter(age >= t) %>%
  mutate(coh = age - t) %>%
  group_by(h2, var.z, s.max, k = coh, t) %>%
  summarise(
    across(contains('_bar'), mean),
    n = n()
  ) %>%
  ungroup() %>%
  mutate(
    long = factor(s.max, labels = paste(c('low', 'medium', 'high'), 'longevity')),
    # long = factor(long, levels = levels(long)[3:1]),
    hert = factor(paste0('h^2 == ', h2)),
    varn = factor(paste0('gamma^2 == ', var.z))
  )
  
head(coh.sum)

  
# Analytical solutions
# Data frame of analytical expectations
coh.analyticals = merge(
  x = pars %>%
    select(h2, sig.z, s.max, s.max, r, sig.0, sig.a, sig.e, lstar, gbar0),
  y = sim.out %>% 
    mutate(sig.z = sqrt(var.z)) %>% 
    group_by(s.max, h2, sig.z) %>% 
    summarise(coh = max(age)) %>% 
    ungroup()
) %>%
  # Add age column
  uncount(weights = coh + 1) %>%
  group_by(s.max, h2, sig.z) %>%
  # k column will denote age
  mutate(k = (1:n()) - 1) %>%
  ungroup() %>%
  # Add time step column
  uncount(weights = pars$timesteps[1] + 1) %>%
  group_by(s.max, h2, sig.z, k) %>%
  mutate(t = (1:n()) - 1) %>%
  ungroup() %>%
  # Add analytical estimations
  mutate(
    # Phenotypic mean components
    bar.b_k = gbar0 * (1 - (t*sig.a^2) / (1 + (k + t) * (sig.a^2 + sig.e^2))),
    bar.e_k = gbar0 * (0 - (t*sig.e^2) / (1 + (k + t) * (sig.a^2 + sig.e^2))),
    bar.z_k = gbar0 * (1 - (t * (sig.a^2 + sig.e^2)) / (1 + (k+t) * (sig.a^2 + sig.e^2))),
    # Phenotypic variance components
    sig.z_k = sig.0^2 / (1 + k*sig.0^2),
    sig.a_k = sig.a^2 * (1 + k*sig.e^2) / (1 + k*(sig.a^2 + sig.e^2)),
    sig.e_k = sig.e^2 * (1 + k*sig.a^2) / (1 + k*(sig.a^2 + sig.e^2)),
    rho_k   = -k * sqrt(sig.a^2 * sig.e^2 / ((1+k*sig.a^2) * (1 + k*sig.e^2)))
  ) %>%
  # Add longevity column (for plotting)
  mutate(
    long = factor(s.max, labels = paste(c('low', 'medium', 'high'), 'longevity')),
    # long = factor(long, levels = levels(long)[3:1]),
    hert = factor(paste0('h^2 == ', h2)),
    varn = factor(paste0('gamma^2 == ', sig.z^2))
  )

head(coh.analyticals)

# Mean breeding value 
coh.sum %>%
  filter(k < 5, k >= t, h2 %in% 0.5, n > 10) %>%
  ggplot(aes(x = t, y = b_i_bar)) +
  geom_line(aes(group = k, colour = k), linewidth = 0.5) +
  geom_point(
    data = coh.analyticals %>% filter(h2 %in% 0.5, k < 5),
    aes(x = t, y = bar.b_k, group = k, colour = k),
    shape = 21, size = 2
  ) +
  scale_colour_gradient(low = 'gray88', high = 'gray11', 'age at\n environmental shift') +
  guides(colour = guide_coloursteps(nbin = 5)) +
  labs(
    x = expression(Time ~ step ~ t), 
    y = expression(Mean ~ age ~ class ~ breeding ~ value ~ bar(b)['k,k+t'])
  ) +
  facet_grid(
    cols = vars(long), rows = vars(varn), 
    labeller = labeller(long = label_value, varn = label_parsed)
  ) +
  theme(
    panel.background = element_blank(),
    legend.position = 'top',
    legend.title.align = 0.5
  )

ggsave('analytical_validation/validation_figs/fig_v11_bbar_change.png',
       width = 8, height = 5)

# Mean environmental portion of phenotype
coh.sum %>%
  filter(k < 5, k >= t, h2 %in% 0.5, n > 10) %>%
  ggplot(aes(x = t, y = e_i_bar)) +
  geom_line(aes(group = k, colour = k), linewidth = 0.5) +
  geom_point(
    data = coh.analyticals %>% filter(h2 %in% 0.5, k < 5),
    aes(x = t, y = bar.e_k, group = k, colour = k),
    shape = 21, size = 2
  ) +
  scale_colour_gradient(low = 'gray88', high = 'gray11', 'age at\n environmental shift') +
  guides(colour = guide_coloursteps(nbin = 5)) +
  labs(
    x = expression(Time ~ step ~ t), 
    y = expression(Mean ~ age ~ class ~ environmental ~ portion ~ of ~ phenotype ~ bar(e)['k,k+t'])
  ) +
  facet_grid(
    cols = vars(long), rows = vars(varn), 
    labeller = labeller(long = label_value, varn = label_parsed)
  ) +
  theme(
    panel.background = element_blank(),
    legend.position = 'top',
    legend.title.align = 0.5
  )

ggsave('analytical_validation/validation_figs/fig_v12_ebar_change.png',
       width = 8, height = 5)

# Mean phenotype
coh.sum %>%
  filter(k < 5, k >= t, h2 %in% 0.5, n > 10) %>%
  ggplot(aes(x = t, y = z_i_bar)) +
  geom_line(aes(group = k, colour = k), linewidth = 0.5) +
  geom_point(
    data = coh.analyticals %>% filter(h2 %in% 0.5, k < 5),
    aes(x = t, y = bar.z_k, group = k, colour = k),
    shape = 21, size = 2
  ) +
  scale_colour_gradient(low = 'gray88', high = 'gray11', 'age at\n environmental shift') +
  guides(colour = guide_coloursteps(nbin = 5)) +
  labs(
    x = expression(Time ~ step ~ t), 
    y = expression(Mean ~ age ~ class ~ phenotype ~ bar(z)['k,k+t'])
  ) +
  facet_grid(
    cols = vars(long), rows = vars(varn), 
    labeller = labeller(long = label_value, varn = label_parsed)
  ) +
  theme(
    panel.background = element_blank(),
    legend.position = 'top',
    legend.title.align = 0.5
  )

ggsave('analytical_validation/validation_figs/fig_v13_zbar_change.png',
       width = 8, height = 5)


# Age-summary
# (to get varinaces of age classes, not cohorts, after env. change)
age.sum = sim.out %>%
  group_by(h2, var.z, s.max, k = age, t) %>%
  summarise(
    across(contains('_var'), mean),
    rho_k = mean(rho),
    n = n()
  ) %>%
  ungroup() %>%
  mutate(
    long = factor(s.max, labels = paste(c('low', 'medium', 'high'), 'longevity')),
    # long = factor(long, levels = levels(long)[3:1]),
    hert = factor(paste0('h^2 == ', h2)),
    varn = factor(paste0('gamma^2 == ', var.z))
  )

# Age analytical solutions
age.analyticals = merge(
  x = pars %>%
    select(h2, sig.z, s.max, s.max, r, sig.0, sig.a, sig.e, lstar, gbar0),
  y = sim.out %>% 
    mutate(sig.z = sqrt(var.z)) %>% 
    group_by(s.max, h2, sig.z) %>% 
    summarise(age = max(age)) %>% 
    ungroup()
) %>%
  # Add age column
  uncount(weights = age + 1) %>%
  group_by(s.max, h2, sig.z) %>%
  # k column will denote age
  mutate(k = (1:n()) - 1) %>%
  ungroup() %>%
  # Add time step column
  uncount(weights = pars$timesteps[1] + 1) %>%
  group_by(s.max, h2, sig.z, k) %>%
  mutate(t = (1:n()) - 1) %>%
  ungroup() %>%
  # Add analytical estimations
  mutate(
    # Phenotypic variance components
    sig.z_k = sig.0^2 / (1 + k*sig.0^2),
    sig.a_k = sig.a^2 * (1 + k*sig.e^2) / (1 + k*(sig.a^2 + sig.e^2)),
    sig.e_k = sig.e^2 * (1 + k*sig.a^2) / (1 + k*(sig.a^2 + sig.e^2)),
    rho_k   = -k * sqrt(sig.a^2 * sig.e^2 / ((1+k*sig.a^2) * (1 + k*sig.e^2)))
  ) %>%
  # Add longevity column (for plotting)
  mutate(
    long = factor(s.max, labels = paste(c('low', 'medium', 'high'), 'longevity')),
    # long = factor(long, levels = levels(long)[3:1]),
    hert = factor(paste0('h^2 == ', h2)),
    varn = factor(paste0('gamma^2 == ', sig.z^2))
  )

# Additive genetic variance after environmental change
age.sum %>%
  filter(t > 0, k < 10, k >= t) %>%
  ggplot(aes(x = k, y = b_i_var)) +
  geom_line(aes(group = t), linewidth = 0.5) +
  geom_point(
    data = age.analyticals %>%  filter(k < 10, t > 0, k >= t),
    aes(x = k, y = sig.a_k),
    shape = 21, size = 2
  ) +
  # scale_colour_manual(values = c('gray77', 'gray44', 'gray11'), expression(h^2)) +
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

ggsave('analytical_validation/validation_figs/fig_v14_bvar_change.png',
       width = 8, height = 5)

# Environmental portion of phenotype variance after environmental change
age.sum %>%
  filter(t > 0, k < 10, k >= t) %>%
  ggplot(aes(x = k, y = e_i_var)) +
  geom_line(aes(group = t), linewidth = 0.5) +
  geom_point(
    data = age.analyticals %>%  filter(k < 10, t > 0, k >= t),
    aes(x = k, y = sig.e_k),
    shape = 21, size = 2
  ) +
  # scale_colour_manual(values = c('gray77', 'gray44', 'gray11'), expression(h^2)) +
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

ggsave('analytical_validation/validation_figs/fig_v15_evar_change.png',
       width = 8, height = 5)

# Phenotypic variance after environmental change
age.sum %>%
  filter(t > 0, k < 10, k >= t) %>%
  ggplot(aes(x = k, y = z_i_var)) +
  geom_line(aes(group = t), linewidth = 0.5) +
  geom_point(
    data = age.analyticals %>%  filter(k < 10, t > 0, k >= t),
    aes(x = k, y = sig.z_k),
    shape = 21, size = 2
  ) +
  # scale_colour_manual(values = c('gray77', 'gray44', 'gray11'), expression(h^2)) +
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

ggsave('analytical_validation/validation_figs/fig_v16_zvar_change.png',
       width = 8, height = 5)

age.sum %>%
  filter(t > 0, k < 10, k >= t, !is.na(rho_k)) %>%
  ggplot(aes(x = k, y = rho_k)) +
  geom_line(aes(group = t), linewidth = 0.5) +
  geom_point(
    data = age.analyticals %>%  filter(k < 10, t > 0, k >= t),
    aes(x = k, y = rho_k),
    shape = 21, size = 2
  ) +
  # scale_colour_manual(values = c('gray77', 'gray44', 'gray11'), expression(h^2)) +
  scale_x_continuous(breaks = (0:3)*3) +
  labs(x = 'Age', y = expression(rho ~ at ~ age)) +
  facet_grid(
    cols = vars(long), rows = vars(varn), 
    labeller = labeller(long = label_value, varn = label_parsed)
  ) +
  theme(
    panel.background = element_blank(),
    legend.position = 'top'
  )

ggsave('analytical_validation/validation_figs/fig_v17_rhok_change.png',
       width = 8, height = 5)


### Session info (6 Feb. 2025)
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
