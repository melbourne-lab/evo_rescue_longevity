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
# SN - init 18 Dec. 2023
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
  # (Equilibrium) size of initial cohort
  s.max = c(0.1, 0.5, 0.9),
  # Heritability of fitness
  h2    = .5,
  # Gamma squared (pheno variance / sel pressure)
  sig.z = sqrt(c(.3, .45, .6))
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
    gbar0 = 1
  ) %>%
  ungroup() %>%
  # Other junk
  mutate(timesteps = 3)

# Run simulations
set.seed(820991)

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
        rho = cor(b_i, e_i)
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

# Simulation summary by cohort
# (to get dynamics of cohort after env. shift)
coh.sum = sim.out %>%
  filter(age >= t) %>%
  mutate(coh = age - t) %>%
  group_by(h2, var.z, p0, k = coh, t) %>%
  summarise(
    across(contains('_bar'), mean),
    n = n()
  ) %>%
  ungroup() %>%
  mutate(
    long = factor(p0),
    long = factor(long, labels = paste(c('high', 'medium', 'low'), 'longevity')),
    long = factor(long, levels = levels(long)[3:1]),
    hert = factor(paste0('h^2 == ', h2)),
    varn = factor(paste0('gamma^2 == ', var.z))
  )
  
head(coh.sum)

  
# Analytical solutions
# Data frame of analytical expectations
coh.analyticals = merge(
  x = pars %>%
    select(h2, sig.z, p0, s.max, r, sig.0, sig.a, sig.e, lstar, gbar0),
  y = sim.out %>% 
    mutate(sig.z = sqrt(var.z)) %>% 
    group_by(p0, h2, sig.z) %>% 
    summarise(coh = max(age)) %>% 
    ungroup()
) %>%
  # Add age column
  uncount(weights = coh + 1) %>%
  group_by(p0, h2, sig.z) %>%
  # k column will denote age
  mutate(k = (1:n()) - 1) %>%
  ungroup() %>%
  # Add time step column
  uncount(weights = pars$timesteps[1] + 1) %>%
  group_by(p0, h2, sig.z, k) %>%
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
    long = factor(p0),
    long = factor(long, labels = paste(c('high', 'medium', 'low'), 'longevity')),
    long = factor(long, levels = levels(long)[3:1]),
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
  guides(colour = guide_colourbar(nbin = 5)) +
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
  guides(colour = guide_colourbar(nbin = 5)) +
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
  guides(colour = guide_colourbar(nbin = 5)) +
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
  group_by(h2, var.z, p0, k = age, t) %>%
  summarise(
    across(contains('_var'), mean),
    rho_k = mean(rho),
    n = n()
  ) %>%
  ungroup() %>%
  mutate(
    long = factor(p0),
    long = factor(long, labels = paste(c('high', 'medium', 'low'), 'longevity')),
    long = factor(long, levels = levels(long)[3:1]),
    hert = factor(paste0('h^2 == ', h2)),
    varn = factor(paste0('gamma^2 == ', var.z))
  )

# Age analytical solutions
age.analyticals = merge(
  x = pars %>%
    select(h2, sig.z, p0, s.max, r, sig.0, sig.a, sig.e, lstar, gbar0),
  y = sim.out %>% 
    mutate(sig.z = sqrt(var.z)) %>% 
    group_by(p0, h2, sig.z) %>% 
    summarise(age = max(age)) %>% 
    ungroup()
) %>%
  # Add age column
  uncount(weights = age + 1) %>%
  group_by(p0, h2, sig.z) %>%
  # k column will denote age
  mutate(k = (1:n()) - 1) %>%
  ungroup() %>%
  # Add time step column
  uncount(weights = pars$timesteps[1] + 1) %>%
  group_by(p0, h2, sig.z, k) %>%
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
    long = factor(p0),
    long = factor(long, labels = paste(c('high', 'medium', 'low'), 'longevity')),
    long = factor(long, levels = levels(long)[3:1]),
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
# [5] purrr_1.0.2        generics_0.1.3     textshaping_0.3.6  glue_1.6.2        
# [9] labeling_0.4.3     colorspace_2.1-0   ragg_1.2.5         scales_1.2.1      
# [13] fansi_1.0.4        grid_4.3.0         munsell_0.5.0      tibble_3.2.1      
# [17] lifecycle_1.0.3    compiler_4.3.0     RColorBrewer_1.1-3 pkgconfig_2.0.3   
# [21] rstudioapi_0.15.0  systemfonts_1.0.4  farver_2.1.1       R6_2.5.1          
# [25] tidyselect_1.2.0   utf8_1.2.3         pillar_1.9.0       magrittr_2.0.3    
# [29] tools_4.3.0        withr_2.5.0        gtable_0.3.4    
