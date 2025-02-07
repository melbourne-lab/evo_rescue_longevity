##########
# Simulated experiment
# Here: recording population-level phenotypic components
# - init 9 Mar 2023
##########

# Packages
library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)
library(parallel)

# Clear namespace
rm(list = ls())

# Load source code
source('model_source/sim_model1_functions.R')

### Load in parameters

# Trials per parameter combo
trys.per = 200

pars = data.frame(
  # Longevity groups defined by s.max
  s.max = c(0.1, 0.5, 0.9)
) %>%
  mutate(
    # Max lifetime fitness
    w.max = 3,
    # Population-wide phenotypic variance (one level here)
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
    timesteps = 50,
    # Initial population size
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

# # Parameters
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
#     # alpha = log(lstar) / n.pop0,
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

# Run simulations
set.seed(4523)

sim.out2 = mclapply(
  pars %>% uncount(trys.per) %>% mutate(try.no = 1:(nrow(.))) %>% split(.$try.no),
  function(pars) {
    sim(pars, theta.t = 0, init.rows = 50 * 30000) %>%
      mutate(e_i = z_i - b_i) %>%
      group_by(t) %>%
      summarise(
        n = n(),
        bbar = mean(b_i),
        bvar = var(b_i),
        bvarw = ifelse(any(age > 0), var(b_i[age > 0]), NA),
        ebar = mean(e_i),
        evar = var(e_i),
        evarw = ifelse(any(age > 0), var(e_i[age > 0]), NA),
        zbar = mean(z_i),
        zvar = var(z_i),
        zvarw = ifelse(any(age > 0), var(z_i[age > 0]), NA)
      )  %>%
      mutate(
        trial = pars$try.no, 
        s.max = pars$s.max,
        h2    = pars$h2
      )
  },
  mc.cores = 16
) %>%
  do.call(rbind, .)

write.csv(
  sim.out2,
  file = 'run_sims/out/sim_results_m1_phtype.csv',
  row.names = FALSE
)

### SesssionInfo (14 Nov 2024 - final run)

# R version 4.4.1 (2024-06-14)
# Platform: x86_64-conda-linux-gnu
# Running under: Ubuntu 18.04.6 LTS
# 
# Matrix products: default
# BLAS/LAPACK: /home/XXX/miniconda3/envs/r-evo_rescue/lib/libopenblasp-r0.3.27.so;  LAPACK version 3.12.0
# 
# locale:
# [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
# [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
# [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
# [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
# [9] LC_ADDRESS=C               LC_TELEPHONE=C            
# [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
# 
# attached base packages:
# [1] parallel  stats     graphics  grDevices utils     datasets  methods  
# [8] base     
# 
# other attached packages:
# [1] mc2d_0.2.1    mvtnorm_1.2-5 purrr_1.0.2   tidyr_1.3.1   dplyr_1.1.4  
# [6] ggplot2_3.5.1
# 
# loaded via a namespace (and not attached):
# [1] vctrs_0.6.5      cli_3.6.3        rlang_1.1.4      car_3.1-2       
# [5] generics_0.1.3   ggpubr_0.6.0     glue_1.7.0       backports_1.5.0 
# [9] colorspace_2.1-1 scales_1.3.0     fansi_1.0.6      grid_4.4.1      
# [13] abind_1.4-5      carData_3.0-5    rstatix_0.7.2    munsell_0.5.1   
# [17] tibble_3.2.1     lifecycle_1.0.4  ggsignif_0.6.4   compiler_4.4.1  
# [21] pkgconfig_2.0.3  R6_2.5.1         tidyselect_1.2.1 utf8_1.2.4      
# [25] pillar_1.9.0     magrittr_2.0.3   withr_3.0.1      gtable_0.3.5    
# [29] broom_1.0.6   
