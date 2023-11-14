##########
# Simulated experiment
# Here: trying to capture age structure for populations
# SN - init 6 Mar 2023
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
trys.per = 200

# Parameters
pars = expand.grid(
  # (Equilibrium) size of initial cohort
  s.max = c(0.1, 0.5, 0.9),
  # Heritability of fitness
  h2    = c(.25, .5, 1)
) %>%
  # Demographic rates
  mutate(
    # Maximum expected lifetime fitness
    w.max = 3,
    # Gamma squared (pheno variance / sel pressure)
    sig.z = sqrt(.4),
    # Equilibrium lifetime fitness
    wstar = w.max * (1 - s.max) / (sqrt(1 + sig.z^2) - s.max),
    # Mean fecundity
    r     = w.max * (1 - s.max) / s.max,
    # Equilibrium population growth rate
    lstar = (s.max + w.max * (1 - s.max)) / (s.max + (w.max/wstar) * (1 - s.max)),
    # Initial population size
    n.pop0 = 20000,
    # Strength of density dependence
    alpha = log(lstar) / n.pop0,
    # Ceciling-type carrying capacity just in case
    kceil = 30000,
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
    # Non-inherited standard dxeviation in new cohorts
    sig.e = sqrt((1-h2) * sig.0^2),
    # Population-wide breeding value standard deviation
    sig.p = sqrt(gamma.a.calc(sig.a^2, s.max / lstar, r, sig.e^2)),
    mu    = 1,
    sig.m = sqrt(wfitn^2 * (sig.a^2 - (sig.p^2 - p0*sig.a^2)/(1-p0))),
    gbar0 = 2
  ) %>%
  ungroup() %>%
  # Other junk
  mutate(timesteps = 50)

# Run simulations
set.seed(4523)

sim.out2 = mclapply(
  pars %>% uncount(trys.per) %>% mutate(try.no = 1:(nrow(.))) %>% split(.$try.no),
  function(pars) {
    sim(pars, theta.t = 0, init.rows = 50 * 30000) %>%
      mutate(
        e_i = z_i - b_i
      ) %>%
      group_by(t, age) %>%
      summarise(
        n = n(),
        r = sum(r_i),
        bbar = mean(b_i),
        bvar = var(b_i),
        ebar = mean(e_i),
        evar = var(e_i),
        zbar = mean(z_i),
        zvar = var(z_i)
      )  %>%
      mutate(
        trial = pars$try.no, 
        p0    = pars$p0,
        lstar = pars$lstar,
        h2    = pars$h2
      )
  },
  mc.cores = 12
) %>%
  do.call(rbind, .)

write.csv(
  sim.out2,
  file = 'run_sims/out/sim_results_m1_ages.csv',
  row.names = FALSE
)

# R Session info for final run (14 Nov. 2023)

# R version 4.1.2 (2021-11-01)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Ubuntu 18.04.6 LTS
# 
# Matrix products: default
# BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.7.1
# LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.7.1
# 
# locale:
# [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8       
# [4] LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
# [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C              
# [10] LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
# 
# attached base packages:
# [1] parallel  stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
# [1] mc2d_0.1-22   mvtnorm_1.1-3 tidyr_1.1.3   dplyr_1.0.7   ggplot2_3.4.0
# 
# loaded via a namespace (and not attached):
# [1] rstudioapi_0.13   magrittr_2.0.1    tidyselect_1.1.1  munsell_0.5.0    
# [5] viridisLite_0.4.0 colorspace_2.0-2  R6_2.5.0          rlang_1.0.6      
# [9] fansi_0.5.0       tools_4.1.2       grid_4.1.2        gtable_0.3.0     
# [13] utf8_1.2.2        cli_3.6.0         withr_2.5.0       ellipsis_0.3.2   
# [17] tibble_3.1.3      lifecycle_1.0.3   crayon_1.4.1      farver_2.1.0     
# [21] purrr_0.3.4       vctrs_0.5.2       glue_1.4.2        labeling_0.4.2   
# [25] compiler_4.1.2    pillar_1.6.2      generics_0.1.0    scales_1.2.1     
# [29] pkgconfig_2.0.3  
