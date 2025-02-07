##########
# Simulated experiment
# Here: trying to capture age structure for populations
# - init 6 Mar 2023
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

pars = expand.grid(
  # (Equilibrium) size of initial cohort
  s.max = c(0.1, 0.5, 0.9)
) %>%
  mutate(
    # Maximum annual growth rate
    l.max = 1.4,
    # Phenotypic variance
    sig.z = sqrt(0.4),
    # Mean fecundity
    r     = (l.max / s.max) - 1
  )

pars = cbind(
  pars,
  pars %>%
    split(~ s.max + sig.z) %>%
    map(\(x) newt.method.2d(1.1, .5, x$s.max, x$r, x$sig.z^2, 1e-5)) %>%
    do.call(rbind, .)
)

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

# Run simulations
set.seed(233024)

sim.out2 = mclapply(
  pars %>% uncount(trys.per) %>% mutate(try.no = 1:(nrow(.))) %>% split(.$try.no),
  function(pars) {
    sim(pars, theta.t = 0, init.rows = 50 * 20000) %>%
      mutate(
        e_i = z_i - b_i
      ) %>%
      group_by(t) %>%
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
        s.max = pars$s.max,
        lstar = pars$lstar,
        h2    = pars$h2
      )
  },
  mc.cores = 16
) %>%
  do.call(rbind, .)

write.csv(
  sim.out2,
  file = 'run_sims/out/equal_lambda/sim_results_m1_phenotypic_components.csv',
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
