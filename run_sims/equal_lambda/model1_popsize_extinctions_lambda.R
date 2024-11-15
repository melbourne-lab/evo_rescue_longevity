##########
# Simulated experiment
# Here: all life histories run with equivalent lambda hat (i.e., equivalent max
# population growth rates)
# - init 6 Mar 2023
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
trys.per = 500

# Parameters

pars = expand.grid(
  # (Equilibrium) size of initial cohort
  s.max = c(0.1, 0.5, 0.9),
  # Phenotypic variance
  sig.z = sqrt(c(0.1, 0.25, 0.4)),
  # Heritability of fitness
  h2    = c(.25, .5, 1)
) %>%
  # Demographic rates
  mutate(
    # Maximum expected lifetime fitness
    l.max = 1.4,
    # Mean fecundity
    r     = (l.max / s.max) - 1,
    # Equilibrium population growth rate
    lstar = l.max / sqrt(1 + sig.z^2),
    # Initial population size
    n.pop0 = 20000,
    # Strength of density dependence
    # alpha = log(lstar) / n.pop0,
    # Ceciling-type carrying capacity just in case
    kceil = 20000, # 30000,
    p0    = r / (1 + r)
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
  mutate(timesteps = 100)

# Run simulations
set.seed(233024)

sim.out = mclapply(
  pars %>% uncount(trys.per) %>% mutate(try.no = 1:(nrow(.))) %>% split(.$try.no),
  function(pars) {
    sim(pars, theta.t = 0, init.rows = 100 * 20000) %>%
      group_by(t) %>%
      summarise(n = n())  %>%
      mutate(
        trial = pars$try.no, 
        p0    = pars$p0,
        h2    = pars$h2,
        var.z = pars$sig.z^2
      )
  },
  mc.cores = 16
) %>%
  do.call(rbind, .)

sim.r = sim.out %>%
  group_by(trial) %>%
  mutate(log.lam = c(diff(log(n)), NA)) %>%
  filter(!is.na(log.lam)) %>%
  group_by(t, p0, var.z, h2) %>%
  summarise(
    llbar = mean(log.lam),
    llvar = var(log.lam),
    n     = mean(n),
    nn = n()
  )

sim.n.all = sim.out %>%
  merge(
    expand.grid(
      t     = 0:pars$timesteps[1],
      trial = 1:(nrow(pars) * trys.per)
    ),
    all.y = TRUE
  ) %>%
  mutate(n = ifelse(is.na(n), 0, n)) %>%
  group_by(trial) %>%
  arrange(t) %>%
  mutate(
    p0 =    p0[1],
    h2 =    h2[1],
    var.z = var.z[1]
  ) %>%
  group_by(p0, var.z, h2, t, trial) %>%
  summarise(n = sum(n)) %>%
  group_by(p0, var.z, h2, t) %>%
  summarise(
    nbar = mean(n),
    nvar = var(n),
    psrv = mean(n > 0),
    nn   = n()
  )

sim.n.surv = sim.out %>%
  group_by(trial) %>%
  filter(max(t) == pars$timesteps[1]) %>%
  group_by(h2, var.z, p0, t) %>%
  summarise(
    nbar = mean(n),
    nvar = var(n),
    nn   = n()
  )

# Get 20 trials per parameter combination
sim.disagg = sim.out %>% filter((trial %% trys.per) < 20)

write.csv(
  sim.r,
  row.names = FALSE,
  file = 'run_sims/out/equal_lambda/sim_results_m1_sizes_r.csv'
)

write.csv(
  sim.n.all,
  row.names = FALSE,
  file = 'run_sims/out/equal_lambda/sim_results_m1_allsizes_n.csv'
)

write.csv(
  sim.n.surv,
  row.names = FALSE,
  file = 'run_sims/out/equal_lambda/sim_results_m1_survsizes_n.csv'
)

write.csv(
  sim.disagg,
  row.names = FALSE,
  file = 'run_sims/out/equal_lambda/sim_results_m1_disaggregated_n.csv'
)

### SesssionInfo (13 Nov 2023 - final run)

# R version 4.1.2 (2021-11-01)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Ubuntu 18.04.6 LTS
# 
# Matrix products: default
# BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.7.1
# LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.7.1
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
# [1] tidyr_1.1.3   dplyr_1.0.7   ggplot2_3.4.0
# 
# loaded via a namespace (and not attached):
# [1] fansi_0.5.0      withr_2.5.0      crayon_1.4.1     utf8_1.2.2      
# [5] grid_4.1.2       R6_2.5.0         lifecycle_1.0.3  gtable_0.3.0    
# [9] magrittr_2.0.1   scales_1.2.1     pillar_1.6.2     rlang_1.0.6     
# [13] cli_3.6.0        generics_0.1.0   vctrs_0.5.2      ellipsis_0.3.2  
# [17] glue_1.4.2       purrr_0.3.4      munsell_0.5.0    compiler_4.1.2  
# [21] pkgconfig_2.0.3  colorspace_2.0-2 tidyselect_1.1.1 tibble_3.1.3
