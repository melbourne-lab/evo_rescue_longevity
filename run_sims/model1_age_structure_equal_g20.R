##########
# Simulated experiment
# Here: recording age structure
# Simulated experiment designed with equal initial genetic variance (gamma^2_0)
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

### Write two more wrapper functions

g1.lambda.calc = function(lam, s, r, g2_0) {
  # Inputs:
  # lam  = lambda* (lstar) guess to estimate with
  # s, r, g2_0 are specified, it's just a matter of estimating lambda
  
  # Initialize sums
  # g1 is age distribution error function:
  #   g1 = sum_k p_k - 1
  sum.g1 = 0
  # Derivative of above function wrt lambda*
  sum.dg1.lam = 0
  
  # Sum over ages
  for (k in 0:1e5) {
    
    # store a value for p_k term in age distribution
    pk = (r / (1 + r)) * (s / lam)^k  / sqrt(1 + k*g2_0)
    
    # g1 (age distribution) terms
    sum.g1 = sum.g1 + pk
    sum.dg1.lam = sum.dg1.lam + (pk * (-k / lam))
    
  }
  
  return(
    list(
      # x are the parameter values *input* to estimate these quantities
      x = lam,
      # error functions (note subtracted terms!)
      g = sum.g1 - 1,
      # Jacobian matrix
      j = sum.dg1.lam
    )
  )
}

### Newton's method function to solve numerically for lambda star, gamma^2_0
newt.method = function(l.init, s, r, g2_0, tol) {
  # l.init (initial guess for lambda*)
  # s, r, g2_0 are fixed values
  # tol is tolerance for routine
  
  # Initialize all values
  cur.iter = g1.lambda.calc(l.init, s, r, g2_0)
  
  while (abs(cur.iter$g) > tol) {
    # Refit while any of the error functions are larger than the tolerance
    cur.x = with(cur.iter, x - (g / j))
    cur.iter = g1.lambda.calc(cur.x[1], s, r, g2_0)
  }
  
  # Return *just the very final set of values*
  return(cur.iter$x)
  
}

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
    # Fecundity per time step
    r = w.max * (1 - s.max) / s.max,
    # Genetic fitness groups defined by 
    sig.0 = sqrt(0.4)
  ) %>%
  # Estimate lambda*
  rowwise() %>%
  mutate(lstar = newt.method(1.1, s.max, r, sig.0^2, 1e-5)) %>%
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

### Run simulations
set.seed(1310)

sim.out2 = mclapply(
  pars %>% uncount(trys.per) %>% mutate(try.no = 1:(nrow(.))) %>% split(.$try.no),
  function(pars) {
    sim(pars, theta.t = 0, init.rows = 50 * 20000) %>%
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
        s.max = pars$s.max,
        h2    = pars$h2
      )
  },
  mc.cores = 16
) %>%
  do.call(rbind, .)

write.csv(
  sim.out2,
  file = 'run_sims/out/sim_results_m1_ages_equalg20.csv',
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
