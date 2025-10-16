##########
# Simulated experiment
# Here, run with newborn cohort phenotypic variance equal (rather than
# population-level phenotypic variance)
# Population size and rates of extionction
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

### Write two mroe wrapper functions

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
trys.per = 1000

# Parameters

pars = expand.grid(
  # Longevity groups defined by s.max
  s.max = c(0.1, 0.5, 0.9),
  # Genetic fitness groups defined by 
  sig.0 = sqrt(c(0.1, 0.25, 0.4))
) %>%
  mutate(
    # Max lifetime fitness
    w.max = 3,
    # Fecundity per time step
    r = w.max * (1 - s.max) / s.max
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
    timesteps = 100,
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
set.seed(1310)

sim.out = mclapply(
  pars %>% uncount(trys.per) %>% mutate(try.no = 1:(nrow(.))) %>% split(.$try.no),
  function(pars) {
    sim(pars, theta.t = 0, init.rows = 100 * 20000) %>%
      group_by(t) %>%
      summarise(n = n())  %>%
      mutate(
        trial = pars$try.no, 
        s.max = pars$s.max,
        var.z = pars$sig.0^2,
        h2    = pars$h2
      )
  },
  # Number of cores to run on in parallel
  mc.cores = 6
) %>%
  do.call(rbind, .)

# Get mean population size for all time steps
sim.n.all = sim.out %>%
  # Need to merge these in to include extinct populations as size zero
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
    s.max = s.max[1],
    var.z = var.z[1],
    h2 =    h2[1]
  ) %>%
  group_by(s.max, var.z, h2, t, trial) %>%
  # (taking sums to include zeros)
  summarise(n = sum(n)) %>%
  group_by(s.max, var.z, h2, t) %>%
  # Get mean and variance of population groups
  summarise(
    nbar = mean(n),
    nvar = var(n),
    # psrv is proportion of populations surviving
    psrv = mean(n > 0),
    # Should be equal to ntrials
    nn   = n()
  )

# Get 20 trials per parameter combination
sim.disagg = sim.out %>% filter((trial %% trys.per) < 20)

# Export population size over time
write.csv(
  sim.n.all,
  row.names = FALSE,
  file = 'run_sims/out/sim_results_m1_allsizes_n_equalg20.csv'
)

write.csv(
  sim.disagg,
  row.names = FALSE,
  file = 'run_sims/out/sim_results_m1_disaggregated_n_equalg20.csv'
)

# Export session info
writeLines(capture.output(sessionInfo()), 'run_sims/out/sessioninfo.txt')

### SesssionInfo (30 Apr 2025)

