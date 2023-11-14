##### Source code for main simulation model

##### Load packages
library(dplyr)
library(tidyr)
library(mc2d) # required for bivariate normal distribution

### Suppress output from summarise()
options(dplyr.summarise.inform = FALSE)

##### Aux functions

### Wrapper for adding data to data frame
# (a pre-allocated data frame to improve performance)
dim.add = function(df, rows, addition) {
  if (sum(is.na(df$i)) < nrow(addition)) {
    df = df %>%
      rbind(df %>%
              sample_n(size = rows) %>%
              mutate_all(function(x) NA))
  }
  if (nrow(addition)) {
    df[which.max(is.na(df$i)) + 0:(nrow(addition)-1),] = addition
  }
  return(df)
}

### Functions for getting gamma^2_0 for a given s, r
# that define the stage distribution

g1.fun = function(g2_0, s, r) {
  # g2_0 is gamma^2_0 = sigma^2_0 / w^2
  # s and r are max survival and fecundity
  sumo = 0
  for (k in 0:1e5) sumo = sumo + (s^k / sqrt(1 + k*g2_0))
  return((r / (r+1)) * sumo - 1)
}

g1.prm = function(g2_0, s, r) {
  # derivative of g1 with respect go gamma^2_0
  # (note: not with respect to gamma)
  sumo = 0
  for (k in 0:1e5) sumo = sumo - (s^k * (k/2) / sqrt(1 + k*g2_0)^3)
  return(r / (r+1) * sumo)
}

newt.method.g1 = function(x0, tol, s, r) {
  # Wrapper for Newton's method
  xold = x0
  while (abs(g1.fun(xold, s, r)) > tol) xold = xold - (g1.fun(xold, s, r) / g1.prm(xold, s, r))
  return(xold)
}

gamma.calc = function(g2_0, s, r) {
  # function for getting overall population phenotypic variance
  # (in terms of gamma)
  sumo = 0
  for (k in 0:1e5) sumo = sumo + ((s^k / sqrt(1 + k*g2_0)) * (g2_0 / (1 + k*g2_0)))
  return(r / (r+1) * sumo)
}

gamma.a.calc = function(g2a0, s, r, g2e0) {
  # function for getting overall population genotypic variance
  # (in terms of gamma)
  # needed for getting mutation rate (to hold population-level variance constant)
  sumo = 0
  for (k in 0:1e5) {
    sumo = sumo + ((s^k / sqrt(1 + k*(g2a0 + g2e0))) * ((g2a0*(1 + k*g2e0)) / (1 + k*(g2e0 + g2a0))))
  }
  return(r / (r+1) * sumo)
}

### Wrapper for initializing sim

init.sim = function(params, theta0) {
  
  # initial population size
  size0 = params$n.pop0
  # optimal phenotypic value
  s.max = params$s.max 
  # mean fecundity per individual
  r     = params$r
  # selection strength (this seems like a good value to use)
  wfitn = params$wfitn
  # genetic variance (should be approx. static)
  sig.a = params$sig.a
  # standard dev. of environmental phenoytpic noise
  sig.e = params$sig.e
  # initial genotype
  gbar0 = ifelse(any(grepl('gbar0', names(params))), params$gbar0, 0)
  # density dependence strength
  alpha = ifelse(any(grepl('alpha', names(params))), params$alpha, 0)
  # ceiling-like carrying capacity term
  kceil = ifelse(any(grepl('ceil' , names(params))), params$kceil, Inf)
  # equilibrium population growth rate 
  lstar = ifelse(any(grepl('lstar', names(params))), params$lstar, 1)
  
  # Get age distribution
  p.age = data.frame(age = 0:500) %>%
    mutate(p.age = (r / (1 + r)) * 
             (s.max / lstar)^age * 
             sqrt(wfitn^2 / (wfitn^2 + age*(sig.a^2 + sig.e^2)))
    )
  
  # Get phenotpyic age variance/covariance matrices
  vcv.a = vector('list', length = nrow(p.age))
  
  for (k in 1:length(vcv.a)) {
    vcvk = c(sig.a^2 * (1 + (k-1)*sig.e^2) / (1 + (k-1)*(sig.a^2 + sig.e^2)),
             -(k-1) * sig.a^2 * sig.e^2    / (1 + (k-1)*(sig.a^2 + sig.e^2)),
             -(k-1) * sig.a^2 * sig.e^2    / (1 + (k-1)*(sig.a^2 + sig.e^2)),
             sig.e^2 * (1 + (k-1)*sig.a^2) / (1 + (k-1)*(sig.a^2 + sig.e^2))
    )
    vcv.a[[k]] = matrix(vcvk, nrow = 1)
  }
  
  popn = data.frame(
    # Unique identifier
    i   = 1:size0,
    # Time step
    t   = 0,
    # Age
    age = sample(size = size0, x = p.age$age, prob = p.age$p.age, replace = TRUE),
    # Sex (TRUE = female)
    fem = as.logical(sample(0:1, size0, replace = TRUE))
  ) %>%
    cbind(
      rmultinormal(size0, mean = c(gbar0, 0), sigma = do.call(rbind, vcv.a[.$age + 1]))
    ) %>%
    rename(b_i = `1`, e_i = `2`) %>%
    mutate(
      # Phenotype
      z_i = b_i + e_i,
      # Survival
      s_i = s.max * exp(-(z_i - theta0)^2 / (2*wfitn^2)) * exp(-alpha * size0),
      # Offspring (0 for males, Poisson draw for females)
      r_i = rpois(size0, lambda = ifelse(fem, 2 * r, 0)),
      # Phenotypic optimum in this time step
      theta_t = theta0
    ) %>%
    select(-e_i)
  
  # If initial size is above the specified ceiling, truncate
  if (nrow(popn) > kceil) { 
    popn = popn %>% 
      slice_sample(kceil, replace = FALSE) 
    # Sort rows to order individuals by index
    arrange(i)
  }
  
  return(popn)
  
}

### Function for porpagating one time step of simulated populations forward
# NOTE: here, survival (and selection) are first
# under this model, replacement rate of individual i should be s_i (1 + r_i)
propagate.sim = function(popn, params, theta, cur.i) {
  
  # Pull out parameters of interest
  
  # Maximum survival
  s.max = params$s.max
  # mean fecundity per individual
  r     = params$r
  # selection strength (this seems like a good value to use)
  wfitn = params$wfitn
  # genetic variance (should be approx. static)
  sig.a = params$sig.a
  # standard dev. of environmental phenoytpic noise
  sig.e = params$sig.e
  # mutation rate per capita
  mu    = ifelse(any(grepl('mu', names(params))), params$mu, 0)
  # standard dev. of mutations
  sig.m = ifelse(any(grepl('sig.m', names(params))), params$sig.m, 0)
  # Strength of density dependence
  alpha = ifelse(any(grepl('alpha', names(params))), params$alpha, 0) 
  # ceiling-like carrying capacity term
  kceil = ifelse(any(grepl('ceil' , names(params))), params$kceil, Inf)
  
  # Current time step (useful for adding time information to finished data frame)
  cur.t = max(popn$t)
  
  # Survival (of all individuals)
  popn.surv = popn %>% 
    # recalculate s_i based on current theta value
    mutate(theta_t = theta) %>%
    mutate(s_i = s.max * exp(-(z_i - theta_t)^2 / (2*wfitn^2)) * exp(-alpha * nrow(.))) %>%
    # simulate survival step using Bernoulli draw for each individual
    filter(as.logical(rbinom(n = nrow(.), size = 1, prob = s_i))) %>%
    # Handle parents - iterating forward age, current time step
    mutate(age = age + 1, t = cur.t + 1)  %>%
    # get new number of offspring per parent
    mutate(r_i = rpois(nrow(.), lambda = ifelse(fem, 2 * r, 0)))
      
    
  # Perform mating to get offspring (if there are fecund females)
  if (sum(popn.surv$r_i) & sum(popn.surv$fem) & sum(!popn.surv$fem)) {
    
    offspring = cbind(
      # Mom info
      popn.surv %>%
        filter(fem) %>%
        select(-c(i, t, age, fem, z_i, s_i, theta_t)) %>%
        rename(mom_b_i = b_i),
      # Dad info
      popn.surv %>%
        sample_n(size = sum(fem),
                 weight = as.numeric(!fem) / sum(as.numeric(!fem)),
                 replace = TRUE) %>%
        select(b_i) %>%
        rename(dad_b_i = b_i)
    ) %>%
      # Get midparent value for each mating pair
      mutate(midp_b_i = (mom_b_i + dad_b_i) / 2) %>%
      # Duplicate midparent value according to r_i, number of offspring per pair
      uncount(r_i) %>%
      mutate(i = max(popn$i) + 1:nrow(.),           # assign ID number
             t = cur.t + 1,                         # add time step
             age = 0,                               # all offspring are age 0
             fem = as.logical(sample(0:1, nrow(.), replace = TRUE)), # assign sex
             b_i = rnorm(nrow(.), midp_b_i, sqrt(var(popn.surv$b_i) / 2)) +
                    rbinom(nrow(.), 1, mu) * rnorm(nrow(.), 0, sig.m),
                                                    # assign breeding value (mean is midparent, sd is sqrt(additive var))
             z_i = rnorm(nrow(.), b_i, sig.e),      # assign phenotype (mean is breding value, sd is sd of env. variance)
             s_i = s.max * exp(-(z_i - theta)^2 / (2*wfitn^2)) * exp(-alpha * nrow(.)), # determine survival probability
             r_i = rpois(nrow(.), lambda = ifelse(fem, 2 * r, 0)), # re-draw number of offspring per mother
             theta_t = theta) %>%                   # phenotypic optimum in this time step
      select(i, t, age, fem, b_i, z_i, s_i, r_i, theta_t)
    
  } else {
    # If sum(r_i) is zero, there are zero offspring from mating this time step
    offspring = popn.surv %>% sample_n(size = 0)
  }
  
  # Combine offspring and parents
  popn.out = rbind(popn.surv, offspring) 
  
  # If offspring + parents is above the ceiling carrying capacity, truncate
  if (nrow(popn.out) > kceil) {
    popn.out = popn.out %>% 
      # If population is above the ceiling, sample just the ceiling
      slice_sample(n = kceil, replace = FALSE) %>%
      # Sort rows to order individuals by index
      arrange(i)
  }
  
  # Return
  return(popn.out)
  
}

### Wrapper function for simulating a population

sim = function(params, theta.t, init.rows, init.popn = NULL) {
  
  ### Declaring variables
  # (must this be done here?)
  
  # length of simulations
  timesteps = params$timesteps
  # initial population size
  size0 = params$n.pop0
  # optimal phenotypic value
  s.max = params$s.max 
  # mean fecundity per individual
  r     = params$r
  # selection strength (this seems like a good value to use)
  wfitn = params$wfitn
  # genetic variance (should be approx. static)
  sig.a = params$sig.a
  # standard dev. of environmental phenoytpic noise
  sig.e = params$sig.e
  # initial genotype
  gbar0 = ifelse(any(grepl('gbar0', names(params))), params$gbar0, 0)
  # density dependence strength
  alpha = ifelse(any(grepl('alpha', names(params))), params$alpha, 0)
  # ceiling-like carrying capacity term
  kceil = ifelse(any(grepl('ceil' , names(params))), params$kceil, Inf)
  
  ### Initialize data frame
  all.data = data.frame(
    i   = rep(NA, init.rows),
    t   = rep(NA, init.rows),
    age = rep(NA, init.rows),
    fem = rep(NA, init.rows),
    b_i = rep(NA, init.rows),
    z_i = rep(NA, init.rows),
    s_i = rep(NA, init.rows),
    r_i = rep(NA, init.rows),
    theta_t = rep(NA, init.rows)
  )
  
  ### Check lengths of theta parameter
  if (!length(theta.t) %in% c(timesteps+1, 1)) {
    stop("Length of theta vector does not match other inputs")
  } else if (length(theta.t) == 1) {
    theta.t = rep(theta.t, timesteps+1)
  }
  
  ### Initialize population
  if (!is.null(init.popn)) {
    init.popn = init.popn %>%
      mutate(theta_t = theta.t[1]) %>%
      mutate(s_i = s.max * exp(-(z_i - theta_t)^2 / (2*wfitn^2)) * exp(-alpha * nrow(.)))
  } else {
    init.popn = init.sim(params = params, theta0 = theta.t[1])
  }
  
  ### Add initial population to data frame
  all.data = dim.add(df = all.data, rows = init.rows, addition = init.popn)
  
  ### Iterate simulation
  prev.t = init.popn
  cur.i  = max(init.popn$i)
  
  for (tt in 1:timesteps) {
    if (nrow(prev.t)) {
      # Propagate sim (one time step)
      popn = propagate.sim(popn = prev.t, params = params, theta = theta.t[tt+1], cur.i = cur.i)
      # Add to data frame
      all.data = dim.add(df = all.data, rows = init.rows, addition = popn)
      # Update "previous time step"
      prev.t = popn
      # Update current ID (for assigning future IDs)
      cur.i = max(all.data$i)
    }
  }
  
  ### Remove empty rows at the end
  all.data = all.data %>% filter(!is.na(i))
  
  return(all.data)
  
}

### Function for porpagating one time step of simulated populations forward
# NOTE: this one has reproduction (r) first, then survival
# this may affect model behavior slightly, including analytical solutions
propagate.sim.r.first = function(popn, params, theta) {
  
  # Pull out parameters of interest
  
  # Maximum survival
  s.max = params$s.max
  # mean fecundity per individual
  r     = params$r
  # selection strength (this seems like a good value to use)
  wfitn = params$wfitn
  # genetic variance (should be approx. static)
  sig.a = params$sig.a
  # standard dev. of environmental phenoytpic noise
  sig.e = params$sig.e
  # Strength of density dependence
  alpha = ifelse(any(grepl('alpha', names(params))), params$alpha, 0) 
  
  # Only iterate if there are both male and females available for mating
  if (sum(popn$fem) & sum(!popn$fem)) {
    
    # Current time step (useful for adding time information to finished data frame)
    cur.t = max(popn$t)
    
    # Perform mating to get offspring (if there are fecund females)
    if (sum(popn$r_i)) {
      
      offspring = cbind(
        # Mom info
        popn %>%
          filter(fem) %>%
          select(-c(i, t, age, fem, z_i, s_i, theta_t)) %>%
          rename(mom_b_i = b_i),
        # Dad info
        popn %>%
          sample_n(size = sum(fem),
                   weight = as.numeric(!fem) / sum(as.numeric(!fem)),
                   replace = TRUE) %>%
          select(b_i) %>%
          rename(dad_b_i = b_i)
      ) %>%
        # Get midparent value for each mating pair
        mutate(midp_b_i = (mom_b_i + dad_b_i) / 2) %>%
        # Duplicate midparent value according to r_i, number of offspring per pair
        uncount(r_i) %>%
        mutate(i = max(popn$i) + 1:nrow(.),           # assign ID number
               t = cur.t + 1,                         # add time step
               age = 1,                               # all offspring are age 1
               fem = as.logical(sample(0:1, nrow(.), replace = TRUE)), # assign sex
               b_i = rnorm(nrow(.), midp_b_i, sig.a) +
                      as.numeric(rbinom(nrow(.), 1, mu) * rnorm(nrow(.), 0, sig.m)),
                                                      # assign breeding value (mean is midparent, sd is sqrt(additive var))
               z_i = rnorm(nrow(.), b_i, sig.e),      # assign phenotype (mean is breding value, sd is sd of env. variance)
               s_i = s.max * exp(-(z_i - theta)^2 / (2*wfitn^2)) * exp(-alpha * nrow(.)), # determine survival probability
               r_i = rpois(nrow(.), lambda = ifelse(fem, 2 * r, 0)), # re-draw number of offspring per mother
               theta_t = theta) %>%                   # phenotypic optimum in this time step
        select(i, t, age, fem, b_i, z_i, s_i, r_i, theta_t)
      
    } else {
      # If sum(r_i) is zero, there are zero offspring from mating this time step
      offspring = popn %>% sample_n(size = 0)
    }
    
    # Survival (of all individuals - includes those just born)
    
    popn.out = popn %>% 
      # Handle parents - iterating forward age, current time step
      mutate(age = age + 1, t = cur.t + 1)  %>%
      # get new number of offspring per parent
      mutate(r_i = rpois(nrow(.), lambda = ifelse(fem, 2 * r, 0))) %>%
      # recalculate s_i based on current theta value
      mutate(theta_t = theta) %>%
      mutate(s_i = s.max * exp(-(z_i - theta_t)^2 / (2*wfitn^2)) * exp(-alpha * nrow(.))) %>%
      # combine with offspring
      rbind(offspring) %>%
      # simulate survival step using Bernoulli draw for each individual
      filter(as.logical(rbinom(n = nrow(.), size = 1, prob = s_i)))
    
    
  } else {
    # If the population passed in doesn't have individuals of both sex, return and empty DF
    popn.out = popn %>% sample_n(size = 0)
  }
  
  # Return
  return(popn.out)
  
}
