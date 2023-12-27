# Evolutionary rescue with overlapping generations

Models of evolutionary rescue with overlapping generations. Code associated with manuscript by Nordstrom and Melbourne (2024).

Contact: scottwatsonnordstrom & gmail : com

### Required packages

Simulations are run in `R` (version 4.3.2).

The following packages are required for simulation code to run. The version numbers used in the final runs are included, but the exact versions are likely not needed for analysis.

* `dplyr` (1.0.7) - used for piping and data frame manipulation

* `tidyr` (1.1.3) - used for tidyselect features and manipulation

* `mc2d` (0.1-22) - used for easy generation of bivariate normal distribution

The following package is used to run simulations in parallel, but is *not required* to run base simulation code.

* `parallel` (associated with R:4.3.2).

The following packages are used for plotting and analysis:

* `ggplot2` (3.4.3)

* `cowplot` (1.1.1)

### Repository structure

##### Simulation code

* `model_source/sim_model1_functions.R` - script with wrapper functions to run simulations. Note that this script includes code to run Newton's method to numerically calculate parameters (specifically, equilibrium phenotypic variances); as such it is sourced not just for running simulations but for any script where parameters are estimated.

##### Model validation

* `analytical_validation/equilibrium_conditoins.R` - script using simulations to validate analytical expressions at equilibrium conditions. Exports manuscript supporting info figures S2 - S10.

* `analytical_validation/changed_conditions.R` - script using simulations to validate analytical expressions following environmental change. Exports  manuscript supporting info figures S11 - S17.

* `analytical_validation/validation_figs` - subdirectory used to house validation figures upon creation

###### *Note about running validation scripts:*

Each of these scripts begins with running a small number of simulations; as a result, they take some time to run (under, but possibly up to, an hour). Currently, the scripts are written to run these simulations in serial, but they can easily be modified to run simulations in parallel. We ran these simulations locally (rather than on a remote server, see below), and as such there may be very slight differences in output or functioning due to minor differences in R or package versions used.

##### Running instances of simulations

* `run_sims/model1_popsize_extinctinos.R` - script that runs main batch of simulations (1000 trials for each of 27 parameter combinations). Exports files `run_sims/out/sim_results_m1_allsizes_n.csv`, `run_sims/out/sim_results_m1_disaggregated_n.csv`, `run_sims/out/sim_results_m1_sizes_r.csv`, and `run_sims_m1_survsizes_n.csv`

* `run_sims/model1_phenotypic_components.R` - script that runs smaller batch of simulations (200 trials for each of nine parameter combinations). Exports aggregated population-level mean and variance of each phenotypic component for each trial in each time step in file `run_sims/out/sim_results_m1_phtype.csv`.

* `run_sims/model1_age_structure.R` - script that runs smaller batch of simulations (200 trials for each of nine parameter combinations). Exports aggregated cohort size for each trial in each generation in each time step in file `run_sims/out/sim_results_m1_ages.csv`. Note that the parameters for these simulations, as well as the random seed used, are identical to those used in `model1_phenotypic_components.R`.

* `run_sims/out` - subdirectory containing all simulation output in .csv form

###### *Notes about running simulations*: 

We ran simulations on a server, running 16 cores in parallel using the function `mclapply()` (from the `parallel` package). The number of cores is hard-coded. To re-run locally, *you will want to change the number of cores to something that matches your device.*

Even with sixteen cores running in parallel, the main batch of simulations (`model1_popsize_extinctions.R`) takes a long time to run (several hours, but less than one day). The smaller batch of simulations (`model1_phenotypic_components.R` and `model1_age_structure.R`) 

##### Analysis of simulation output and figure generation

* `analyze_results/fig_lambda_n.R` - recreates manuscript figure 1 (note: this figure does not require any simulation data).

* `analyze_results/fig_popsize_extinctions.R` - recreates manuscript figure 2 (population size over time) and figure 3 (extinction rates over time), reading in the simulation output file `run_sims/out/sim_results_m1_allsizes_n.csv`.

* `analyze_results/fig_phenotypic_components.R` - recreates manuscript figure 4 (mean of each phenotypic component over time) and figure 6 (variance of each phenotypic component over time), reading in the simulation output file `run_sims/out/sim_results_m1_phtype.csv`.

* `analyze_results/fig_obs_age_dist.R` - recreates manuscript figure 5 (age structure for high longevity treatments) and supporting info figure S18 (age structure for medium and low longevity treatments), reading in the simulation output file `run_sims/out/sim_results_m1_ages.csv`.

* `analyze_results/fig_ex_age_dist.R` - recreates manuscript supporting info figure S1 (note: this figure does not require any simulation data).

* `analyze_results/generation_time_scaling.R` - recreates manuscript supporting info figure S19, which shows mean phenotypic components over time, but with time scaled by the generation time of each treatment, reading in data from `run_sims/out/sim_results_m1_phtype.csv`. This script also features code to estimate generation time. 

* `analyze_results/figs_out` - subdirectory used to house figures upon creation

###### *Note about recreating figures*

Script for analyses were run on a local machine, i.e., *not on the same machine used to run simulations.* As such there may be very slight discrepancies in R and package versions used between the two devices. While we imagine this is unlikely to change results, we can not guarantee that results will be identical to those we found for this reason.

