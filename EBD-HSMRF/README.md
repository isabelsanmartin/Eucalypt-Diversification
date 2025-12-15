HSMRF-EBD (Horseshoe Markov Random Field EBD) analysis

This autocorrelated Brownian model, implemented in RevBayes (Höhna et al. 2016) assumes the median rate within an interval is centred on that of the previous interval, making it well-suited to capture gradual changes. To mitigate excessive autocorrelation, a global scale hyperparameter  allows rare, extreme rate shifts. 

Preliminary runs with tutorial defaults (https://revbahes.github.io/tutorials/divrate.html) underfit extreme variation detected by 496 COMET, so we multiplied the global scale parameter for speciation and extinction rates by 100, to allow for greater rate variation. Four independent MCMC chains of 10.000 generations each were run. Convergence was assessed using the R package Convenience68, with ESS diagnostics for precision and Kolmogorov-Smirnov (KS) tests 500 for reproducibility. 

Scripts to reproduce the HSMRF-EBD analysis:"mcmc_EBD_HSMRF_env.Rev" (main scripts to source in RevBayes); "Convenience_Convergence.R" (script to assess convergence); "plot_EBD_env_corr_mod.R" (R script to plot Bayes Factor statistical test for significant correlations between speciation/extinction and global temperature); "plot_EBD_env_mod.R" (R script to plot speciation, extinction, diversification and relative extinction rates over time, with global temperature curve overlaid on top).
