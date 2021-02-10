source("simulate_cds.R")
simulate_synonymous_divergence(N=1000, basename="highrate", rate=0.1)
simulate_synonymous_divergence(N=1000, basename="mediumrate", rate=0.01)
simulate_synonymous_divergence(N=1000, basename="lowrate", rate=0.001)
