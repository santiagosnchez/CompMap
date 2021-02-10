source("simulate_cds.R")
# N, number of genes
# basename, base name for simulated dataset
# rate, synonymous site divergence between two sequences
simulate_synonymous_divergence(N=1000, basename="highrate", rate=0.1)
simulate_synonymous_divergence(N=1000, basename="mediumrate", rate=0.01)
simulate_synonymous_divergence(N=1000, basename="lowrate", rate=0.001)
