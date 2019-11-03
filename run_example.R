library("diversitree")
library("picante")
source("my_path/mcmc-SSE-lib.R")

# tree file (NEXUS format, can contain multiple trees)
tree_file  = "my_path/example_files/bromelioieae_consensus.tre"

# trait table
trait_file = "my_path/example_files/traitGeoSSE.txt"

# Diversificaiton model. Available options: "musse", "geosse", "classe"
model      = "geosse" 

### OTHER AVAILABLE OPTIONS
# sampling_freq = 100 
# print_freq    = 1000
# burnin        = 0
# iterations    = 50000
# rho           = NULL # state-specific taxon sampling
# constraint    = NULL # Parameters to be constrained. Introduce "lamdas", "mus" or "qs"
# nTREES        = 1    # number of trees to be analyzed
# update_freq   = c(0.2,0.4,0.99) # Given the proportion of lambdas, mus, qs and gammas
# Window size for lambda, mu, q and gamma (shape parameters for the hyperpriors)
# must be > 1
# tuning_prm = c(1.5, 1.5, 2, 2)


run_mcmc_SSE(tree_file, trait_file, model)

