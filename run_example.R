library("diversitree")
library("picante")

source("/Users/danielesilvestro/Software/mcmc-diversitree/mcmc-SSE-lib.R")
setwd("/Users/danielesilvestro/Software/mcmc-diversitree/example_files/")


# tree file (NEXUS format, can contain multiple trees)
tree_file  = "bromelioieae_consensus.tre"

# trait table
trait_file = "traitGeoSSE.txt"

# Diversificaiton model. Available options: "musse", "geosse", "classe"
model      = "geosse" 

# state-specific taxon sampling (A, B, AB)
rho        = c(1,1,1) # 


### OTHER AVAILABLE OPTIONS
# sampling_freq = 100 
# print_freq    = 1000
# burnin        = 0
# iterations    = 50000
# constraint    = NULL # Parameters to be constrained. Introduce "lamdas", "mus" or "qs"
# nTREES        = 1    # number of trees to be analyzed
# update_freq   = c(0.2,0.4,0.99) # Given the proportion of lambdas, mus, qs and gammas
# Window size for lambda, mu, q and gamma (shape parameters for the hyperpriors)
# must be > 1
# tuning_prm = c(1.5, 1.5, 2, 2)


run_mcmc_SSE(tree_file, trait_file, model, outfile = "bromeliad_geosse.log", iterations = 10000)

# read output file
mcmc.log = "bromeliad_geosse.log"

post = read.table(mcmc.log,header=T)
plot(post$likelihood, type="l")

# plot area-specific speciation rates
hist(post$sA)


# plot speciation rates per area
par(mfrow=c(1,3))
boxplot(post[,c("sA","sB","sAB")], col = "blue", main = "speciation rates")
boxplot(post[,c("xA","xB" )], col = "red",main = "extinction rates")
boxplot(post[,c("dA","dB" )], col = "gray", main = "dispersal rates")



difference_in_dispersal_rate = post$dA - post$dB
length(difference_in_dispersal_rate[difference_in_dispersal_rate>0]) / length(difference_in_dispersal_rate)
hist(difference_in_dispersal_rate)

difference_in_extinction_rate = post$xA - post$xB
length(difference_in_extinction_rate[difference_in_extinction_rate>0]) / length(difference_in_extinction_rate)
hist(difference_in_extinction_rate)



# run on multiple trees
burnin        = 0
nTREES=5

tree_file  = "bromelioieae_100.trees"

# trait table
trait_file = "traitGeoSSE.txt"

# Diversificaiton model. Available options: "musse", "geosse", "classe"
model      = "geosse" 

run_mcmc_SSE(tree_file, trait_file, model, outfile = "bromeliad_geosse.log", iterations = 10000, nTREES=5, burnin=100)

