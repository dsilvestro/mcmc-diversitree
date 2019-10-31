## mcmc-diversitree
This repository includes Bayesian implementations of State Speciation and Extinction models, including MuSSE, GeoSSE and ClaSSE.
The use of (half-Cauchy or Exponential) hyper-priors on the rate parameters allows to control for over-parameterization. 
The `mcmc-SSE.R` script uses exponential priors on the rate parameters with a gamma hyper-prior.

### Main commands
The script requires a tree file (NEXUS format), a text file with the trait data (see [example files](https://github.com/dsilvestro/mcmc-diversitree/tree/master/example_files)), and specifying a model. The models currently available are "musse", "classe", "geosse".  

`RScript mcmc-SSE.R example_files/bromelioieae_consensus.tre example_files/traitGeoSSE.txt geosse`

The taxon sampling (fraction of sampled species out of the total) is provided for each character state using the flag `--rho`. For example `--rho "0.5 0.4 1" ` in a geosse model specifies a 50%, 40% sampling for species in area 1 and 2, respectively, and a complete sampling for widespread species.

`RScript mcmc-SSE.R example_files/bromelioieae_consensus.tre example_files/traitGeoSSE.txt geosse --rho "0.5 0.4 1" `


Additional options are:

`--i`: number of MCMC iterations  
`--s`: sampling frequency  
`--p`: print frequency





## Requirements
R libraries: `ape`, `optparse`, `picante`, `diversitree` all available on [CRAN](https://cran.r-project.org).


## References
### `mcmc-diversitree`:  
[Silvestro et al. 2014 Evolution](http://onlinelibrary.wiley.com/doi/10.1111/evo.12236/abstract)  
[Burin et al. 2016 Nature Comm](http://www.nature.com/ncomms/2016/160407/ncomms11250/full/ncomms11250.html)

### `diversitree` library:  
[FitzJohn et al. 2008 Syst Biol](http://sysbio.oxfordjournals.org/content/58/6/595)  
[GitHub repository](https://github.com/richfitz/diversitree)

