library("diversitree")
library("picante")
library("phytools")

comb_areas <- function(data_SPGC, columns_areaA, columns_areaB){
  x1 = rowSums(data_SPGC[, columns_areaA])
  x1[x1>0] <- 1
  x2 = rowSums(data_SPGC[, columns_areaB])
  x2[x2>0] <- 2	
  out <- x1 + x2	
  names(out) <- NULL
  out <- data.frame(label = rownames(data_SPGC), trait = out)
  out <- out[out$trait > 0, ]
  out$trait <- out$trait - 1	
  return(out)
}




run_mcmc_SSE <- function( tree, trait.data, model = "musse", rho = NULL, 
                          sampling_freq = 100, print_freq = 1000, iterations = 50000, 
                          constraint = NULL, nTREES = 1, burnin = 0, outfile = "",
                          update_freq = c(0.2,0.4,0.99), tuning_prm = c(1.5,1.5,2,2)){
  
  # # Reading Input Files
  # tree=read.nexus(file=tree_file)
  # # MLS = change from read.tree
  #
  #
  if (class(tree)=="multiPhylo"){current_tree=tree[[1]]
  }else{current_tree=tree}
  
  # states_frame <- read.table(file=trait_file, h=T, row.names=1)
  states_frame <- trait.data
  matched_tree <- match.phylo.data(current_tree,states_frame) # Drop tips with no info
  matched_tree$data <- as.data.frame(matched_tree$data)
  states_vector <- as.data.frame(matched_tree$data)[,1] # vector of states
  names(states_vector) <- row.names(matched_tree$data) # adding names to each state
  no_states <- length(unique(states_vector))
  if(is.null(rho)) # by default complete sampling
  {
    rhos <- rep(1,no_states)
  } else
  {
    rhos <- rho
  }
  
  current_tree <- matched_tree$phy
  
  #####################################
  # Defining Variables and Objects
  print(model)
  if (model == "musse"){
    initial_parameters <- starting.point.musse(current_tree, no_states, q.div=5, yule=FALSE)   # Starting conditions
  } else if (model == "geosse") {
    initial_parameters <- starting.point.geosse(current_tree, eps=0.5)
  } else if (model == "classe") {
    initial_parameters <- starting.point.classe(current_tree, no_states, eps=0.5)
  }
  
  acc=0
  
  d <- NULL
  # Window size for lambda
  d[c(1:no_states)] <- tuning_prm[1]
  # Window size for mu
  d[c((no_states+1):(no_states*2))] <- tuning_prm[2]
  # Window size for q
  d[c((no_states*2+1):(length(initial_parameters)))] <- tuning_prm[3]
  # Window size for gamma_exp_rate_hp
  d_exp_rate_hp <- rep(tuning_prm[4],3)
  #print(d_exp_rate_hp)
  
  
  ######################################
  # Defining Functions
  
  # Sliding window
  update_parameters_sw<- function(accepted_state_values, indexes, d, M) {
    i = accepted_state_values[indexes]
    ii = NULL
    ii= abs(i+(runif(length(i),0,1)-.5)*d[indexes])
    #for (j in 1:length(i)){
    #  ii[j] = abs(i[j]+(runif(1,0,1)-.5)*d[indexes[j]])
    #}
    if (max(ii) > M) {ii=abs((M-(ii-M)))}
    accepted_state_values[indexes]=ii
    return(accepted_state_values)
  }
  
  # Updating parameters through multiplier proposal
  update_parameters_mp <- function(accepted_state_values, indexes, d) {
    u = runif(length(indexes),0,1)
    lambda =  2*log(d[indexes])
    m = exp(lambda*(u-0.5))
    accepted_state_values[indexes] = accepted_state_values[indexes]*m
    return(list(accepted_state_values,sum(log(m))))
  }
  
  if (model == "musse"){
    
    indexes_lamb= c(1:no_states)
    indexes_mu= c((no_states+1):(no_states*2))
    indexes_q= c((no_states*2+1):(length(initial_parameters)) )
  } else if (model == "geosse") {
    
    indexes_lamb= c(1:no_states)
    indexes_mu= c((no_states+1):(no_states+2))
    indexes_q= c((no_states+3):(length(initial_parameters)) )
  } else if (model == "classe") {
    x = no_states*no_states*(no_states+1)/2
    indexes_lamb= c(1:x)
    indexes_mu= c((x+1):(x+no_states))
    indexes_q= c((x+no_states+1):(length(initial_parameters)) )
  } else if (model == "quasse") {
    x = no_states*no_states*(no_states+1)/2
    indexes_lamb= c(1:x)
    indexes_mu= c((x+1):(x+no_states))
    indexes_q= c((x+no_states+1):(length(initial_parameters)) )}
  
  
  
  # gibbs sampler exp rate parameter with Gamma prior
  sample_rate_exp <-function(x,hp_G_alpha=2,hp_G_beta=2){
    g_shape= hp_G_alpha + length(x)
    g_rate=  hp_G_beta  + sum(x)
    lam = rgamma(1, shape= g_shape, rate= g_rate)
    return (lam)		    	
  }
  
  
  
  ##################################
  # Main Program
  if (outfile==""){
    output_file = paste(tree_file, "_", model, "_mcmcEXP.log", sep="")		
  }else{
    output_file = outfile
  }
  
  cat(paste("Writing to File: ", output_file, "\n", sep=""))
  
  
  real_it=0
  for (J in 1:nTREES){
    if (class(tree)=="multiPhylo"){current_tree=tree[[J]]}
    else{current_tree=tree}
    
    if (!is.ultrametric(current_tree)){
        current_tree <- force.ultrametric(current_tree)
    }
    
    
    states_frame <- trait.data
    matched_tree <- match.phylo.data(current_tree,states_frame) # Drop tips with no info
    matched_tree$data <- as.data.frame(matched_tree$data)
    states_vector <- as.data.frame(matched_tree$data)[,1] # vector of states
    names(states_vector) <- row.names(matched_tree$data) # adding names to each state
    current_tree <- matched_tree$phy
    current_exp_rate_hp <- rep(1,3) # starting values for hyperpriors parameters
    
    if (model == "musse") {
      LIKELIHOOD <- make.musse(current_tree, states_vector, no_states, strict=T, sampling.f=rhos)
      initial_parameters <- starting.point.musse(current_tree, no_states, q.div=5, yule=FALSE)   # Starting conditions
      
    } else if (model == "classe") { 
      LIKELIHOOD <- make.classe(current_tree, states_vector, no_states, strict=T, sampling.f=rhos)
      initial_parameters <- starting.point.classe(current_tree, no_states)   # Starting conditions
    } else if (model == "geosse") {
      LIKELIHOOD <- make.geosse(current_tree, states_vector, strict=T, sampling.f=rhos)
    }
    
    initial_parameters <- initial_parameters + 0.1
    
    current_state=initial_parameters
    accepted_state=current_state
    accepted_lik=-Inf
    accepted_prior=-Inf
    
    if (real_it==0) {
      cat(c("Iteration", "likelihood", "prior", "acceptance", "tree", names(accepted_state), "g_lamb", "g_mu", "g_q"), "\n", append=FALSE, file=output_file, sep="\t")
    }
    
    acc=0
    
    for (it in 1:iterations){
      
      # Updating the parameters
      #
      flag <- runif(1,0,1)
      
      # if (flag < update_freq[1]) {
      #     ind = sample(indexes_lamb,1)
      #     temp <- update_parameters_mp(accepted_state, ind, d)
      # } else if (flag < update_freq[2]) {
      #     ind = sample(indexes_mu,1)
      #     temp <- update_parameters_mp(accepted_state, ind, d)
      # } else if (flag < update_freq[3]) {
      #     ind = sample(indexes_q,1)
      #     temp <- update_parameters_mp(accepted_state, ind, d)
      
      if (flag < update_freq[3] | it==1) {
        ind = sample(indexes_lamb,1)
        temp1 <- update_parameters_mp(accepted_state, ind, d)
        ind = sample(indexes_mu,1)
        temp2 <- update_parameters_mp(temp1[[1]], ind, d)
        ind = sample(indexes_q,1)
        temp3 <- update_parameters_mp(temp2[[1]], ind, d)
        
        current_state <- temp3[[1]]
        hasting <- temp1[[2]]+temp2[[2]]+temp3[[2]]
        gibbs_s=0
        
        
        # Lik Calculation
        current_lik <- try(LIKELIHOOD(c(current_state)))
        if (is.na(current_lik) | (class(current_lik) == "try-error" )) {current_lik = -Inf}
        
      } else {
        # GIBBS SAMPLER
        current_exp_rate_hp[1] = sample_rate_exp(current_state[indexes_lamb])
        current_exp_rate_hp[2] = sample_rate_exp(current_state[indexes_mu]  )
        current_exp_rate_hp[3] = sample_rate_exp(current_state[indexes_q]   )
        current_lik=accepted_lik
        current_state=accepted_state
        hasting = 0
        gibbs_s=1
      }
      
      # Prior calculation
      #
      current_prior <- sum(dexp(current_state[indexes_lamb],rate=current_exp_rate_hp[1],log=TRUE)) +
        sum(dexp(current_state[indexes_mu]  ,rate=current_exp_rate_hp[2],log=TRUE)) +
        sum(dexp(current_state[indexes_q]   ,rate=current_exp_rate_hp[3],log=TRUE))
      
      
      # Acceptance Conditions
      # print(c(current_lik, accepted_lik))
      if ( (current_lik-accepted_lik) + (current_prior-accepted_prior) + hasting >= log(runif(1,0,1)) | it==0 | gibbs_s==1){
        accepted_lik=current_lik
        accepted_prior=current_prior
        accepted_state=current_state
        acc = acc + 1
      }
      
      if (it==0) {acc=1}
      
      # Appending % of accepted states
      if ( (it %% sampling_freq == 0) & (it >= burnin) ) {
        cat(c(real_it, accepted_lik, accepted_prior, acc/it, J, accepted_state, current_exp_rate_hp,"\n"),sep="\t", append=TRUE, file=output_file)
        real_it=real_it+sampling_freq
      }
      if ( (it %% print_freq == 0)) {
        cat(sprintf("%s\t", round(c(it, accepted_lik, accepted_prior,acc/it),3)), "\n", sep="")
      }
      
      
    }
    
  }
}