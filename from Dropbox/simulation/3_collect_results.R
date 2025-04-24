main_dir <- '~/Dropbox/RESEARCH/FCTemplateICA/simulation/'
#main_dir <- '~/Dropbox (Brown)/FCTemplateICA/simulation/'
setwd(main_dir)
source('0_setup.R')

##############################################################################
### LOAD IN PREVIOUSLY CREATED RESULTS

load(file = 'GICA.RData') #ICs, Q, GICA, N

subjICs <- readRDS('subjICs.RDS')

load(file='TCs_sim.RData')
TC_sim <- TC1 

template <- readRDS('templates/template.RDS')

##############################################################################
### COMPILE RESULTS OVER SUBJECTS

n_test <- 50
inds_test <- 500 + (1:n_test)
template_mean <- newdata_xifti(GICA, template$template$mean)
template_FC_mean <- template$template$FC$psi/(template$template$FC$nu - Q - 1)
FC_mean <- cov2cor(template_FC_mean) #actually it's already a correlation


### SIMULATION PART 1 --------------------------------------------------------

algos <- c('tICA','VB1','VB2') 
num_alg <- length(algos)

result_dir <- 'results'


# IC ESTIMATES
S_true <- subjICs[,,inds_test]
S_est <- array(NA, dim=c(N, Q, n_test, num_alg + 1)) # + 1 for DR

# FC ESTIMATES 
FC_true <- abind(apply(TC_sim[601:1200,,inds_test], #truth based on hold-out second half of scan 
                       3, cor, simplify = FALSE), along=3) # list of FC matrices, abind into array
FC_true0 <- abind(apply(TC_sim[1:600,,inds_test], #truth from first half of scan (the part used for model estimation)
                       3, cor, simplify = FALSE), along=3) # list of FC matrices, abind into array
FC_est <- array(NA, dim=c(Q, Q, n_test, num_alg + 2)) # + 2 for DR and DR2
FC_LB <- FC_UB <- array(NA, dim=c(Q, Q, n_test, num_alg - 1)) # - 1 to exclude tICA

# COMPUTATION TIME, ITERATIONS
comptime <- matrix(NA, n_test, num_alg+1)
numiter <- matrix(NA, n_test, num_alg)
  
for(ii in 1:n_test){

  print(paste0('~~~~~~~~~~~~~~~~ SUBJECT ',ii,' ~~~~~~~~~~~~~~~~'))
  
  load(file=file.path(result_dir, paste0('testsubj',ii,'.RData'))) #results_ii

  #loop over algorithms
  for(aa in 2:num_alg){ 
    result_aa <- results_ii[[aa]]

    S_est[,,ii,aa] <- as.matrix(result_aa$subjICmean)
    FC_est[,,ii,aa] <- cor(result_aa$A) #VB
    FC_LB[,,ii,aa-1] <- (result_aa$FC$LB) #VB
    FC_UB[,,ii,aa-1] <- (result_aa$FC$UB) #VB
    
    numiter[ii,aa] <- result_aa$numiter
    comptime[ii,aa] <- result_aa$comptime["FC-tICA"]

    if(aa==2){
      S_est[,,ii,1] <- as.matrix(result_aa$result_tICA$subjICmean) #tICA
      S_est[,,ii,4] <- t(as.matrix(result_aa$result_DR$S)) #DR
      FC_est[,,ii,1] <- cor(result_aa$result_tICA$A) #tICA
      FC_est[,,ii,4] <- cor(result_aa$result_DR$A) #DR
      FC_est[,,ii,5] <- cor(result_aa$result_DR$A2) #DR
      numiter[ii,1] <- result_aa$result_tICA$numiter
      comptime[ii,1] <- result_aa$comptime["tICA"]
      comptime[ii,4] <- result_aa$comptime["DR"]
    }
  }
} #end loop over subjects

#save results from all test subjects
save(S_true, S_est, template_mean, 
     FC_true, FC_true0, FC_est, FC_mean, FC_LB, FC_UB, 
     comptime, numiter, 
     file = file.path(result_dir, 'test_set_results.RData'))

#save results of this analysis for part 2 (by duration analysis)
FC_est_600 <- FC_est
save(FC_est_600, file = 'results/by_duration/FC_600.RData')


  
