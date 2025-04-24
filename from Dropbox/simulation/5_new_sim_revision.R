library(fMRItools)

main_dir <- '~/Dropbox/RESEARCH/FCTemplateICA/simulation/'
setwd(main_dir)
source('0_setup.R')

first_run <- FALSE

####################################################################
### GENERATE SCANS WITH LONG DURATION
####################################################################

ICs <- c(1,3,4,2,16)
Q <- length(ICs)

# GICA_fname <- '../data/melodic_IC_25.dscalar.nii' #avoid cloud-based files!
# GICA <- read_cifti(GICA_fname, brainstructures = c('left','right'), resamp_res = 10000)
# GICA <- newdata_xifti(GICA, as.matrix(GICA)[,ICs])
# save(GICA, file='GICA_longsim.RData')
load(file = 'GICA_longsim.RData') 

template_mean <- as.matrix(GICA/20)
template_mean[template_mean > 1] <- 1
template_mean <- newdata_xifti(GICA, template_mean)
template_var <- (template_mean/1.5)^2
N <- nrow(template_mean)


###################################################################################
# GENERATE SUBJECT ICS (need high-resolution since T is now big?)

FWHM <- 8 #smoothing kernel
n <- 500 + 50 #500 training subjects (use to re-estimate template) + 50 test subjects

if(first_run){
  subjICs <- array(dim=c(N, Q, n))
  tmean_mat <- as.matrix(template_mean)
  tvar_mat <- as.matrix(template_var)
  for(ii in 1:n){
    print(ii)
    for(q in 1:Q){
      dev_iq <- rnorm(N, rep(0,N), sqrt(tvar_mat[,q])) #generate subject effects
      dev_iq <- smooth_xifti(newdata_xifti(GICA, dev_iq), surf_FWHM = FWHM)
      subjICs[,q,ii] <- scale(tmean_mat[,q] + as.matrix(dev_iq), scale=FALSE)
    }
  }
  saveRDS(subjICs, 'subjICs_longsim.RDS')
}

subjICs <- readRDS('subjICs_longsim.RDS')
template_mean_est <- apply(subjICs[,,1:500], c(1,2), mean)
template_var_est <- apply(subjICs[,,1:500], c(1,2), var)


###################################################################################
# RE-ESTIMATE TEMPLATES

if(first_run){
  
  #generate and write out simulated fMRI data for each TRAINING subject
  BOLD_fnames <- paste0('data_longsim/subj',1:500,'.dtseries.nii')
  ntime <- 1200 #length of fMRI timeseries

  #first, determine signal variance for SNR to fix error variance
  var_TC <- 1 #pre-set based on scaling
  tvar <- rep(0,Q)
  #identify intensity of peak voxels for each IC (top 1% of each IC)
  peak_vals <- apply(template_mean_est, 2, function(x) {
    mean(x[x > quantile(x, 0.99)])
  })
  #signal variance = variance of signal time courses scaled by peak IC intensity, averaged over ICs
  sd_sig <- sqrt(mean(var_TC*peak_vals)) #0.93 approx
  sd_err <- 2*sd_sig #1.86 approx
  
  #time courses from original simulation
  load(file='TCs_sim.RData')
  TC_sim <- TC1 
  
  set.seed(734859) #set seed before generating seeds for each subject
  seeds <- runif(n, min=0, max=1000000) #generate a seed for each subject's noise
  xifti0 <- ciftiTools:::convert_to_dtseries(newdata_xifti(GICA, matrix(0, N, ntime))) #template xifti object
  for(ii in 1:n){
    print(ii)
    TC_ii <- TC_sim[,,ii]
    set.seed(seeds[ii])
    errs_ii <- matrix(rnorm(N*ntime, mean=0, sd=sd_err), nrow=N, ncol=ntime) #random noise
    data_ii <- subjICs[,,ii] %*% t(TC_ii) + errs_ii
    xifti_ii <- newdata_xifti(xifti0, data_ii)
    write_xifti(xifti_ii, cifti_fname = BOLD_fnames[ii])
  }

  #ESTIMATE AND SAVE TEMPLATES
  template <- estimate_template(BOLD = BOLD_fnames[1:500], #training subjects
                                GICA = GICA/20,
                                scale = 'none',
                                Q2 = 0,
                                brainstructures = c('left','right'),
                                FC = TRUE)
  saveRDS(template, 'templates_longsim/template.RDS')
  
  #PLOT FC TEMPLATE (Inverse-Wishart)
  template_FC_mean <- template$template$FC$psi/(template$template$FC$nu - Q - 1)
  template_FC_var <- template_FC_mean*0
  for(q1 in 1:Q){
    for(q2 in 1:Q){
      template_FC_var[q1,q2] <- templateICAr:::IW_var(template$template$FC$nu, Q, template_FC_mean[q1,q2], template_FC_mean[q1,q1], template_FC_mean[q2,q2])
    }
  }
  pdf('templates_longsim/FCtemplate.pdf', height=5, width=5.5)
  zlim_FC <- c(-0.8, 0.8)
  diag(template_FC_mean) <- diag(template_FC_var) <- NA
  plot_FC(template_FC_mean, zlim=zlim_FC, break_by = 0.2, title = "FC Template Mean", cor=FALSE)
  plot_FC(sqrt(template_FC_var), zlim=c(0.1, 0.3), break_by = 0.1, cols = viridis(12), title = "FC Template SD", cor=FALSE)
  #plot_FC(template_FC_mean/TC_sim_cor_avg, zlim=c(0.50,1.0), break_by = 0.05, digits_legend = 2, cols=c(heat.colors(12), 'white'), title = "FC Template Mean / Population Mean", cor=FALSE)
  #plot_FC(sqrt(template_FC_var)/sqrt(TC_sim_cor_var), zlim=c(1,2), break_by = 0.2, digits_legend = 2, cols=rev(c(heat.colors(12), 'white')), title = "FC Template SD / Population SD", cor=FALSE)
  dev.off()
  
  #PLOT FC TEMPLATE (Cholesky)
  template_FC_mean2 <- template$template$FC_Chol$FC_samp_mean
  template_FC_var2 <- template$template$FC_Chol$FC_samp_var
  pdf('templates_longsim/FCtemplate2.pdf', height=5, width=5.5)
  zlim_FC <- c(-0.8, 0.8)
  diag(template_FC_mean2) <- diag(template_FC_var2) <- NA
  plot_FC(template_FC_mean2, zlim=zlim_FC, break_by = 0.2, title = "Cholesky FC Template Mean", cor=FALSE)
  plot_FC(sqrt(template_FC_var2), zlim=c(0.1, 0.3), break_by = 0.1, cols = viridis(12), title = "Cholesky FC Template SD", cor=FALSE)
  dev.off()
  
  #PLOT IC TEMPLATE
  #built-in templateICAr function plots IC template means/vars
  plot(template, fname = 'templates_longsim/template_est', idx=1:5, stat='mean', zlim=c(-0.5,0.5))
  plot(template, fname = 'templates_longsim/template_est', idx=1:5, stat='var', zlim=c(0,0.025))

}



###################################################################################
# GENERATE LONG TIME COURSES for TEST SUBJECTS

#number of test subjects
n <- 50
inds_test <- 501:(500+n)

#delete this once fMRItools::UT2mat is fixed.  Then change UT2mat to fMRItools::UT2mat
UT2mat <- function(x, diag = 1){
  x <- fMRItools::UT2mat(x) #form 5x5 matrix
  x <- x + t(x)
  diag(x) <- diag
  x
}

nT_long <- 1200*4 #4 ~15 min sessions = one hour
if(first_run){
  
  #real fMRI data time courses based on dual regression
  TCs <- readRDS('TCs_real/TCs.RDS') #TxQxn (n=subjects) -- 1113 sessions
  TCs <- TCs[,,!is.na(TCs[1,1,])] #exclude missing subjects -- 1068 sessions remaining
  nS <- dim(TCs)[3]
  nQ <- dim(TCs)[2]
  
  #define true FC for each subject
  FC_true <- apply(TCs, 3, function(x) fMRItools::mat2UT(cor(x)))
  
  #function to generate long multivariate time courses with a given QxQ correlation structure
  TC_fun <- function(FC){
    #get the FC correlation matrix
    FC_mat <- UT2mat(FC) #form 5x5 matrix
    #generate long time courses with the desired FC
    Z <- matrix(rnorm(nT_long*5), 5, nT_long) #generate independent N(0,1) data (no correlation between networks yet)
    X <- expm::sqrtm(FC_mat) %*% Z #induce correlation
  }
  
  #for each real FC matrix, generate a long TC matrix (TxQ) as MVN variates
  TCs_FC <- array(apply(FC_true, 2, TC_fun), dim = c(nQ, nT_long, nS))
  save(FC_true, TCs_FC, file='TCs_FC.RData')
}
load(file='TCs_FC.RData')



# ---- Plot mean and var of FC
if(first_run){
  TC_sim_cor_avg <- apply(FC_true, 1, mean, na.rm=TRUE) #similar to TCs_cor_avg
  TC_sim_cor_var <- apply(FC_true, 1, var, na.rm=TRUE) #similar to TCs_cor_var
  TC_sim_cor_avg <- UT2mat(TC_sim_cor_avg, diag = 1)
  TC_sim_cor_var <- UT2mat(TC_sim_cor_var, diag = 0)
  
  #compute and plot FC mean & var
  pdf('plots_longsim/FC_mean_and_SD.pdf', height=5, width=5.5)
  plot_FC(TC_sim_cor_avg, zlim=c(-0.8, 0.8), break_by = 0.2, title = "True Mean of FC", cor=TRUE)
  plot_FC(sqrt(TC_sim_cor_var), zlim=c(0.1, 0.3), break_by = 0.1, cols = viridis(12), title = "True SD of FC", cor=TRUE)
  dev.off()
}


################################################################
### GENERATE TEST SUBJECT DATA & RUN MODELS (VARY DURATION)


## Setup for data generation

#for reproducibility
set.seed(734859) #set seed before generating seeds for each subject
seeds <- runif(n, min=0, max=1000000) #generate a seed for each subject's noise

#first, determine signal variance for SNR to fix error variance
var_TC <- 1 #pre-set based on scaling
tvar <- rep(0,Q)
#identify intensity of peak voxels for each IC (top 1% of each IC)
peak_vals <- apply(template_mean_est, 2, function(x) {
  mean(x[x > quantile(x, 0.99)])
})
#signal variance = variance of signal time courses scaled by peak IC intensity, averaged over ICs
sd_sig <- sqrt(mean(var_TC*peak_vals)) #0.93 approx
sd_err <- 2*sd_sig #1.86 approx

## Setup for running models

#Note: In the original simulation setting, V = 2652 which is much smaller than in typical real data analyses.
#Therefore, we only consider durations up to 1200 which is roughly half of the "observations" when estimating the TxQ "A" matrix
#If we allow T to grow to the size of V or beyond, while V remains fixed, we do not expect the estimates to converge

durations <- c(200, 400, 800, 1200, 1600, 2400) 
num_dur <- length(durations)

#what is the difference between the true FC and the FC of the generated time series?
#this is a lower bound for the accuracy of the algorithm
FC_MSE_min <- rep(NA, num_dur)
for(dd in 1:num_dur){
  dur <- durations[dd]
  FC_dd <- apply(TCs_FC[,1:dur,], 3, function(x) mat2UT(cor(t(x))))
  FC_MSE_min[dd] <- mean((FC_dd - FC_true)^2)
}

algos <- c('VB1','VB2') 
num_alg <- length(algos)

template <- readRDS('../simulation/templates_longsim/template.RDS')
FC_mean_VB1 <- template$template$FC$psi/(template$template$FC$nu - Q - 1)
FC_mean_VB2 <- template$template$FC_Chol$FC_samp_mean
FC_mean_VB1_rep <- abind(rep(list(FC_mean_VB1), num_dur), along=3)
FC_mean_VB2_rep <- abind(rep(list(FC_mean_VB2), num_dur), along=3)
FC_mean_rep <- abind(FC_mean_VB1_rep, FC_mean_VB2_rep, along=4)

MSE_true <- MSE_prior <- MSE2_true <- MSE2_prior <- CI_cover <- CI_width <- array(0, dim = c(Q, Q, num_dur, num_alg))
MSE_TC <- array(0, dim = c(num_dur, num_alg)) #MSE of TxQ "A" matrices (compute iteratively)

xifti0 <- ciftiTools:::convert_to_dtseries(newdata_xifti(GICA, matrix(0, N, nT_long))) #blank xifti object
NN <- 0
for(ii in 1:n){
  
  print(paste0('~~~~~~ Subject: ',ii,' ~~~~~~~~~'))

  # ################################################################
  # ### GENERATE DATA (maximum duration)
  # 
  # TC_ii <- t(TCs_FC[,,ii])
  # IC_ii <- subjICs[,,inds_test[ii]]
  # set.seed(seeds[ii])
  # errs_ii <- matrix(rnorm(N*nT_long, mean=0, sd=sd_err), nrow=N, ncol=nT_long) #random noise
  # data_ii <-  IC_ii %*% t(TC_ii) + errs_ii
  # xifti_ii <- newdata_xifti(xifti0, data_ii)
  # FC_true_ii <- UT2mat(FC_true[,ii])
  # 
  # ################################################################
  # ### RUN MODELS, VARY DURATION
  # 
  # FC_est_ii <- FC_est2_ii <- FC_LB_ii <- FC_UB_ii <- array(NA, dim = c(Q, Q, num_dur, num_alg))
  # MSE_TC_ii <- array(NA, dim = c(num_dur, num_alg)) #subject-level squared error of estimated vs. true time courses
  # 
  # 
  # for(dd in 1:num_dur){
  # 
  #   dur <- durations[dd]
  #   print(paste0('~~~~~~ duration: ',dur,' ~~~~~~~~~'))
  # 
  #   #loop over algorithms
  #   for(aa in 1:num_alg){
  # 
  #     print(alg <- algos[aa])
  # 
  #     result_ad <- templateICA(BOLD = xifti_ii,
  #                              template = template,
  #                              scale = 'none',
  #                              Q2 = 0,
  #                              epsilon = 1e-6,
  #                              brainstructures = c('left','right'),
  #                              verbose=FALSE,
  #                              method_FC = alg,
  #                              reduce_dim = FALSE,
  #                              time_inds = 1:dur)
  # 
  #     FC_est_ii[,,dd,aa] <- cor(result_ad$A)
  #     FC_est2_ii[,,dd,aa] <- result_ad$FC$mean
  #     FC_LB_ii[,,dd,aa] <- result_ad$FC$LB
  #     FC_UB_ii[,,dd,aa] <- result_ad$FC$UB
  # 
  #     MSE_TC_ii[dd,aa] <- mean((result_ad$A - TC_ii[1:dur,])^2)
  # 
  #   } #end loop over algorithms
  # 
  # } #end loop over durations
  # 
  # save(FC_est_ii, FC_est2_ii, FC_LB_ii, FC_UB_ii, #MSE_TC_ii, 
  #       file=paste0('results/long_duration/FC_',ii,'.RData'))
  
  #DELETE ME
  load(file=paste0('results/long_duration/FC_',ii,'_VB1.RData')) #this includes MSE_TC_ii but the VB2 results do not!
  # FC_est_ii1 <- FC_est_ii
  # FC_est2_ii1 <- FC_est2_ii
  # FC_LB_ii1 <- FC_LB_ii
  # FC_UB_ii1 <- FC_UB_ii
  # 
  # load(file=paste0('results/long_duration/FC_',ii,'.RData'))
  # FC_est_ii[,,,1] <- FC_est_ii1[,,,1]
  # FC_est2_ii[,,,1] <- FC_est2_ii1[,,,1]
  # FC_LB_ii[,,,1] <- FC_LB_ii1[,,,1]
  # FC_UB_ii[,,,1] <- FC_UB_ii1[,,,1]

  ## MSE of time courses
  
  #average of MSE over subjects (compute sequentially)
  MSE_TC <- MSE_TC + MSE_TC_ii
  
  # ## difference between estimate and truth
  # 
  # FC_true_ii_rep <- abind(rep(list(FC_true_ii), num_dur), along=3)
  # FC_true_ii_rep <- abind(rep(list(FC_true_ii_rep), num_alg), along=4)
  # #MSE over subjects (compute sequentially)
  # MSE_true <- MSE_true + (FC_est_ii - FC_true_ii_rep)^2 #FC_est = cor(A)
  # MSE2_true <- MSE2_true + (FC_est2_ii - FC_true_ii_rep)^2 #FC_est = mean of samples of G
  # 
  # ## difference between estimate and prior mean
  # 
  # #MSE over subjects (compute sequentially)
  # MSE_prior <- MSE_prior + (FC_est_ii - FC_mean_rep)^2 #FC_est = cor(A)
  # MSE2_prior <- MSE2_prior + (FC_est2_ii - FC_mean_rep)^2 #FC_est = mean of samples of G
  # 
  # ## coverage and width of credible intervals
  # 
  # CI_cover_ii <- (FC_true_ii_rep >= FC_LB_ii) & (FC_true_ii_rep <= FC_UB_ii)
  # CI_cover <- CI_cover + CI_cover_ii #coverage rate of intervals (compute sequentially)
  # CI_width_ii <- FC_UB_ii - FC_LB_ii
  # CI_width <- CI_width + CI_width_ii #width of intervals (compute sequentially)
  # 
  NN <- NN + 1
  # 
  # save(MSE_true, MSE_prior, MSE2_true, MSE2_prior, CI_cover, CI_width, NN, #MSE_TC, 
  #      file = 'results/long_duration/MSE_true_prior.RData')
  # 
} #end loop over subjects

MSE2_true <- MSE2_true/NN
MSE_true <- MSE_true/NN
MSE_prior <- MSE_prior/NN
CI_cover <- CI_cover/NN
MSE_TC <- MSE_TC/NN


save(MSE_true, MSE_prior, CI_cover, MSE_TC, 
     file = 'results/long_duration/MSE.RData')
load(file = 'results/long_duration/MSE.RData')

## Plot MSE of Time Courses

#MSE_TC <- apply(MSE_TC, 1:2, mean) #average MSE over subjects
pdf(file.path('plots_longsim','MSE_TCs.pdf'), width=6, height=4)
#par(mfrow=c(1,2))
#for(aa in 1:2){
aa <- 1
  plot(durations, MSE_TC[,aa], type='l', ylab='MSE of Mixing Matrix', main=paste0('VB',aa), xlab = 'Duration', ylim=c(0.01,0.018))
  abline(h = min(MSE_TC[,aa]), lty=2)
#}
dev.off()

## Plot MSE of FC Matrices

set.seed(1234567)
cols <- brewer_pal(pal = 'Paired')(10)

#compute mean over edges
MSE_true_avg <- apply(MSE_true, 3:4, function(x){ mean(x[upper.tri(x)]) })
MSE_prior_avg <- apply(MSE_prior, 3:4, function(x){ mean(x[upper.tri(x)]) })
MSE2_true_avg <- apply(MSE2_true, 3:4, function(x){ mean(x[upper.tri(x)]) })

#loop over VB1, VB2
pdf(file.path('plots_longsim','MSE_truth_prior.pdf'), width=8, height=6)
par(mfrow=c(1,2))
for(aa in 1:2){
  plot(durations, MSE2_true[1,2,,aa], type='l', col='gray', lwd=1.5, ylim=c(0,0.12), main=paste0('VB',aa), ylab='MSE of FC', lty=1, xlab='Duration')
  lines(durations, MSE_prior[1,2,,aa], col='gray', lty=2, lwd=1.5)
  c <- 1
  for(q1 in 1:4){
    for(q2 in (q1+1):5){
      if(q1==1 & q2==2) next()
      lines(durations, MSE2_true[q1,q2,,aa], col='gray', lty=1, lwd=1.5)
      lines(durations, MSE_prior[q1,q2,,aa], col='gray', lty=2, lwd=1.5)
      c <- c+1
    }
  }
  lines(durations, MSE2_true_avg[,aa], col='black', lty=1, lwd=3) #avg over edges
  lines(durations, MSE_prior_avg[,aa], col='black', lty=2, lwd=3) #avg over edges
  abline(h=0, col='red')
    legend('topleft',legend = c('FC Estimate vs. Truth', 'FC Estimate vs. Prior Mean'), lty=c(1,2)) 
}
dev.off()

# ## Make matrix of colors corresponding to previous plot
# 
# cols_mat <- UT2mat(1:10, diag = 0)
# cols_mat[lower.tri(cols_mat)] <- NA
# pdf(file.path('plots_longsim','colors.pdf'), width=5, height=4.5)
# fMRItools::plot_FC(cols_mat, cols = c('black',cols), zlim = c(0,10), lines = 1:5, lines_col = 'black', lines_lwd = 4)
# dev.off()

apply(CI_cover, 3:4, function(x){ mean(x[upper.tri(x)])})



