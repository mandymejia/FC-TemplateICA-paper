devtools::install_github('mandymejia/templateICAr', ref='8.0')
library(templateICAr) #0.8.0
library(ciftiTools) #0.14.0
#### REINSTALL THIS ONE #library(fMRItools)
#### ANI CHANGE next two lines
#ciftiTools.setOption('wb_path', '/Applications')
ciftiTools.setOption('wb_path', '/Applications/workbench/bin_macosxub/')

#library(matrixStats) #colVars
library(reshape2) #1.4.4
library(ggplot2) #3.5.1
library(ggthemes) #5.1.0
library(scales) #alpha() #1.3.0
library(RColorBrewer) #1.1-3
library(abind) #1.4-5
library(viridis) #0.6.5
library(dplyr) #1.1.4

###ANI CHANGE next two lines
#main_dir <- '~/Dropbox/RESEARCH/FCTemplateICA/simulation/'
main_dir <- '~/Dropbox (Brown)/FCTemplateICA/simulation/'
setwd(main_dir)
source('sim_funs.R')

first_run <- FALSE

run_part1 <- TRUE
run_part2 <- FALSE

####################################################################
### GENERATE SIMULATED DATA
####################################################################

ICs <- c(1,3,4,2,16)
Q <- length(ICs)

GICA_fname <- '../data/melodic_IC_25.dscalar.nii' #avoid cloud-based files!
GICA <- read_cifti(GICA_fname, resamp_res = 3000, brainstructures='left')
GICA <- newdata_xifti(GICA, as.matrix(GICA)[,ICs])

template_mean <- as.matrix(GICA/20)
template_mean[template_mean > 1] <- 1
template_mean <- newdata_xifti(GICA, template_mean)
template_var <- (template_mean/1.5)^2
N <- nrow(template_mean)

if(first_run){
  plot(template_mean, idx=1:5, zlim=c(-0.5,0.5),
       title = paste0('Generating Mean, IC ', 1:5),
       fname = paste0('templates/generating_mean', 1:5))

  plot(template_var, idx=1:5, zlim=c(0,0.05),
       title = paste0('Generating Var, IC ', 1:5),
       fname = paste0('templates/generating_var', 1:5))
}

###################################################################################
# GENERATE SUBJECT ICS

# 1. Draw subject effects (deviations) from N(0,sd_v), where sd_v is from template_var at voxel v
# 2. Spatially smooth with FWHM = 5
# 3. Add to template_mean to get subject ICs
# 4. Re-estimate templates (variance is reduced due to smoothing)

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
    if(ii %in% 1:3){
      subjICs_xifti <- newdata_xifti(GICA, subjICs[,,ii])
      dev_xifti <- newdata_xifti(GICA, subjICs[,,ii] - tmean_mat)
      plot(subjICs_xifti, idx=1:5, zlim=c(-0.5,0.5),
           title = paste0('Subject ',ii,', IC ', 1:5),
           fname = paste0('exampleICs/subject',ii,'_IC', 1:5),
           hemisphere = 'left')
      plot(dev_xifti, idx=1:5, zlim=c(-0.1,0.1),
           title = paste0('Subject ',ii,', Deviation ', 1:5),
           fname = paste0('exampleICs/subject',ii,'_dev', 1:5),
           hemisphere = 'left')
    }
  }
  saveRDS(subjICs, 'subjICs.RDS')
}
subjICs <- readRDS('subjICs.RDS')

# David : This is a VxQ "S" matrix and surface mesh
subjICs_example <- newdata_xifti(GICA, subjICs[,,1])
surfL <- load_surf(hemisphere = 'left', name = 'midthickness')
subjICs_example <- add_surf(subjICs_example, surfL = surfL)

# ESTIMATE ORACLE TEMPLATES (mean should be similar to template_mean above, but variance will be reduced due to smoothing)
if(first_run){

  template_mean_est <- apply(subjICs[,,1:500], c(1,2), mean)
  template_var_est <- apply(subjICs[,,1:500], c(1,2), var)
  tmean_est_xifti <- newdata_xifti(GICA, template_mean_est)
  tvar_est_xifti <- newdata_xifti(GICA, template_var_est)

  plot(tmean_est_xifti, idx=1:5, zlim=c(-0.5,0.5),
       title = paste0('Template Mean, IC ', 1:5),
       fname = paste0('templates/template_mean', 1:5))

  plot(tvar_est_xifti, idx=1:5, zlim=c(0,0.03),
       title = paste0('Template Var, IC ', 1:5),
       fname = paste0('templates/template_var', 1:5))

  plot(sqrt(tvar_est_xifti), idx=1:5, zlim=c(0,0.15),
       title = paste0('Template SD, IC ', 1:5),
       fname = paste0('templates/template_SD', 1:5))
}


################################################################
### GENERATE TIME COURSES

if(first_run){
  #real fMRI data time courses based on dual regression
  TCs <- readRDS('TCs_real/TCs.RDS') #TxQxn (n=subjects) -- 1113 sessions
  TCs <- TCs[,,!is.na(TCs[1,1,])] #exclude missing subjects -- 1068 sessions remaining
  TCs <- apply(TCs, c(2,3), base::scale) #make variance = 1 for identifiability

  #plot three subject's time courses
  TCs_subj1 <- melt(TCs[300:600,,1], varnames = c('time','IC')); TCs_subj1$subject <- 1
  TCs_subj2 <- melt(TCs[300:600,,2], varnames = c('time','IC')); TCs_subj2$subject <- 2
  TCs_subj3 <- melt(TCs[300:600,,4], varnames = c('time','IC')); TCs_subj3$subject <- 3 #actual subject 3 atypical due to head motion
  TCs_example <- rbind(TCs_subj1, TCs_subj2, TCs_subj3)
  pdf('TCs_real/TCs.pdf', width=8, height=5)
  ggplot(TCs_example, aes(x=time, y=value)) + geom_line() +
    facet_grid(IC ~ subject, labeller='label_both') + theme_few()
  dev.off()

  #compute mean and variance
  TCs_cov <- abind(apply(TCs, 3, cov, simplify = FALSE), along=3) # list of FC matrices, abind into array
  TCs_cov_avg <- apply(TCs_cov, c(1,2), mean, na.rm=TRUE)
  TCs_cov_var <- apply(TCs_cov, c(1,2), var, na.rm=TRUE)

  # #estimate parameters of IW distribution based on true FC matrices (1068 subjects) & resulting IW var
  # nu_est <- templateICAr:::estimate_nu(TCs_cov_var, TCs_cov_avg) #28.163
  # psi_est <- TCs_cov_avg*(nu_est - Q - 1)
  # IW_var <- TCs_cov_avg*0
  # for(q1 in 1:Q){ for(q2 in 1:Q) IW_var[q1,q2] <- templateICAr:::IW_var(nu_est, Q, TCs_cov_avg[q1,q2], TCs_cov_avg[q1,q1], TCs_cov_avg[q2,q2]) }

  #plot true FC matrices
  pdf('TCs_real/TCs_real_FC.pdf', height=5, width=5.5)
  plot_FC(TCs_cov[,,1], zlim=c(-0.8, 0.8), break_by = 0.2, title = "True FC, Subject 1", cor=TRUE)
  plot_FC(TCs_cov[,,2], zlim=c(-0.8, 0.8), break_by = 0.2, title = "True FC, Subject 2", cor=TRUE)
  plot_FC(TCs_cov[,,4], zlim=c(-0.8, 0.8), break_by = 0.2, title = "True FC, Subject 3", cor=TRUE)
  plot_FC(TCs_cov_avg, zlim=c(-0.8, 0.8), break_by = 0.2, title = "Empirical Mean of True FC", cor=TRUE)
  plot_FC(sqrt(TCs_cov_var), zlim=c(0.15, 0.30), break_by = 0.05, cols = viridis(12), digits_legend = 2, title = "Empirical SD of True FC", cor=TRUE)
  #plot_FC(sqrt(IW_var), zlim=c(0.15, 0.30), break_by = 0.05, cols = viridis(12), digits_legend = 2, title = "Fitted Inverse Wishart SD", cor=TRUE)
  dev.off()

  #save time courses
  TC1 <- TCs[,,1:n]
  save(TC1, file='TCs_sim.RData')
}
load(file='TCs_sim.RData')
TC_sim <- TC1 # David : These are the TxQ "A" matrices for 550 subjects

# ---- Plot mean and var of FC
if(first_run){
  TC_sim_cor <- abind(apply(TC_sim, 3, cor, simplify = FALSE), along=3) # list of FC matrices, abind into array
  TC_sim_cor_avg <- apply(TC_sim_cor, c(1,2), mean, na.rm=TRUE) #similar to TCs_cor_avg
  TC_sim_cor_var <- apply(TC_sim_cor, c(1,2), var, na.rm=TRUE) #similar to TCs_cor_var
  #compute and plot FC mean & var
  pdf('plots/FC_mean_and_SD.pdf', height=5, width=5.5)
  plot_FC(TC_sim_cor_avg, zlim=c(-0.8, 0.8), break_by = 0.2, title = "True Mean of FC", cor=TRUE)
  plot_FC(sqrt(TC_sim_cor_var), zlim=c(0.1, 0.3), break_by = 0.1, cols = viridis(12), title = "True SD of FC", cor=TRUE)
  dev.off()
}

################################################################
### SETUP FOR GENERATING FMRI DATA

#generate and write out simulated fMRI data for each subject
BOLD_fnames <- paste0('data/subj',1:n,'.dtseries.nii')
ntime <- 1200 #length of fMRI timeseries
if(first_run){
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
}

  ################################################################
  ### ESTIMATE AND PLOT TEMPLATES

if(first_run){
  #ESTIMATE AND SAVE TEMPLATES
  BOLD_fnames_train <- BOLD_fnames[1:500]
    template <- estimate_template(BOLD = BOLD_fnames_train,
                                GICA = GICA/20,
                                scale = 'none',
                                Q2 = 0,
                                brainstructures = 'left',
                                FC = TRUE)
  saveRDS(template, 'templates/template.RDS')

  #PLOT FC TEMPLATE (Inverse-Wishart)
  template_FC_mean <- template$template$FC$psi/(template$template$FC$nu - Q - 1)
  template_FC_var <- template_FC_mean*0
  for(q1 in 1:Q){
    for(q2 in 1:Q){
      template_FC_var[q1,q2] <- templateICAr:::IW_var(template$template$FC$nu, Q, template_FC_mean[q1,q2], template_FC_mean[q1,q1], template_FC_mean[q2,q2])
    }
  }
  pdf('templates/FCtemplate.pdf', height=5, width=5.5)
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
  pdf('templates/FCtemplate2.pdf', height=5, width=5.5)
  zlim_FC <- c(-0.8, 0.8)
  diag(template_FC_mean2) <- diag(template_FC_var2) <- NA
  plot_FC(template_FC_mean2, zlim=zlim_FC, break_by = 0.2, title = "Cholesky FC Template Mean", cor=FALSE)
  plot_FC(sqrt(template_FC_var2), zlim=c(0.1, 0.3), break_by = 0.1, cols = viridis(12), title = "Cholesky FC Template SD", cor=FALSE)
  dev.off()

  #PLOT IC TEMPLATE
  #built-in templateICAr function plots IC template means/vars
  plot(template, fname = 'templates/template_est', idx=1:5, stat='mean', zlim=c(-0.5,0.5))
  plot(template, fname = 'templates/template_est', idx=1:5, stat='var', zlim=c(0,0.03))
  plot(newdata_xifti(GICA, template$template$varNN/template_var_est), zlim=c(1,5), fname = 'templates/template_var_bias')

  #SMALL BIAS IN TEMPLATE MEAN? probably due to dual regression
  pdf('plots/template_mean_bias.pdf', width=5, height=5)
  subjICs_mean <- apply(subjICs[,,1:500], c(1,2), mean, na.rm=TRUE)
  template_mean <- newdata_xifti(GICA, template$template$mean)
  plot(subjICs_mean[,1], as.matrix(template_mean)[,1], xlab='Mean of True ICs', ylab='Estimated Template Mean', main='IC 1')
  abline(a=0, b=1, col='red')
  plot(subjICs_mean[,2], as.matrix(template_mean)[,2], xlab='Mean of True ICs', ylab='Estimated Template Mean', main='IC 2')
  abline(a=0, b=1, col='red')
  plot(subjICs_mean[,3], as.matrix(template_mean)[,3], xlab='Mean of True ICs', ylab='Estimated Template Mean', main='IC 3')
  abline(a=0, b=1, col='red')
  plot(subjICs_mean[,4], as.matrix(template_mean)[,4], xlab='Mean of True ICs', ylab='Estimated Template Mean', main='IC 4')
  abline(a=0, b=1, col='red')
  plot(subjICs_mean[,5], as.matrix(template_mean)[,5], xlab='Mean of True ICs', ylab='Estimated Template Mean', main='IC 5')
  abline(a=0, b=1, col='red')
  dev.off()
}
template <- readRDS('templates/template.RDS')

##############################################################################
### ANALYZE TEST SUBJECTS

inds_test <- 501:550
BOLD_fnames_test <- BOLD_fnames[inds_test]
template_mean <- newdata_xifti(GICA, template$template$mean)
template_noFC <- template; template_noFC$template$FC <- NULL
template_FC_mean <- template$template$FC$psi/(template$template$FC$nu - Q - 1)
FC_mean <- cov2cor(template_FC_mean)
zlim_dev <- c(-0.3,0.3)

### SIMULATION PART 1 --------------------------------------------------------

if(run_part1){

  algos <- c('tICA','VB1','VB2') #Note: tICA run to initialize FC-tICA, but with dimension reduction. For a fair comparison, here we run tICA separately with no dimension reduction.
  num_alg <- length(algos)
  results_list <- vector('list', num_alg)
  names(results_list) <- algos

  # IC ESTIMATES
  S_true <- subjICs[,,inds_test]
  S_est <- array(NA, dim=c(N, Q, 50, num_alg))

  # IC DEVIATIONS
  Sdev_true <- S_true - array(as.matrix(template_mean), dim=dim(S_true))
  Sdev_est <- array(NA, dim=c(N, Q, 50, num_alg))

  # FC ESTIMATES
  FC_true <- abind(apply(TC_sim[,,inds_test], 3, cor, simplify = FALSE), along=3) # list of FC matrices, abind into array
  FC_est <- array(NA, dim=c(Q, Q, 50, num_alg + 1)) # +1 for dual regression
  FC_LB <- FC_UB <- array(NA, dim=c(Q, Q, 50, num_alg - 1)) # exclude tICA

  # COMPUTATION TIME, ITERATIONS
  comptime <- matrix(NA, 50, num_alg)
  numiter <- matrix(NA, 50, num_alg)

  # INTERMEDIATE CONVERGENCE LEVELS
  # eps_inter <- c(1e-3,1e-4,1e-5)
  # S_est_inter <- array(NA, dim=c(N, Q, 50, num_alg, length(eps_inter)))
  # FC_est_inter <- array(NA, dim=c(Q, Q, 50, num_alg, length(eps_inter)))
  #
  # CONVERGENCE CRITERIA (up to 100 iterations)
  #ELBO <- array(NA, dim=c(50, 100)) #for VB
  #logPost <- array(NA, dim=c(50, 100)) #for EM

  for(ii in 50:50){
    print(paste0('~~~~~~~~~~~~~~~~ SUBJECT ',ii,' ~~~~~~~~~~~~~~~~'))
    #load(file=paste0('results/testsubj',ii,'.RData'))

    results_ii <- results_list
    comptime_ii <- rep(NA, num_alg)
    for(aa in 1:num_alg){ #### COMMENT OUT THIS LINE AND RUN FOR aa=1
      print(alg <- algos[aa])
      if(alg=='tICA') {
        template_aa <- template_noFC
        method_FC_aa <- "none"
        } else {
          template_aa <- template
          method_FC_aa <- alg
        }
      time_aa <- system.time(
        result_aa <- templateICA(BOLD = BOLD_fnames_test[ii],
                                 template = template_aa,
                                 scale = 'none',
                                 Q2 = 0,
                                 epsilon = 0.001, #1e-6 -- found that this was worse based on eps_inter.
                                 #eps_inter = eps_inter,
                                 brainstructures = 'left', # For tICA at least, it looks like the areas of engagement stay very similar with additional convergence (0.001 to 1e-6), but the background areas get a little less accurate, maybe due to over-fitting noise?
                                 verbose=TRUE,
                                 method_FC = method_FC_aa,
                                 reduce_dim = FALSE))

      #save subject i results
      results_ii[[aa]] <- result_aa
      comptime_ii[aa] <- time_aa[3]

      #collect results across subjects
      comptime[ii,aa] <- comptime_ii[aa] #<- time_aa[1]
      numiter[ii,aa] <- result_aa$numiter
      S_est[,,ii,aa] <- as.matrix(result_aa$subjICmean)
      Sdev_est[,,ii,aa] <- as.matrix(result_aa$subjICmean) - template$template$mean


      if(alg == 'tICA') {
        FC_est[,,ii,num_alg+1] <- cor(result_aa$result_DR$A) #DR
        #### ADD NEW LINE FC_est[,,ii,num_alg+1] <- cor(result_aa$result_DR$A2) #DR
        FC_est[,,ii,aa] <- cor(result_aa$A) #tICA
      } else {
        FC_est[,,ii,aa] <- result_aa$FC$mean
        FC_LB[,,ii,aa-1] <- result_aa$FC$LB
        FC_UB[,,ii,aa-1] <- result_aa$FC$UB
      }
    }
    save(results_ii, comptime_ii, file=paste0('results/testsubj',ii,'.RData'))
    ##### save(results_ii, comptime_ii, file=paste0('results_dr2/testsubj',ii,'.RData'))
    
    
    ### Visualize Results 

    #IC Estimates
    S_true_xifti <- newdata_xifti(results_ii[[aa]]$subjICmean, S_true[,,ii])
    plot(S_true_xifti, idx=1:Q, zlim=c(-0.5,0.5), title='True IC', fname=paste0('images/testsubj',ii,'_S_true'))
    plot(S_true_xifti - template_mean, idx=1:Q, zlim=c(-0.25,0.25), title='True Deviation', fname=paste0('images/testsubj',ii,'_Sdev_true'))
    for(aa in 1:num_alg){
      alg <- gsub('_','',algos[aa])
      IC_est_aa <- results_ii[[aa]]$subjICmean
      #IC estimates
      plot(IC_est_aa, idx=1:Q, zlim=c(-0.5,0.5),
           title=paste0('IC Estimate (', alg, ')'),
           fname=paste0('images/testsubj',ii,'_S_',alg))
      #Deviation estimates
      plot(IC_est_aa - template_mean, idx=1:Q, zlim=c(-0.25,0.25),
           title=paste0('Deviation Estimate (', alg, ')'),
           fname=paste0('images/testsubj',ii,'_Sdev_',alg))
      if(alg=='tICA'){
        IC_est_DR <- newdata_xifti(results_ii[[aa]]$subjICmean, t(results_ii[[aa]]$result_DR$S))
        plot(IC_est_DR - template_mean, idx=1:Q, zlim=c(-0.25,0.25),
                   title=paste0('Deviation Estimate (DR)'),
                   fname=paste0('images/testsubj',ii,'_Sdev_DR'))
      }
    }

    #FC Estimates
    zlim_cor <- c(-0.6,0.6)
    pdf(paste0('images/testsubj',ii,'FC.pdf'), width=5,height=5)
    plot_FC(FC_true[,,ii], zlim=zlim_cor, break_by = 0.3, cor=TRUE, title='True FC') 
    for(aa in c(4,1:num_alg)) {
      alg <- c(algos, 'Dual Regression')[aa]
      plot_FC(FC_est[,,ii,aa], zlim=zlim_cor, break_by = 0.3, cor=TRUE, title=paste0('FC Estimate (',alg,')'))
    }
    dev.off()

    #FC Deviations
    pdf(paste0('images/testsubj',ii,'FCdev.pdf'), width=5,height=5)
    plot_FC(FC_true[,,ii] - FC_mean, zlim=zlim_cor/2, break_by = 0.1, cor=TRUE,  title='True FC Deviation') 
    for(aa in c(4,1:num_alg)) {
      alg <- c(algos, 'Dual Regression')[aa]
      plot_FC(FC_est[,,ii,aa] - FC_mean, zlim=zlim_cor/2, title=paste0('FC Deviation Estimate (',alg,')')) #break_by = 0.1, cor=TRUE, 
    }
    dev.off()

  } #end loop over test subjects

  #save results from all test subjects
  save(S_true, S_est, template_mean, Sdev_true, Sdev_est, FC_true, FC_est, FC_mean, FC_LB, FC_UB, comptime, numiter, file = 'results/test_set_results.RData')
  #### save(S_true, S_est, template_mean, Sdev_true, Sdev_est, FC_true, FC_est, FC_mean, FC_LB, FC_UB, comptime, numiter, file = 'results_dr2/test_set_results.RData')
  #load(file = 'results/test_set_results.RData')

  # #save results for part 2
  # FC_true_1200 <- FC_true
  # FC_est_1200 <- FC_est
  # save(FC_true_1200, FC_est_1200, file = 'results/by_duration/FC_1200.RData')
  #
  # #MAE of IC Spatial Maps -- using median rather than mean due to outlier subjects (high motion?)
  # MAE_tICA <- sqrt(apply((S_est[,,,1] - S_true)^2, c(1,2), median, na.rm=TRUE))
  # for(aa in 1:num_alg){
  #   alg <- gsub('_','',algos[aa])
  #   MAE_aa <- sqrt(apply((S_est[,,,aa] - S_true)^2, c(1,2), median, na.rm=TRUE)) #equals MAE
  #   plot(newdata_xifti(GICA, MAE_aa), zlim=c(0, 0.05), idx=1:5,
  #        title=paste0('rMSE (', alg, ')'), fname=paste0('images/MSE_',alg))
  #   if(aa>1)
  #     plot(newdata_xifti(GICA, (MAE_aa - MAE_tICA)/MAE_tICA), zlim=c(-0.3, 0.3), idx=1:5,
  #          title=paste0('%Difference in rMSE vs. tICA (', alg, ')'), fname=paste0('images/MSE_diff_',alg))
  # }
  #
  # #MAE of FC (including DR) -- using median rather than mean due to outlier subjects
  # pdf('plots/MAE_FC.pdf', height=5, width=5.5)
  # for(aa in c(4,1:num_alg)){
  #   alg <- gsub('_','',c(algos, 'Dual Regression')[aa])
  #   MAE_aa <- sqrt(apply((FC_est[,,,aa] - FC_true)^2, c(1,2), median, na.rm=TRUE)) #equals MAE
  #   plot_FC(MAE_aa, zlim=c(0, 0.050001), break_by = 0.01, cor=TRUE, digits_legend = 2, cols=viridis_pal(option = 'plasma')(10),
  #           title = paste0('MAE (', alg, ')'), box=TRUE, box_col='red')
  # }
  # dev.off()
  #
  # #Change in MAE of FC
  # MAE_tICA <- sqrt(apply((FC_est[,,,1] - FC_true)^2, c(1,2), median, na.rm=TRUE))
  # MAE_DR <- sqrt(apply((FC_est[,,,4] - FC_true)^2, c(1,2), median, na.rm=TRUE))
  # pdf('plots/MAE_FC_diff.pdf', height=5, width=5.5)
  # for(aa in c(1:num_alg)){
  #   alg <- gsub('_','',algos[aa])
  #   MAE_aa <- sqrt(apply((FC_est[,,,aa] - FC_true)^2, c(1,2), median, na.rm=TRUE))
  #   if(aa > 1) plot_FC(MAE_aa - MAE_tICA, zlim=c(-0.03, 0.030001), break_by = 0.01, cor=TRUE, digits_legend = 2,
  #           title = paste0('MAE Diff. (', alg, ' vs. tICA)'), box=TRUE, box_col='red')
  #   plot_FC(MAE_aa - MAE_DR, zlim=c(-0.03, 0.030001), break_by = 0.01, cor=TRUE, digits_legend = 2,
  #           title = paste0('MAE Diff. (', alg, ' vs. DR)'), box=TRUE, box_col='red')
  # }
  # dev.off()
  #
  # #Computation Time (in seconds)
  # summary(comptime)
  # comptime_df <- as.data.frame(comptime); names(comptime_df) <- algos; comptime_df$subject <- factor(1:50)
  # comptime_df_long <- melt(comptime_df, variable.name = 'algorithm', value.name = 'comptime', id.vars = 'subject')
  # comptime_df_long$algorithm <- factor(comptime_df_long$algorithm, levels=c('tICA','VB','EM'))
  # pdf('plots/comptime.pdf', height=5, width=5)
  # ggplot(comptime_df_long, aes(x=algorithm, y=comptime/60)) +
  #   geom_line(aes(group=subject), alpha=0.1, color='purple') +
  #   geom_boxplot(aes(group=algorithm), alpha=0.5) + theme_few() +
  #   theme(legend.position='none') +
  #   ylab('Computation Time (min)') + xlab('')
  # dev.off()
  #
  # #Number of Iterations
  # numiter_df <- as.data.frame(numiter); names(numiter_df) <- algos; numiter_df$subject <- factor(1:50)
  # numiter_df_long <- melt(numiter_df, variable.name = 'algorithm', value.name = 'numiter', id.vars = 'subject')
  # numiter_df_long$algorithm <- factor(numiter_df_long$algorithm, levels=c('tICA','VB','EM'))
  # pdf('plots/numiter.pdf', height=5, width=5)
  # ggplot(numiter_df_long, aes(x=algorithm, y=numiter)) +
  #   geom_line(aes(group=subject), alpha=0.1, color='purple') +
  #   geom_boxplot(aes(group=algorithm), alpha=0.5) + theme_few() +
  #   theme(legend.position='none') +
  #   ylab('Number of Iterations (Tolerance = 0.001)') + xlab('')
  # dev.off()
  #

} #end run part 1



#   ### PART 2: VARY SCAN DURATION ---------------------------------------------
#
#   if(run_part2){
#
#     #vary duration of timeseries
#     durations <- c(100, 200, 400, 800)
#     num_dur <- length(durations)
#
#     # FC ESTIMATES (truth differs by duration)
#     FC_true <- array(NA, dim=c(Q, Q, 50, num_dur))
#     for(dd in 1:num_dur){
#       dur <- durations[dd]
#       FC_true[,,,dd] <- abind(apply(TC_sim[1:dur,,inds_test], 3, cor, simplify = FALSE), along=3) # list of FC matrices, abind into array
#     }
#
#     algos <- c('tICA', 'VB')
#     num_alg <- length(algos)
#     FC_est <- array(NA, dim=c(Q, Q, 50, num_dur, num_alg+1)) # + 1 for DR
#
#     results_list_dur <- vector('list', num_dur)
#     names(results_list_dur) <- durations
#     results_list_alg <- vector('list', num_alg)
#     names(results_list_alg) <- algos
#
#     for(ii in 1:50){
#       print(paste0('~~~~~~~~~~~~~~~~ SUBJECT ',ii,' ~~~~~~~~~~~~~~~~'))
#
#       #load(file=paste0('results/by_duration/testsubj',ii,'_bydur.RData'))
#       results_ii <- results_list_dur
#       comptime_ii <- matrix(NA, num_alg, num_dur)
#
#       for(dd in 1:num_dur){
#         dur <- durations[dd]
#         print(paste0('~~~~~~ duration: ',dur,' ~~~~~~~~~'))
#
#         results_ii[[dd]] <- results_list_alg
#
#         for(aa in 1:num_alg){
#           print(alg <- algos[aa])
#           if(alg=='tICA') {
#             template_aa <- template_noFC
#             method_FC_aa <- NULL
#           } else {
#             template_aa <- template
#             method_FC_aa <- alg
#           }
#
#             time_ad <- system.time(
#               result_ad <- templateICA(BOLD = BOLD_fnames_test[ii],
#                                      template = template_aa,
#                                      scale = 'none',
#                                      Q2 = 0,
#                                      epsilon = 0.001,
#                                      brainstructures = 'left',
#                                      verbose=FALSE,
#                                      method_FC = method_FC_aa,
#                                      reduce_dim = FALSE,
#                                      time_inds = 1:dur))
#           print(time_ad)
#           results_ii[[dd]][[aa]] <- result_ad
#           comptime_ii[aa,dd] <- time_ad[3]
#
#           FC_est[,,ii,dd,aa] <- cor(result_ad$A) #for tICA and FC-tICA (VB)
#           if(alg == 'tICA') FC_est[,,ii,dd,num_alg+1] <- cor(result_ad$result_DR$A) #DR
#
#         } #end loop over algorithms
#       } #end loop over durations
#
#       #save(results_ii, comptime_ii, file=paste0('results/by_duration/testsubj',ii,'_bydur.RData'))
#
#     } #end loop over subjects
#
#     #save(FC_true, FC_est, file='results/by_duration/FC_all_bydur.RData')
#     #load(file='results/by_duration/FC_all_bydur.RData')
#
#     #incorporate results with T=1200
#     load(file = 'results/by_duration/FC_1200.RData') #FC_true_1200, FC_est_1200
#     FC_true <- abind(FC_true, FC_true_1200, along=4)
#     FC_est_1200 <- FC_est_1200[,,,c(1,3,4)] #exclude EM (contains tICA, EM, VB, DR)
#     FC_est <- abind(FC_est, FC_est_1200, along=4)
#     durations <- c(durations, 1200) #100  200  400  800 1200
#     num_dur <- num_dur + 1
#
#     #MAE of FC (including DR) -- using median rather than mean due to outlier subject 4
#     pdf('plots/MAE_FC_bydur.pdf', height=5, width=5.5)
#     for(aa in 1:(num_alg+1)){
#       alg <- gsub('_','',c(algos, 'Dual Regression')[aa])
#       for(dd in 1:num_dur){
#         MAE_ad <- sqrt(apply((FC_est[,,,dd,aa] - FC_true[,,,dd])^2, c(1,2), median, na.rm=TRUE))
#         plot_FC(MAE_ad, zlim=c(0, 0.050001), break_by = 0.01, cor=TRUE, digits_legend = 2, cols=viridis_pal(option = 'plasma')(10),
#                 title = paste0('MAE (', alg, '), T=',durations[dd]), box=TRUE, box_col='red')
#       }
#     }
#     dev.off()
#
#     IC_names <- c('V1','V2','V3','DMN','M')
#     rows <- matrix(IC_names, nrow=Q, ncol=Q)
#     cols <- matrix(IC_names, nrow=Q, ncol=Q, byrow=TRUE)
#     FC_pairs <- matrix(paste0(rows, '-', cols), nrow=Q, ncol=Q)
#     FC_pairs_UT <- FC_pairs[upper.tri(FC_pairs)]
#
#     FC_est_list <- list(DR = FC_est[,,,,3],
#                    tICA = FC_est[,,,,1],
#                    VB = FC_est[,,,,2])
#     FC_est_err <- lapply(FC_est_list, function(x) abs(x - FC_true)) #absolute error of FC
#     FC_est_err_UT <- lapply(FC_est_err, function(myarray){ #extract upper triangles
#       apply(myarray, 3:4, function(x) x[upper.tri(x)])
#     })
#     FC_est_err_UT_long <- lapply(FC_est_err_UT, function(x){ #reorganize each list element into a data frame
#       df_FC_err <- NULL
#       for(ii in 1:50){
#         for(dd in 1:num_dur){
#           df_id <- data.frame(subj = ii, dur = durations[dd], FC_pair = FC_pairs_UT, FC_err = x[,ii,dd])
#           df_FC_err <- rbind(df_FC_err, df_id)
#         }}
#       return(df_FC_err)}
#       )
#     FC_est_err_UT_long[[1]]$method <- 'DR'
#     FC_est_err_UT_long[[2]]$method <- 'tICA'
#     FC_est_err_UT_long[[3]]$method <- 'FC-tICA (VB)'
#     FC_err_df <- do.call(rbind, FC_est_err_UT_long) #rbind list elements together to make one big data frame
#     FC_err_df$sim <- sim
#     saveRDS(FC_err_df, file='results/by_duration/FC_err_bydur.RData')
#
#     FC_err_df_summ <- FC_err_df %>% group_by(dur, FC_pair, method) %>% summarize(mean = mean(FC_err, na.rm=TRUE), SD = sd(FC_err, na.rm=TRUE))
#     #FC_err_df_summ_aug <- data.frame(dur=1, FC_pair=FC_pairs[!upper.tri(FC_pairs)], method='DR', mean=1, SD=1)
#     #FC_err_df_summ <- rbind(FC_err_df_summ, FC_err_df_summ_aug)
#     #FC_err_df_summ$FC_pair <- factor(FC_err_df_summ$FC_pair, levels = c(FC_pairs))
#     FC_err_df_summ$IC1 <- gsub('-.+','',FC_err_df_summ$FC_pair)
#     FC_err_df_summ$IC2 <- gsub('.+-','',FC_err_df_summ$FC_pair)
#     FC_err_df_summ$IC1 <- factor(FC_err_df_summ$IC1, levels=IC_names)
#     FC_err_df_summ$IC2 <- factor(FC_err_df_summ$IC2, levels=IC_names)
#     FC_err_df_summ$method <- factor(FC_err_df_summ$method, levels=c('DR','tICA','FC-tICA (VB)'))
#
#     pdf('plots/MAE_FC_bydur_lines.pdf')
#     ggplot(FC_err_df_summ, aes(x=dur, y=mean, color=method)) +
#       geom_line(aes(group=method)) +
#       geom_hline(yintercept=0, linetype=2) +
#       facet_grid(IC1 ~ IC2) +
#       scale_color_manual(values = c('darkgrey','orange','royalblue')) +
#       theme_few() + theme(legend.position='bottom') +
#       scale_x_continuous(breaks=c(400,800,1200)) + xlab('Time Series Duration') + ylab('Mean Absolute Error of FC')
#     dev.off()
#
#   }#end simulation part 2
# }#end loop over simulations
#
# # Plots of MAE vs. duration for both scenarios
#
# setwd(main_dir)
# FC_err_df_sim1 <- readRDS('sim1/results/by_duration/FC_err_bydur.RData')
# FC_err_df_sim2 <- readRDS('sim2/results/by_duration/FC_err_bydur.RData')
# FC_err_df <- rbind(FC_err_df_sim1, FC_err_df_sim2)
# FC_err_df_summ <- FC_err_df %>% group_by(sim, dur, FC_pair, method) %>% summarize(mean = median(FC_err, na.rm=TRUE), SD = sd(FC_err, na.rm=TRUE), SE = SD/sqrt(sum(!is.na(FC_err))))
# FC_err_df_summ$sim <- factor(FC_err_df_summ$sim, levels=c(1,2), labels=c('Scenario A', 'Scenario B'))
# FC_err_df_summ$method <- factor(FC_err_df_summ$method, levels=c('DR','tICA','FC-tICA (VB)'))
# FC_err_df_summ_V <- filter(FC_err_df_summ, FC_pair %in% c('V1-V2','V1-V3','V2-V3'))
#
# pdf('MAE_FC_bydur_bothsim.pdf', width=7, height=5)
# ggplot(FC_err_df_summ_V, aes(x=dur, y=mean, color=method)) +
#   geom_line(aes(group=method)) +
#   geom_ribbon(aes(group=method, ymin=mean - 2*SE, ymax = mean + 2*SE, fill=method), alpha=0.5, color=NA) +
#   #geom_hline(yintercept=0) +
#   geom_hline(data=filter(FC_err_df_summ_V, dur==200, method=='FC-tICA (VB)'), aes(yintercept = mean), color='royalblue', linetype=2) +
#   facet_grid(sim ~ FC_pair) +
#   scale_color_manual(values = c('darkgrey','orange','royalblue')) +
#   scale_fill_manual(values = c('darkgrey','orange','royalblue')) +
#   theme_few() + theme(legend.position='bottom') +
#   scale_x_continuous(breaks=durations[-1]) + xlab('Time Series Duration') + ylab('MAE of FC')
# dev.off()
#
# pdf('MAE_FC_bydur_bothsim_forR01.pdf', width=6, height=4)
# FC_err_df_summ_V <- filter(FC_err_df_summ_V, dur >= 400, sim == 'Scenario A')
# ggplot(FC_err_df_summ_V, aes(x=dur, y=mean, color=method)) +
#   geom_line(aes(group=method)) +
#   geom_ribbon(aes(group=method, ymin=mean - 2*SE, ymax = mean + 2*SE, fill=method), alpha=0.5, color=NA) +
#   #geom_hline(yintercept=0) +
#   geom_hline(data=filter(FC_err_df_summ_V, dur==400, method=='FC-tICA (VB)'), aes(yintercept = mean), color='royalblue', linetype=2) +
#   facet_wrap( ~ FC_pair, nrow=1, scales='free') +
#   scale_color_manual(values = c('darkgrey','orange','royalblue')) +
#   scale_fill_manual(values = c('darkgrey','orange','royalblue')) +
#   theme_few() + theme(legend.position='bottom') +
#   scale_x_continuous(breaks=c(400,800,1200)) + xlab('Time Series Duration') + ylab('MAE of FC')
# dev.off()
#
#
