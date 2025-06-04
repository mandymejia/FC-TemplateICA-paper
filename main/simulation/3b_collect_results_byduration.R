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

inds_test <- 501:550
n_test <- length(inds_test)
template_mean <- newdata_xifti(GICA, template$template$mean)
template_FC_mean <- template$template$FC$psi/(template$template$FC$nu - Q - 1)
FC_mean <- cov2cor(template_FC_mean)

#read in results based on T=1200 (created in 3_collect_results.R)
load(file = 'results/by_duration/FC_600.RData') #FC_est_600

#vary duration of timeseries
durations <- c(100, 200, 400, 600)
num_dur <- length(durations)

#algorithms
algos <- c('tICA','VB1','VB2') 
num_alg <- length(algos)

# True FC (held out)
FC_true <- abind(apply(TC_sim[601:1200,,inds_test], #truth based on hold-out second half of scan 
                       3, cor, simplify = FALSE), along=3) # list of FC matrices, abind into array

FC_est <- array(NA, dim=c(Q, Q, n_test, num_alg+1, num_dur)) 

for(ii in 1:n_test){
  
  print(paste0('~~~~~~~~~~~~~~~~ SUBJECT ',ii,' ~~~~~~~~~~~~~~~~'))
  
  load(file=paste0('results/by_duration/testsubj',ii,'_bydur.RData')) #results_ii
  
  for(dd in 1:(num_dur-1)){
  
    results_id <- results_ii[[dd]]
    
    #loop over algorithms
    for(aa in 2:3){ 
      
      result_aa <- results_id[[aa]]
      
      FC_est[,,ii,aa,dd] <- cor(result_aa$A) #VB
      if(aa==2){
        FC_est[,,ii,1,dd] <- cor(result_aa$result_tICA$A) #tICA
        FC_est[,,ii,4,dd] <- cor(result_aa$result_DR$A) #DR
      }
    }
  }
} #end loop over subjects

FC_est[,,,,4] <- FC_est_600[,,,1:4]

#save results from all test subjects
save(FC_est, FC_true, durations,
     file = file.path('results', 'by_duration','test_set_results.RData'))

#################################################################
### CODE FROM BY-DURATION ANALYSIS TO COLLECT RESULTS AND MAKE PLOTS

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
 
