library(templateICAr)
library(ciftiTools)
ciftiTools.setOption('wb_path', '/Applications/')
library(dplyr)
library(matrixStats)
library(viridis)
library(ggplot2)
library(ggthemes)

main_dir <- '~/Dropbox/RESEARCH/FCTemplateICA'
setwd(main_dir)

source('code/functions.R')


### FOR nQ=25 ONLY

nQ <- nL <- 25
result_dir <- file.path('results','GICA25')

load(file=paste0('data/inds_',nQ,'.RData')) #inds_old, inds_new, reorder, networks, change, labels, networks_names_df, network_ICs_df_ordered

#this is to correct the ordering between template estimation and "inds_old" when we manually reassigned one network
inds_orig <- c(inds_old[1:5],inds_old[25],inds_old[6:24])
#inds_fix <- c(1:5,7:25,6)

template <- readRDS(file.path(main_dir,'templates',paste0('template_',nQ,'.rds')))
template_exp <- export_template(template)
xii <- template_exp$mean
xii$meta$cifti$names <- paste0('IC',inds_orig) # NOTE THAT THE XIFTI TEMPLATE HAS THE OLD ORDERING, SO NO NEED TO REORDER SPATIAL MAPS WHEN PLOTTING

#mask for ICC summary (USING OLD IC ORDERING)
sd_Q <- sqrt(colVars(as.matrix(template_exp$mean))) #thresholds to use in activations (x2, x3, ...)
thr_mat <- matrix(sd_Q, nrow=nrow(template_exp$mean), ncol=ncol(template_exp$mean), byrow=TRUE)
template_mean_bin <- abs(as.matrix(template_exp$mean)) > thr_mat
template_mean_bin[template_mean_bin==0] <- NA


###############################################################
## RELIABILITY & UNIQUENESS OF FC
###############################################################

FC1 <- readRDS(file=file.path(result_dir,'FC_sess1.rds'))
FC2 <- readRDS(file=file.path(result_dir,'FC_sess2.rds'))
FC1 <- FC1[reorder,reorder,,]
FC2 <- FC2[reorder,reorder,,]
mat2UT <- function(x){
  x[upper.tri(x)]
}

### ICC of FC

#compute & summarize
FC_var_tot1 <- apply(FC1, c(1,2,4), var)
FC_var_tot2 <- apply(FC2, c(1,2,4), var)
FC_var_tot <- (FC_var_tot1 + FC_var_tot2)/2
FC_var_noise <- 1/2 * apply(FC1 - FC2, c(1,2,4), var)
FC_var_sig <- FC_var_tot - FC_var_noise
FC_var_sig[FC_var_sig < 0] <- 0
FC_ICC <- FC_var_sig/FC_var_tot
FC_ICC_change <- (FC_ICC[,,1] - FC_ICC[,,2])/FC_ICC[,,2]
summary(mat2UT(FC_ICC_change))
mean(mat2UT(FC_ICC_change) > 0, na.rm=TRUE)


#visualize
pdf(file.path(main_dir,'plots','GICA25',paste0('FC_ICC.pdf')), width=5, height=5)
cols <- viridis_pal(option = 'D')(7)
plot_FC(FC_ICC[,,1], zlim=c(0, 0.5), break_by = 0.05, digits_legend = 2, title = "FC Template ICA", cor=TRUE, lines = change, labels = labels, lwd_lines=1.5, cols=cols)#, col_lines = 'white')
plot_FC(FC_ICC[,,2], zlim=c(0, 0.5), break_by = 0.05, digits_legend = 2, title = "Template ICA", cor=TRUE, lines = change, labels = labels, lwd_lines=1.5, cols=cols)
plot_FC(FC_ICC[,,1] - FC_ICC[,,2], zlim=c(-0.05, 0.0500001), break_by = 0.01, digits_legend = 2, title = "Change vs. tICA", cor=TRUE, lines = change, labels = labels, lwd_lines=1.5, cols_rev = TRUE)
#plot_FC(FC_ICC[,,3], zlim=c(0, 0.55), break_by = 0.05, digits_legend = 2, title = "Dual Regression", cor=TRUE, lines = change, lwd_lines=1.5, cols=cols)
#note: DR NOT comparable since based on group average maps.  Higher ICC, but this could be due to functional topology differences!
dev.off()

### Size of Activations / Deviations

measures <- c('act','dev_pos','dev_neg')
ylabs <- c('Activations', 'Deviations', 'Neg. Deviations')
med_change <- rep(NA, 3)
pct_better <- rep(NA, 3)
for(m in 1:3){
  meas <- measures[m]
  S_act1 <- readRDS(file=file.path(result_dir,paste0('S_',meas,'_sess1.rds')))
  S_act2 <- readRDS(file=file.path(result_dir,paste0('S_',meas,'_sess2.rds')))
  S_act1 <- S_act1[,reorder,,]
  S_act2 <- S_act2[,reorder,,]

  #size of activations
  size1 <- apply(S_act1, 2:4, sum)
  size2 <- apply(S_act2, 2:4, sum)
  size_avg <- (size1 + size2) / 2 #avg over sessions
  size_change <- (size_avg[,,1] - size_avg[,,2])/size_avg[,,2] #%change in size (FC-tICA - tICA) for each subject and IC
  med_change[m] <- median(size_change, na.rm=TRUE)
  pct_better[m] <- mean(size_change > 0, na.rm=TRUE)

  size_df <- rbind(data.frame(size = c(size_avg[,,1]), IC = rep(inds_new, times=100), network = rep(networks, times=100), subj = rep(1:100, each=nQ), method='FC-tICA'),
                   data.frame(size = c(size_avg[,,2]), IC = rep(inds_new, times=100), network = rep(networks, times=100), subj = rep(1:100, each=nQ), method='tICA'))
  size_df$IC <- factor(size_df$IC, levels = inds_new)
  size_df$method <- factor(size_df$method, levels=c('tICA','FC-tICA'))
  size_df$network <- factor(size_df$network, levels=unique(networks_names_df$network2))
  if(m==1) breaks <- seq(0,9000,3000) else breaks <- c(0,9000)
  pdf(file.path(main_dir,'plots','GICA25',paste0('size_',meas,'.pdf')), width=10, height=3)
  print(ggplot(size_df, aes(x=IC, y=size)) +
          geom_boxplot(aes(group=interaction(method, IC), fill=method), alpha=0.5) +
          scale_fill_manual(values=c('white','turquoise')) + scale_y_continuous(breaks = breaks) +
          facet_grid(. ~ network, scales='free', space='free') + ylab(paste0('Size of ',ylabs[m])) +
          theme_few() )
  dev.off()

  #size of deviations
  size_diff_df <- data.frame(size_diff = c(size_avg[,,1] - size_avg[,,2]),
                             IC = rep(inds_new, times=100),
                             network = rep(networks, times=100),
                             subj = rep(1:100, each=nQ))
  size_diff_df$IC <- factor(size_diff_df$IC, levels = inds_new)
  size_diff_df$network <- factor(size_diff_df$network, levels=unique(networks_names_df$network2))
  size_diff_df$method <- 'FC-tICA' #this is just for plot sizing
  pdf(file.path(main_dir,'plots','GICA25',paste0('sizediff_',meas,'.pdf')), width=10, height=3)
  print(ggplot(size_diff_df, aes(x=IC, y=size_diff)) +
          geom_hline(yintercept=0, linetype=2) +
          geom_boxplot(aes(group=IC, fill=method), alpha=0.5) +
          scale_fill_manual(values='purple') + #this is just for plot sizing
          facet_grid(. ~ network, scales='free', space='free') + ylab('Difference (FC-tICA - tICA)') +
          theme_few())
  dev.off()

}
med_change # what is the median percentage size change across subjects and ICs?
# [1] 0.1117819 0.2361677 0.1552817
pct_better # what percentage of subject-level maps got larger?
# [1] 0.9971978 0.9984000 0.8816000


###############################################################
## RELIABILITY & UNIQUENESS OF SPATIAL IC MAPS
###############################################################


### ICC of Spatial IC Maps

S1 <- readRDS(file=file.path(result_dir,'S_sess1.rds'))
S2 <- readRDS(file=file.path(result_dir,'S_sess2.rds'))

# #compute & summarize
# S_err <- abs(S1 - S2)
# S_MAE <- apply(S_err, c(1,2,4), median)
# #Calculate change in FC MAE vs. tICA
# S_MAE_change <- (S_MAE[,,1] - S_MAE[,,2])/S_MAE[,,2]
#
# plot(newdata_xifti(xii, S_MAE[,,1]), idx=1:nL, zlim=c(0,0.1), fname=file.path(main_dir,'images','GICA25','S_MAE_FCtICA'))
# plot(newdata_xifti(xii, S_MAE[,,2]), idx=1:nL, zlim=c(0,0.1), fname=file.path(main_dir,'images','GICA25','S_MAE_tICA'))
# plot(newdata_xifti(xii, S_MAE[,,3]), idx=1:nL, zlim=c(0,0.1), fname=file.path(main_dir,'images','GICA25','S_MAE_DR'))
# plot(newdata_xifti(xii, S_MAE_change), idx=1:nL, zlim=c(-0.1,0.1), fname=file.path(main_dir,'images','GICA25','S_MAE_change_vs_tICA'))
# #plot(newdata_xifti(template_exp$mean, rowMeans(S_MAE[,,1])), zlim=c(0,0.1), fname=file.path(main_dir,'images','GICA25','S_MAE_FCtICA_avg'))
# #plot(newdata_xifti(template_exp$mean, rowMeans(S_MAE[,,2])), zlim=c(0,0.1), fname=file.path(main_dir,'images','GICA25','S_MAE_tICA_avg'))
# #plot(newdata_xifti(template_exp$mean, rowMeans(S_MAE[,,3])), zlim=c(0,0.2), fname=file.path(main_dir,'images','GICA25','S_MAE_DR_avg'))
# plot(newdata_xifti(xii, rowMeans(S_MAE_change)), zlim=c(-0.1,0.1), fname=file.path(main_dir,'images','GICA25','S_MAE_change_vs_tICA_avg'))


#compute & summarize
var_tot1 <- apply(S1, c(1,2,4), var)
var_tot2 <- apply(S2, c(1,2,4), var)
var_tot <- (var_tot1 + var_tot2)/2
var_noise <- 1/2 * apply(S1 - S2, c(1,2,4), var)
var_sig <- var_tot - var_noise
var_sig[var_sig < 0] <- 0
ICC <- var_sig/var_tot


#mean ICC within binarized areas
ICC_avg_FCtICA <- colMeans(ICC[,,1]*template_mean_bin, na.rm=TRUE); mean(ICC_avg_FCtICA) #0.4052276
ICC_avg_tICA <- colMeans(ICC[,,2]*template_mean_bin, na.rm=TRUE); mean(ICC_avg_tICA) #0.4036792
ICC_avg_DR <- colMeans(ICC[,,3]*template_mean_bin, na.rm=TRUE); mean(ICC_avg_DR) #0.3605865
mean((ICC_avg_FCtICA - ICC_avg_tICA)/ICC_avg_FCtICA) #0.005696175
mean((ICC_avg_FCtICA - ICC_avg_DR)/ICC_avg_FCtICA) #0.1055397


plot(newdata_xifti(xii, ICC[,,1]), idx=1:nL, zlim=c(0,0.6), colors='viridis', fname=file.path(main_dir,'images','GICA25','S_ICC_FCtICA'))
plot(newdata_xifti(xii, ICC[,,2]), idx=1:nL, zlim=c(0,0.6), colors='viridis', fname=file.path(main_dir,'images','GICA25','S_ICC_tICA'))
plot(newdata_xifti(xii, ICC[,,3]), idx=1:nL, zlim=c(0,0.6), colors='viridis', fname=file.path(main_dir,'images','GICA25','S_ICC_DR'))


### Computation Time

comptime1_25 <- readRDS(file=file.path('results','GICA25','comptime_sess1.rds'))
comptime2_25 <- readRDS(file=file.path('results','GICA25','comptime_sess2.rds'))
comptime_25 <- (comptime1_25 + comptime2_25)/2 #avg over sessions
mean(comptime_25, na.rm=TRUE) #around 2192 seconds or 37 minutes

comptime_100 <- readRDS(file=file.path('results','GICA100','comptime_sess1.rds'))
comptime_100 #10836 seconds or 3 hours

### Make table of RSN assignments for paper

library(xtable)

load(file=paste0('data/IC_networks_ordered_25.RData')) #IC_order, badinds, network_ICs_df_ordered
df25 <- network_ICs_df_ordered %>% group_by(network) %>% summarize(count25 = n())

load(file=paste0('data/IC_networks_ordered_100.RData')) #IC_order, badinds, network_ICs_df_ordered
df100 <- network_ICs_df_ordered %>% group_by(network) %>% summarize(count100 = n())

df <- full_join(df25, df100)
df$count25[is.na(df$count25)] <- 0
df$count100[is.na(df$count100)] <- 0

dfx <- xtable(df, digits=0)
print(dfx, include.rownames=FALSE)

### Make full table of RSN assignments for paper

# Q = 25

df <- network_ICs_df_ordered[,2:1]
dfx <- xtable(df, digits=0)
print(dfx, include.rownames=FALSE)

# Q = 100

load(file=paste0('data/IC_networks_ordered_100.RData')) #IC_order, badinds, network_ICs_df_ordered
network_ICs_df_ordered <- network_ICs_df_ordered[network_ICs_df_ordered$IC != 74,]
df <- cbind(network_ICs_df_ordered[1:33,2:1],network_ICs_df_ordered[34:66,2:1],network_ICs_df_ordered[67:99,2:1])
dfx <- xtable(df, digits=0)
print(dfx, include.rownames=FALSE)

