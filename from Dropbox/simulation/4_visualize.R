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

# results from all test subjects
#load(file = 'results/test_set_results.RData') #S_true, S_est, template_mean, Sdev_true, Sdev_est, FC_true, FC_true0, FC_est, FC_mean, FC_LB, FC_UB, comptime, numiter

# # IC DEVIATIONS
# Sdev_true <- S_true - array(as.matrix(template_mean), dim=dim(S_true))
# Sdev_est <- array(NA, dim=c(N, Q, 50, num_alg + 1)) # + 1 for DR

algos <- c('tICA','VB1','VB2') 
algos2 <- c(algos, 'DR') #for spatial maps
algos3 <- c(algos2, 'DR2') #for FC matrices
num_alg <- length(algos)

load(file='results/testsubj1.RData') #results_ii
xii <- results_ii[[2]]$subjICmean

colors <- c('darkblue','turquoise','white','pink','red')
palfun <- grDevices::colorRampPalette(colors)
pal_FC <- palfun(100-1)

##############################################################################
### VISUALIZE GROUND TRUTH ICs (first 5 subjects)

load(file = file.path('results', 'test_set_results.RData')) #S_true, S_est, template_mean, Sdev_true, Sdev_est, FC_true, FC_est, FC_mean, FC_LB, FC_UB, comptime, numiter

#first 5 subjects only
for(ii in 1:5){
  
  print(paste0('~~~~~~~~~~~~~~~~ SUBJECT ',ii,' ~~~~~~~~~~~~~~~~'))
  
  #True IC maps
  S_true_xifti <- newdata_xifti(xii, S_true[,,ii])
  plot(S_true_xifti, idx=1:Q, zlim=c(-0.5,0.5), title='True IC', fname=paste0('images/examples/testsubj',ii,'_S_true'))
  plot(S_true_xifti - template_mean, idx=1:Q, zlim=c(-0.25,0.25), title='True Deviation', fname=paste0('images/testsubj',ii,'_Sdev_true'))
}

CI_width_VB1 <- CI_width_VB2 <- CI_coverage_VB1 <- CI_coverage_VB2 <- NULL

result_dir <- 'results'

load(file = file.path(result_dir, 'test_set_results.RData')) #S_true, S_est, template_mean, Sdev_true, Sdev_est, FC_true, FC_true0, FC_est, FC_mean, FC_LB, FC_UB, comptime, numiter
    
##############################################################################
### VISUALIZE IC AND FC ESTIMATES (first 5 subjects)

#first 5 subjects only
for(ii in 1:5){

  print(paste0('~~~~~~~~~~~~~~~~ SUBJECT ',ii,' ~~~~~~~~~~~~~~~~'))

  #IC Estimates
  for(aa in 1:(num_alg+1)){
    alg <- gsub('_','',algos2[aa])
    IC_est_aa <- newdata_xifti(xii, S_est[,,ii,aa])
    #IC estimates
    plot(IC_est_aa, idx=1:Q, zlim=c(-0.5,0.5),
         title=paste0('IC Estimate (', alg, ')'),
         fname=file.path('images','examples',paste0('testsubj',ii,'_S_',alg)))
    #Deviation estimates
    plot(IC_est_aa - template_mean, idx=1:Q, zlim=c(-0.25,0.25),
         title=paste0('Deviation Estimate (', alg, ')'),
         fname=file.path('images','examples',paste0('testsubj',ii,'_Sdev_',alg)))
  }

  #FC Estimates
  zlim_cor <- c(-0.6,0.6)
  pdf(file.path('plots',paste0('testsubj',ii,'FC.pdf')), width=5,height=5)
  fMRItools::plot_FC(FC_true[,,ii], zlim=zlim_cor, cols = pal_FC, lines=1:5, title='True FC', diag_val = NA)
  for(aa in c(1:(num_alg+2))) {
    alg <- algos3[aa]
    fMRItools::plot_FC(FC_est[,,ii,aa], zlim=zlim_cor, cols = pal_FC, lines=1:5, title=paste0('FC Estimate (',alg,')'), diag_val = NA)
  }
  dev.off()

  #FC Deviations
  pdf(file.path('plots',paste0('testsubj',ii,'FCdev.pdf')), width=5,height=5)
  fMRItools::plot_FC(FC_true[,,ii] - FC_mean, zlim=zlim_cor/2, cols = pal_FC, title='True FC Deviation', lines=1:5)
  for(aa in c(1:(num_alg+2))) {
    alg <- algos3[aa]
    fMRItools::plot_FC(FC_est[,,ii,aa] - FC_mean, zlim=zlim_cor/2, cleg_ticks_by = 0.1, cols = pal_FC, title=paste0('FC Deviation (',alg,')'), lines=1:5)
  }
  dev.off()

} #end loop over test subjects

##############################################################################
### VISUALIZE ACCURACY MEASURES 

#MAE of IC Spatial Maps -- using median rather than mean due to outlier subjects (high motion?)
MAE_tICA <- sqrt(apply((S_est[,,,1] - S_true)^2, c(1,2), median, na.rm=TRUE))
MAE_DR <- sqrt(apply((S_est[,,,4] - S_true)^2, c(1,2), median, na.rm=TRUE))
for(aa in 1:(num_alg+1)){
  alg <- gsub('_','',algos2[aa])
  MAE_aa <- (apply(abs(S_est[,,,aa] - S_true), c(1,2), median, na.rm=TRUE)) 
  plot(newdata_xifti(xii, MAE_aa), zlim=c(0, 0.05), idx=1:5,
       title=paste0('MAE (', alg, ')'), fname=file.path('images','MAE',paste0('MAE_',alg)))
  if(aa %in% 2:3)
    plot(newdata_xifti(xii, (MAE_aa - MAE_tICA)/MAE_tICA), zlim=c(-0.5, 0.5), idx=1:5,
         title=paste0('%Diff in MAE vs. tICA (', alg, ')'), fname=file.path('images','MAE',paste0('MAE_diff_',alg)))
  if(aa %in% 1:3)
    plot(newdata_xifti(xii, (MAE_aa - MAE_DR)/MAE_DR), zlim=c(-0.5, 0.5), idx=1:5,
         title=paste0('%Diff in MAE vs. DR (', alg, ')'), fname=file.path('images','MAE',paste0('MAE_diff2_',alg)))
}

#MAE of FC (including DR and DR2) -- using median rather than mean due to outlier subjects
pdf(file.path('plots','MAE_FC.pdf'), height=5, width=5.5)
for(aa in 1:(num_alg+2)){
  alg <- algos3[aa]
  MAE_aa <- (apply(abs(FC_est[,,,aa] - FC_true), c(1,2), median, na.rm=TRUE)) 
  plot_FC(MAE_aa, zlim=c(0, 0.18001), break_by = 0.06, cor=TRUE, digits_legend = 2, cols=viridis_pal(option = 'magma')(10),
          title = paste0('MAE (', alg, ')'), box=TRUE, box_col='springgreen')
}
dev.off()

#Change in MAE of FC
pdf(file.path('plots','MAE_FC_diff_DR.pdf'), height=5, width=5.5)
MAE_DR <- apply(abs(FC_est[,,,4] - FC_true), c(1,2), median, na.rm=TRUE)
for(aa in c(1:num_alg)){
  alg <- gsub('_','',algos[aa])
  MAE_aa <- apply(abs(FC_est[,,,aa] - FC_true), c(1,2), median, na.rm=TRUE)
  plot_FC(100*(MAE_aa - MAE_DR)/MAE_DR, zlim=c(-50, 50), break_by = 50, cor=TRUE, digits_legend = 2,
          title = paste0('MAE %Diff. (', alg, ' vs. DR)'), box=TRUE, box_col='springgreen')
}
dev.off()

# pdf(file.path('plots','MAE_FC_diff_DR2.pdf'), height=5, width=5.5)
# MAE_DR2 <- apply(abs(FC_est[,,,5] - FC_true), c(1,2), median, na.rm=TRUE)
# for(aa in c(1:num_alg)){
#   alg <- gsub('_','',algos[aa])
#   MAE_aa <- apply(abs(FC_est[,,,aa] - FC_true), c(1,2), median, na.rm=TRUE)
#   plot_FC(100*(MAE_aa - MAE_DR2)/MAE_DR2, zlim=c(-50, 50), break_by = 50, cor=TRUE, digits_legend = 2,
#           title = paste0('MAE %Diff. (', alg, ' vs. DR+)'), box=TRUE, box_col='springgreen')
# }
# dev.off()

pdf(file.path('plots','MAE_FC_diff_tICA.pdf'), height=5, width=5.5)
MAE_tICA <- apply(abs(FC_est[,,,1] - FC_true), c(1,2), median, na.rm=TRUE)
for(aa in 2:3){
  alg <- gsub('_','',algos[aa])
  if(aa==1) { plot(1:3); next() }
  MAE_aa <- apply(abs(FC_est[,,,aa] - FC_true), c(1,2), median, na.rm=TRUE)
  plot_FC(100*(MAE_aa - MAE_tICA)/MAE_tICA, zlim=c(-50, 50), break_by = 50, cor=TRUE, digits_legend = 2,
         title = paste0('MAE %Diff. (', alg, ' vs. tICA)'), box=TRUE, box_col='springgreen')
}
dev.off()

### MAE line plots

#FC pair names
IC_names <- c('V1','V2','V3','DMN','M')
rows <- matrix(IC_names, nrow=Q, ncol=Q)
cols <- matrix(IC_names, nrow=Q, ncol=Q, byrow=TRUE)
FC_pairs <- matrix(paste0(rows, '-', cols), nrow=Q, ncol=Q)
FC_pairs_UT <- FC_pairs[upper.tri(FC_pairs)]

MAE_df <- NULL
FC_true_avg <- apply(FC_true, 1:2, mean, na.rm=TRUE)
for(aa in 1:(num_alg+2)){
  alg <- algos3[aa]
  MAE_aa <- (apply(abs(FC_est[,,,aa] - FC_true), c(1,2), median, na.rm=TRUE)) 
  MAE_aa_df <- data.frame(algo = alg, 
                          FC_val = FC_true_avg[upper.tri(FC_true_avg)], 
                          MAE = MAE_aa[upper.tri(MAE_aa)], 
                          edge = FC_pairs_UT)
  MAE_df <- rbind(MAE_df, MAE_aa_df)
}

MAE_df$algo <- factor(MAE_df$algo, levels = c("DR","DR2","tICA","VB1","VB2"))
MAE_df$edge <- factor(MAE_df$edge, levels = FC_pairs_UT)
MAE_df$type <- ifelse(MAE_df$edge %in% c('V1-V2','V1-V3','V2-V3'),
                      'Visual-Visual',
                      ifelse(MAE_df$edge %in% c('V1-M','V2-M','V3-M'), 
                             'Visual-Motor',
                             'DMN-Visuomotor'))
MAE_df$type <- factor(MAE_df$type, 
                      levels = c('Visual-Visual', 'Visual-Motor', 'DMN-Visuomotor'),
                      labels = c('Visual-Visual (strong)', 'Visual-Motor (moderate)', 'DMN-Visuomotor (weak)'))
  
pdf(file.path('plots','MAE_FC_lines.pdf'), height=5, width=9)
print(ggplot(MAE_df, aes(x=algo, y=(MAE), color=abs(FC_val), group=edge)) + 
  geom_line() + geom_point(aes(shape = edge), size=2) +
  scale_shape_manual(name = 'FC Pair', values = c(15:17,0:6)) +
  scale_color_gradient2(name = 'FC Magnitude', low='black',mid='blue',high='hotpink',midpoint = 0.3) +
  facet_wrap( ~ type, scales = "free", nrow=1) +
  theme_few() + xlab('Algorithm') + ylab('Median Absolute Error'))
dev.off()

##############################################################################
### FC CREDIBLE INTERVALS

FC_true_best <- (FC_true0 + FC_true)/2

# FC estimates
df_FC <- NULL
for(ii in 1:50){
  FC_est_ii <- apply(FC_est[,,ii,], 3, function(x) x[upper.tri(x)])
  FC_LB_ii <- apply(FC_LB[,,ii,], 3, function(x) x[upper.tri(x)]); FC_LB_ii <- cbind(NA, FC_LB_ii, NA, NA) 
  FC_UB_ii <- apply(FC_UB[,,ii,], 3, function(x) x[upper.tri(x)]); FC_UB_ii <- cbind(NA, FC_UB_ii, NA, NA)
  FC_true_ii <- FC_true[,,ii][upper.tri(FC_true[,,ii])]
  FC_true0_ii <- FC_true0[,,ii][upper.tri(FC_true0[,,ii])]
  FC_true_best_ii <- FC_true_best[,,ii][upper.tri(FC_true_best[,,ii])]
  df_FC_ii <- data.frame(subject = ii, 
                         FC_est = c(FC_est_ii), 
                         FC_LB = c(FC_LB_ii),
                         FC_UB = c(FC_UB_ii),
                         FC_true = FC_true_ii,
                         FC_true0 = FC_true0_ii,
                         FC_true_best = FC_true_best_ii,
                         algo = rep(algos3, each = 10), 
                         edge = rep(1:10, times = 5),
                         row = rows[upper.tri(rows)],
                         col = cols[upper.tri(cols)])
  df_FC <- rbind(df_FC, df_FC_ii)
}

df_FC <- subset(df_FC, algo %in% c('VB1','VB2'))
df_FC$row <- factor(df_FC$row, levels = IC_names)
df_FC$col <- factor(df_FC$col, levels = IC_names)

df_FC_10 <- subset(df_FC, subject %in% 1:10)
df_FC_10_true <- df_FC_10[,c('subject','edge','row','col','FC_true')]; df_FC_10_true$truth <- 'held out data'
df_FC_10_true0 <- df_FC_10[,c('subject','edge','row','col','FC_true0')]; df_FC_10_true0$truth <- 'model data'
names(df_FC_10_true0)[5] <- 'FC_true'
df_FC_10_true <- rbind(df_FC_10_true, df_FC_10_true0)

pdf(file.path('plots','FC_CIs.pdf'), height=6, width=6)
print(ggplot(df_FC_10, aes(x = subject, color=algo)) +
  #geom_point(aes(y = FC_true), size=2, shape=8) +
  #geom_point(aes(y = FC_est)) +
  geom_errorbar(data=subset(df_FC_10, algo=='VB2'), aes(ymin = FC_LB, ymax = FC_UB)) +
  geom_errorbar(data=subset(df_FC_10, algo=='VB1'), aes(ymin = FC_LB, ymax = FC_UB)) +
  geom_point(data=df_FC_10_true, aes(x = subject, y = FC_true, shape=truth), color='black', size=0.8) +
  scale_shape_manual(values=c(4, 19)) +
  scale_color_manual(values=c('royalblue','orange')) +
  theme_few() + theme(legend.position='bottom') +
  ggtitle('FC Credible Intervals') + ylab('FC Value') + 
  ylim(-1,1) + scale_x_continuous(breaks=1:10) +
  facet_grid(row ~ col))
dev.off()

#quantify coverage 

LB_check <- array(apply(FC_LB, 4, function(x){ x < FC_true0 } ), dim = dim(FC_LB))
UB_check <- array(apply(FC_UB, 4, function(x){ x > FC_true0 } ), dim = dim(FC_UB))
CI_check <- LB_check & UB_check
CI_coverage <- apply(CI_check, c(1:2,4), mean)
print(mean(CI_coverage[,,1][upper.tri(CI_coverage[,,1])])) #VB1 0.112 
print(mean(CI_coverage[,,2][upper.tri(CI_coverage[,,2])])) #VB2 0.73 


# CI width

CI_width <- FC_UB - FC_LB
CI_width <- apply(CI_width, c(1:2,4), mean)
mean(CI_width[,,1][upper.tri(CI_width[,,1])]) #VB1 0.0141
mean(CI_width[,,2][upper.tri(CI_width[,,2])]) #VB2 0.120

colnames(comptime) <- algos2
print(colMeans(comptime)) #in seconds
# > colMeans(comptime)
# tICA         VB1         VB2          DR 
# 2.1414350   2.9944548 100.7714554   0.4830997 
print(apply(comptime, 2, sd))
# > apply(comptime, 2, sd)
# tICA        VB1        VB2         DR 
# 0.1423187  0.2898549 20.3561383  0.0149910 
colnames(numiter) <- algos
print(colMeans(numiter))
# > colMeans(numiter)
# tICA  VB1  VB2 
# 9.50 8.86 9.18

##############################################################################
### FC ACCURACY BY DURATION

load(file = file.path('results', 'by_duration','test_set_results.RData')) #FC_est, FC_true, durations

algos <- c('tICA','VB1','VB2','DR') 
num_alg <- length(algos)

durations <- c(100, 200, 400, 600)
num_dur <- length(durations)

#FC pair names
IC_names <- c('V1','V2','V3','DMN','M')
rows <- matrix(IC_names, nrow=Q, ncol=Q)
cols <- matrix(IC_names, nrow=Q, ncol=Q, byrow=TRUE)
FC_pairs <- matrix(paste0(rows, '-', cols), nrow=Q, ncol=Q)
FC_pairs_UT <- FC_pairs[upper.tri(FC_pairs)]


#MAE of FC -- using median rather than mean due to outlier subjects
MAE_FC_df <- NULL
for(aa in 1:num_alg){
  alg <- algos[aa]
  pdf(paste0('plots/MAE_FC_byduration_',algos[aa],'.pdf'), height=5, width=5.5)
  for(dd in 1:num_dur){
    MAE_ad <- sqrt(apply((FC_est[,,,aa,dd] - FC_true)^2, c(1,2), median, na.rm=TRUE))
    plot_FC(MAE_ad, zlim=c(0, 0.3001), break_by = 0.1, cor=TRUE, digits_legend = 2, cols=viridis_pal(option = 'magma')(10),
            title = paste0('MAE (', alg, '), T=',durations[dd]), box=TRUE, box_col='springgreen')
    
    MAE_FC_df_ad <- data.frame(Algorithm = alg, 
                               Duration = durations[dd],
                               MAE = MAE_ad[upper.tri(MAE_ad)],
                               FC_pair = FC_pairs_UT)
    MAE_FC_df <- rbind(MAE_FC_df, MAE_FC_df_ad)
  }
  dev.off()
}

MAE_FC_df <- MAE_FC_df %>% 
  filter(Duration >= 200) %>%
  group_by(Algorithm, Duration) %>% 
  summarize(MAE = mean(MAE))

MAE_FC_df$Algorithm <- factor(MAE_FC_df$Algorithm, 
                              levels = c('VB1','VB2','tICA','DR'),
                              labels = c('FC-tICA (VB1)','FC-tICA (VB2)','tICA','DR'))

pdf('plots/MAE_FC_byduration.pdf', width=6, height=3.5)
ggplot(MAE_FC_df, aes(x=Duration, y=MAE, color=Algorithm)) +
  geom_point() + geom_line(aes(group=Algorithm)) +
  scale_color_manual(values=c('royalblue','orange','black','gray')) +
  theme_few() + ylab('MAE of FC') + xlab('Scan Duration (Volumes)')
dev.off()

#     
