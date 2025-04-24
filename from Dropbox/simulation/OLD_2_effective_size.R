xlibrary(ciftiTools)
ciftiTools.setOption('wb_path', '/Applications')
library(coda)
library(viridis)

main_dir <- '~/Dropbox/RESEARCH/FCTemplateICA/simulation/'
setwd(main_dir)
source('sim_funs.R')

load(file='TCs_sim.RData')
TCs <- TC_both[[1]]
subjICs <- readRDS('subjICs.RDS')

load(paste0("sim1/results/testsubj1.RData"))
nsamp <- dim(result1[[1]])[3]
nvox <- dim(result1[[1]])[1]
nICs <- dim(result1[[1]])[2]

FC_within <- FC_within_oracle <- matrix(0, nrow=nICs, ncol=nICs)
S_within <- S_within_oracle <- matrix(0, nrow=nvox, ncol=nICs)

nsubj <- 50
for(ii in 1:nsubj){

  print(ii)
  # Look at the effective size of the first (S) and sixth (FC) result element
  # and calculate the standard error
  if(ii > 1) load(paste0("sim1/results/testsubj",ii,".RData"))

  # Effective sample size for S
  neff_S_ii <- apply(result1[[1]], 1:2, coda::effectiveSize)
  if(ii==1) neff_S <- neff_S_ii else neff_S <- neff_S + neff_S_ii

  # Credible interval and coverage for S
  CI_S <- apply(result1[[1]],c(1,2),quantile, probs=c(0.025, 0.975))
  CI_S_LB <- CI_S[1,,]
  CI_S_UB <- CI_S[2,,]
  S_true <- subjICs[,,500+ii]
  S_within_ii <- (S_true >= CI_S_LB) & (S_true <= CI_S_UB)
  S_within <- S_within + S_within_ii

  # using oracle FC template
  CI_S <- apply(result2[[1]],c(1,2),quantile, probs=c(0.025, 0.975))
  CI_S_LB <- CI_S[1,,]
  CI_S_UB <- CI_S[2,,]
  S_within_ii <- (S_true >= CI_S_LB) & (S_true <= CI_S_UB)
  S_within_oracle <- S_within_oracle + S_within_ii

  # Effective sample size for FC
  neff_FC_ii <- apply(result1[[6]],1,coda::effectiveSize)
  neff_FC_ii <- matrix(neff_FC_ii, nICs, nICs)
  if(ii==1) neff_FC <- neff_FC_ii else neff_FC <- neff_FC + neff_FC_ii

  # Credible interval and coverage for FC
  FC_true <- cor(TCs[,,500+ii])
  CI_FC <- apply(result1[[6]],1,quantile, probs=c(0.025, 0.975))
  CI_FC_LB <- matrix(CI_FC[1,], nICs, nICs)
  CI_FC_UB <- matrix(CI_FC[2,], nICs, nICs)
  FC_within_ii <- (FC_true >= CI_FC_LB) & (FC_true <= CI_FC_UB)
  FC_within <- FC_within + FC_within_ii

  # using oracle FC template
  CI_FC <- apply(result2[[6]],1,quantile, probs=c(0.025, 0.975))
  CI_FC_LB <- matrix(CI_FC[1,], nICs, nICs)
  CI_FC_UB <- matrix(CI_FC[2,], nICs, nICs)
  FC_within_ii <- (FC_true >= CI_FC_LB) & (FC_true <= CI_FC_UB)
  FC_within_oracle <- FC_within_oracle + FC_within_ii

}

#plot efficiency for S, averaged across subjects
effic_S <- neff_S / nsubj / nsamp # There were 950 samples taken from the posterior
effic_S <- newdata_xifti(result1$subjICmean, effic_S)
neff_S <- newdata_xifti(result1$subjICmean, neff_S / nsubj)
plot(effic_S, zlim=c(0.5,1), idx=1:5, fname='sim1/images/effic')
plot(neff_S, zlim=c(5000,10000), idx=1:5, fname='sim1/images/neff')

#plot coverage rate for S, averaged across subjects
S_within <- S_within / nsubj
S_within_oracle <- S_within_oracle / nsubj
summary(S_within_oracle - S_within)
S_within <- newdata_xifti(result1$subjICmean, S_within)
S_within_oracle <- newdata_xifti(result1$subjICmean, S_within_oracle)
plot(S_within, zlim=c(0.5,1), idx=1:5, fname='sim1/images/coverage')
plot(S_within_oracle, zlim=c(0.5,1), idx=1:5, fname='sim1/images/coverage_oracle')

#plot FC effiency, averaged across subjects
effic_FC <- neff_FC / nsubj / nsamp
pdf('sim1/plots/FC_effic.pdf', height=5, width=5.5)
plot_FC(effic_FC, zlim=c(0,0.3), break_by=0.1, cols = viridis(12), title = 'Efficiency', cor=TRUE)
plot_FC(neff_FC / nsubj, zlim=c(200,300), break_by=50, cols = viridis(12), title = 'Effective Sample Size', cor=TRUE)
dev.off()

#plot FC coverage rate, averaged across subjects
FC_within <- FC_within / nsubj
FC_within_oracle <- FC_within_oracle / nsubj
pdf('sim1/plots/FC_coverage.pdf', height=5, width=5.5)
plot_FC(FC_within, zlim=c(0,1), break_by=0.1, cols = viridis(12), title = 'FC Coverage Rate', cor=TRUE)
plot_FC(FC_within_oracle, zlim=c(0,1), break_by=0.1, cols = viridis(12), title = 'FC Coverage Rate, Oracle', cor=TRUE)
dev.off()

#plot sample chains of FC for one subject
library(ggthemes)
library(ggplot2)
library(reshape2)
UT_entries <- as.vector(upper.tri(diag(5)))
pair_names <- c('IC1-IC2','IC1-IC3','IC2-IC3','IC1-IC4','IC2-IC4','IC3-IC4','IC1-IC5','IC2-IC5','IC3-IC5','IC4-IC5')
FC_samples_ii <- as.data.frame(result2$FC_samples[UT_entries,1:1000])
FC_samples_ii$PairName <- pair_names
FC_samples_long <- melt(FC_samples_ii, id.vars = 'PairName')
names(FC_samples_long)[2] <- 'Sample'
FC_samples_long$Sample <- as.numeric(FC_samples_long$Sample)
FC_true_long <- data.frame(Pair = 1:10, PairName = pair_names, Truth = FC_true[upper.tri(FC_true)])
pdf('sim1/plots/FC_chains.pdf', height=7, width=6)
ggplot(FC_samples_long, aes(x=Sample, y=value, group=PairName, color=PairName)) +
  geom_line() +
  geom_hline(data=FC_true_long, aes(yintercept=Truth, color=PairName)) +
  scale_color_brewer(palette='Paired') +
  theme_few() + theme(legend.position='bottom', legend.title=element_blank())
dev.off()

#plot 10,000th sample of A for each IC
A_samples1 <- result1$A_samples[,,5000]
A_samples2 <- result1$A_samples[,,10000]
A_true_ii <- TCs[,,500+ii]
A_samples1_long <- melt(A_samples1, varnames = c('Time','IC'))
A_samples2_long <- melt(A_samples2, varnames = c('Time','IC'))
A_true_long <- melt(A_true_ii, varnames = c('Time','IC'))
ggplot(A_true_long, aes(x=Time, y=value)) + geom_line() +
  geom_line(data=A_samples1_long, aes(x=Time, y=value), color='red', alpha=0.5) +
  geom_line(data=A_samples2_long, aes(x=Time, y=value), color='blue', alpha=0.5) +
  facet_grid(IC ~ ., labeller='label_both') + xlim(0,500) +
  theme_few()



