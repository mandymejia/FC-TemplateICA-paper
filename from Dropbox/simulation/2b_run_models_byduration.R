### PART 2: VARY SCAN DURATION ---------------------------------------------
main_dir <- '~/Dropbox/RESEARCH/FCTemplateICA/simulation/'
#main_dir <- '~/Dropbox (Brown)/FCTemplateICA/simulation/'
setwd(main_dir)
source('0_setup.R')

##############################################################################
### LOAD IN PREVIOUSLY CREATED RESULTS

load(file = 'GICA.RData') #ICs, Q, GICA, N

subjICs <- readRDS('subjICs.RDS')

template <- readRDS('templates/template.RDS')

n <- 500 + 50 #500 training subjects (use to re-estimate template) + 50 test subjects
BOLD_fnames <- paste0('data/subj',1:n,'.dtseries.nii')

##############################################################################
### ANALYZE TEST SUBJECTS

inds_test <- 501:550
BOLD_fnames_test <- BOLD_fnames[inds_test]

#vary duration of timeseries
durations <- c(100, 200, 400)
num_dur <- length(durations)

algos <- c('tICA','VB1','VB2') 
num_alg <- length(algos)

#for saving subject-level results
results_list_dur <- vector('list', num_dur)
names(results_list_dur) <- durations
results_list_alg <- vector('list', num_alg)
names(results_list_alg) <- algos

for(ii in 1:50){
  print(paste0('~~~~~~~~~~~~~~~~ SUBJECT ',ii,' ~~~~~~~~~~~~~~~~'))

  results_ii <- results_list_dur
  comptime_ii <- matrix(NA, num_alg, num_dur)

  #loop over varying durations
  for(dd in 1:num_dur){
    dur <- durations[dd]
    print(paste0('~~~~~~ duration: ',dur,' ~~~~~~~~~'))

    results_ii[[dd]] <- results_list_alg
    
    #loop over algorithms
    for(aa in 2:num_alg){
      
      print(alg <- algos[aa])
      
      # if(alg=='tICA') {
      #   template_aa <- template_noFC
      #   method_FC_aa <- "none"
      #} else {
        # template_aa <- template
        # method_FC_aa <- alg
      #}
      
      print(system.time(result_ad <- templateICA(BOLD = BOLD_fnames_test[ii],
                                 template = template,
                                 scale = 'none',
                                 Q2 = 0,
                                 epsilon = 0.001, 
                                 brainstructures = 'left', 
                                 verbose=FALSE,
                                 method_FC = alg,
                                 reduce_dim = FALSE,
                                 time_inds = 1:dur)))

      results_ii[[dd]][[aa]] <- result_ad
      #comptime_ii[aa,dd] <- time_ad[3]

      } #end loop over algorithms
    
    } #end loop over durations
  
  save(results_ii, file=paste0('results/by_duration/testsubj',ii,'_bydur.RData'))

} #end loop over subjects




