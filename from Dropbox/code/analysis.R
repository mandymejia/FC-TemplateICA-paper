#roxygen2::roxygenize('~/Documents/Github/templateICAr/')
#devtools::install_github('mandymejia/templateICAr', ref='8.0')
library(templateICAr) #0.8.5
library(fMRItools) #0.5.0
library(ciftiTools) #0.16.0
ciftiTools.setOption('wb_path', '/Applications/')
library(dplyr) #1.1.4
library(viridis) #0.6.5
library(RColorBrewer) #1.1-3
library(matrixStats) #1.3.0
library(ggplot2) #3.5.1
library(ggthemes) #5.1.0

main_dir <- '~/Dropbox/RESEARCH/FCTemplateICA'
setwd(main_dir)

#need to transfer the data from Carbonate to MacPro? (one time only)
run_carbonate <- FALSE

#need to assign ICs to networks?
run_ordering <- FALSE

#need to estimate template? (one time only)
run_template <- FALSE

#make images with ciftiTools
make_images <- FALSE

source('code/functions.R') #plot_FC()

###############################################################
## COLLECT AND TRANSFER THE DATA
###############################################################

fname_ts <- '_Atlas_hp2000_clean.dtseries.nii'
session_names <- c('rfMRI_REST1_LR','rfMRI_REST2_LR')
cifti_fnames<- file.path(session_names, paste0(session_names, fname_ts))

if(run_carbonate){
  data_dir <- '/N/dcwan/projects/hcp' #as of June 2023 this data was moved to /N/project/hcp_dcwan but will eventually be retired
  data_dir_local <- '/Volumes/Lab_Data_Drive/data/HCP_Resting_State'

  ### RANDOMLY SELECT SUBJECTS FOR TRAINING AND TESTING

  ## RUN THIS ON CARBONATE
  library(stringr) #str_which
  #select training and test subjects
  all_subjects <- list.files(data_dir) #each subject is a directory
  all_subjects <- all_subjects[nchar(all_subjects)==6] #subject IDs have 6 characters
  all_subjects <- all_subjects[str_which(all_subjects, "[0-9]+")] #subject IDs are numeric
  set.seed(782439)
  samp_subjects <- sample(all_subjects, 501, replace = FALSE)
  save(samp_subjects, file='subjects.RData')
  load(file='subjects.RData') #samp_subjects

  ### TRANSFER DATA FROM CARBONATE TO MAC PRO
  for(ii in 1:501){
    print(paste0('Subject ', ii))
    subi <- samp_subjects[ii]
    for(jj in 1:2){
      print(jj)
      dir_ij <- file.path(data_dir,subi,cifti_path[jj])
      from_ij <- file.path(dir_ij, cifti_fnames[jj])
      to_ij <- file.path(data_dir_local, paste0(subi,'_',cifti_fnames[jj]))
      cmd <- paste0('scp afmejia@carbonate.uits.iu.edu:',from_ij,' ',to_ij)
      system(cmd)
    }
  }

  #count number of files for each subject
  numfiles <- table(substr(list.files(data_dir_local), 1, 6))
  #use subjects who have two scans
  subjects_2scans <- names(numfiles)[which(numfiles==2)]

  template_subjects <- subjects_2scans[1:(length(subjects_2scans)-1)]
  test_subjects <- subjects_2scans[length(subjects_2scans)]
  save(template_subjects, test_subjects, file='subjects2.RData')

}

nQ <- 25

###############################################################
# ASSIGN ICs to CORTICAL NETWORKS
###############################################################

GICA_fname <- paste0('data/melodic_IC_',nQ,'.dscalar.nii')

if(run_ordering){

  GICA <- read_cifti(GICA_fname, brainstructures = 'all')
  mask <- c(GICA$meta$cortex$medial_wall_mask$left, GICA$meta$cortex$medial_wall_mask$right)
  #visualize GICA maps
  GICA$meta$cifti$names <- paste0('IC',GICA$meta$cifti$names)
  if(make_images) plot(GICA, idx=1:nQ, fname = file.path(main_dir,'templates',paste0('GICA',nQ),'GICA'))

  #read in cortical parcellation
  Yeo17 <- ciftiTools::load_parc('Yeo_17')
  if(make_images) plot(Yeo17, fname = file.path(main_dir,'templates','Yeo17'))
  Yeo17_vec <- as.matrix(Yeo17)[mask]
  #extract Yeo network names
  Yeo17_names <- gsub('_.+','',gsub('17Networks_.H_','',rownames(Yeo17$meta$cifti$labels$parcels), fixed = FALSE))
  Yeo17_names <- gsub('[ABC]$','',Yeo17_names)
  Yeo17_names_df <- data.frame(parcel_num = Yeo17$meta$cifti$labels$parcels$Key, parcel_name = Yeo17_names)
  Yeo17_colors_df <- Yeo17$meta$cifti$labels$parcels
  Yeo17_colors_df$names <- Yeo17_names
  Yeo17_colors_df$colorRGB <- rgb(red=Yeo17_colors_df$Red, green = Yeo17_colors_df$Green, blue = Yeo17_colors_df$Blue)

  #read in subcortical parcellation
  xii <- read_cifti(ciftiTools.files()$cifti["dscalar_ones"], brainstructures="all")
  subcort_names <- xii$meta$subcort$labels
  subcort_names <- gsub('-[LR]$','',xii$meta$subcort$labels)
  subcort_names_unique <- unique(subcort_names)
  subcort_names_df <- data.frame(parcel_num = 1000 + 1:length(subcort_names_unique), #add 1k to avoid overlap with Yeo parcel nums
                                 parcel_name = subcort_names_unique)
  subcort_vec_df <- data.frame(parcel_name = subcort_names)
  subcort_vec_df <- left_join(subcort_vec_df, subcort_names_df)

  xii$data$subcort[,1] <- as.numeric(as.factor(subcort_names))
  xii <- ciftiTools:::convert_to_dlabel(xii)
  rownames(xii$meta$cifti$labels$ones) <- levels(as.factor(subcort_names))
  xii$meta$cifti$labels$ones[1,2:4] <- rep(0.3,3) #recolor accumbens to dark gray
  view_xifti_volume(xii)

  #combine cortical and subcortical parcels
  parcel_vec <- c(Yeo17_vec, subcort_vec_df$parcel_num) #contains Yeo sub-parcels that will be combined later
  parcel_names_df <- rbind(Yeo17_names_df, subcort_names_df) #contains Yeo sub-parcels that will be combined later
  parcel_names_df$parcel_size <- table(parcel_vec)
  #all.equal(as.character(parcel_names_df$parcel_num), names(table(parcel_vec))) #TRUE

  network_ICs <- rep(NA, nQ)
  for(q in 1:nQ){
    #threshold IC using SD after median-centering
    IC_q <- as.matrix(GICA)[,q]
    IC_q <- IC_q - median(IC_q)
    sd_q <- sd(IC_q)
    if(nQ==25) thr <- sd_q
    if(nQ==100) thr <- 0.5*sd_q
    IC_q_thr <- (abs(IC_q) > thr)

    #calculate the number of IC vertices/voxels overlapping with each parcel
    parcel_num <- IC_q_thr * parcel_vec
    parcel_df <- data.frame(parcel_num = parcel_num, values = IC_q_thr * abs(IC_q)) #second column is the intensity at each vertex (0 if thresholded)
    parcel_df_summ <- parcel_df %>% group_by(parcel_num) %>% summarize(count = n(), total = sum(values))
    parcel_df_summ <- left_join(parcel_df_summ, parcel_names_df, by = 'parcel_num') #bring in parcel name and size
    parcel_df_summ <- parcel_df_summ %>% group_by(parcel_name) %>% summarize(total = sum(total), count = sum(count), parcel_size = sum(parcel_size)) #group parcels (e.g. DMN1, DMN2)
    parcel_df_summ$score <- parcel_df_summ$total/sqrt(parcel_df_summ$parcel_size) #sum up the values and divide by parcel size
    parcel_df_summ$score0 <- parcel_df_summ$total/(parcel_df_summ$parcel_size) #sum up the values and divide by parcel size
    network_ICs[q] <- parcel_df_summ$parcel_name[which.max(parcel_df_summ$score)]
    network_ICs0 <- parcel_df_summ$parcel_name[which.max(parcel_df_summ$score0)]
    if(network_ICs0=="Cerebellum" | network_ICs[q] == "Cerebellum"){
      print(q)
      print(network_ICs0)
      print(network_ICs[q])
    }

    #consider first runner up for Cerebellar assignments (results in one reassignment for IC 50)
    if(network_ICs[q]=='Cerebellum'){
      parcel_df_summ2 <- parcel_df_summ
      parcel_df_summ2$score[parcel_df_summ2$parcel_name==network_ICs[q]] <- 0
      diff <- max(parcel_df_summ2$score)/max(parcel_df_summ$score)
      if(diff > 0.75) network_ICs[q] <- parcel_df_summ2$parcel_name[which.max(parcel_df_summ2$score)]
    }

    # #old methodology (saved below since we used this ordering for template estimation and model fitting)
    # counts_df <- data.frame(table(parcel_num)[-1]) #count number of vertices per parcel, excluding IC background vertices
    # counts_df$parcel_num <- as.numeric(as.character(counts_df$parcel_num)) #convert factor to numeric
    # counts_df <- left_join(counts_df, parcel_names_df, by = 'parcel_num') #bring in parcel name (grouped)
    # counts_df <- counts_df %>% group_by(parcel_name) %>% summarize(count = sum(Freq), parcel_size = sum(parcel_size))
    # #compute what percentage of the parcel the overlap with the IC represents, identify maximum
    # counts_df$pct <- counts_df$count/counts_df$parcel_size*100
    # network_ICs[q] <- counts_df$parcel_name[which.max(counts_df$pct)]
    # network_ICs_pct[q] <- counts_df$pct[which.max(counts_df$pct)]
    # #re-classify this one based on visual inspection <-- this was fixed with the new methodology!
    #Note: this was done after estimating the templates, so we have to fix the indexing of the final results with `inds_fix` below
    # if(nQ==25 & q==23) network_ICs[q] <- 'Cerebellum'

    #visualize thresholded GICA map
    IC_q[IC_q_thr==FALSE] <- NA
    xii_q <- newdata_xifti(GICA, IC_q)
    title_q <- paste0('IC ',q,' (',network_ICs[q],')')
    if(make_images) plot(xii_q, title=title_q, fname=file.path('templates',paste0('GICA',nQ),paste0('GICA_IC',q,'_thr')), legend_fname=NULL)
  }
  network_ICs_df <- data.frame(IC = 1:nQ,network = network_ICs)

  sort(table(network_ICs), decreasing = TRUE)
  # nQ = 25:
  #   Cerebellum        Cont     Default    DorsAttn SalVentAttn      SomMot     TempPar     VisCent     VisPeri
  #.     3           5           3           1           2           3           1           5           2


  # nQ = 100:
  # Cerebellum     Default     VisCent        Cont      SomMot  Brain Stem    DorsAttn SalVentAttn     TempPar     VisPeri     Putamen     Caudate Hippocampus      Limbic    Thalamus
  #.  27          13          13          12           7           6           6           4           3           3           2           1           1           1           1

  # #Bad ICs among 100
  # #badinds <- c(74,96:100) #Identified these because they exhibited high between subject variance. Inspection of GICA maps showed these to contain high levels of noise. Later we realized many of these are subcortical ICs, so they appeared noisy in the cortex.
  # #badinds <- which(network_ICs_pct < 25) #for nQ=100, IC 74 has 10% overlap with brainstem
  # if(nQ==100) badinds <- 74 else badinds <- c()

  #reorder ICs in order of parcels
  parcel_order <- data.frame(network = c(unique(Yeo17_names[-1]), subcort_names_unique))
  network_ICs_df_ordered <- right_join(parcel_order, network_ICs_df)
  IC_order <- network_ICs_df_ordered$IC
  IC_order <- setdiff(IC_order, badinds)
  save(IC_order, badinds, network_ICs_df_ordered, network_ICs, file=paste0('data/IC_networks_ordered_',nQ,'.RData'))

  #group subcortical regions
  if(length(badinds) > 0) {
    badrow <- which(network_ICs_df_ordered$IC %in% badinds)
    network_ICs_df_ordered <- network_ICs_df_ordered[-badrow,]
  }
  networks <- network_ICs_df_ordered$network
  subcort <- c('Amygdala','Hippocampus','Accumbens','Putamen','Caudate','Thalamus')
  network_ICs_df_ordered$network[network_ICs_df_ordered$network %in% subcort] <- 'Subcortical'

  #make short names for networks
  networks_names_df <- data.frame(network = c('VisCent','VisPeri','SomMot','DorsAttn','SalVentAttn','Limbic','Cont','Default','TempPar','Brain Stem','Cerebellum','Subcortical'),
                                  network2 = c('V','V','M','A','A','L','C','D','TP','B','CB','SUB'))
  network_ICs_df_ordered <- left_join(network_ICs_df_ordered, networks_names_df)

  #find change in network labels (where to draw lines in FC matrices)
  networks <- network_ICs_df_ordered$network2
  change <- which(!(networks[1:(nL-1)] == networks[2:nL]))
  labels <- unique(networks)

  #save(inds_orig, inds, reorder, networks, change, labels, networks_names_df, network_ICs_df_ordered, file=paste0('data/inds_',nQ,'.RData'))
  save(networks, change, labels, networks_names_df, network_ICs_df_ordered, file=paste0('data/inds_',nQ,'.RData'))
} # end run_orderin

# updated assignments with weighted overlap in Sep 2023
load(file=paste0('data/IC_networks_ordered_',nQ,'.RData')) #IC_order, badinds, network_ICs_df_ordered
nL <- length(IC_order)
cex <- 0.8 #for FC matrix network labels
#if(nQ==100) cex <- 0.6
load(file=paste0('data/inds_',nQ,'.RData')) #networks, change, labels, networks_names_df, network_ICs_df_ordered

###############################################################
## ESTIMATE TEMPLATES
###############################################################

if(run_template){

  # #re-assign the subjects to hold out 100 test subjects
  # load(file='subjects2.RData')
  # subjects <- c(template_subjects, test_subjects)
  # set.seed(2384723)
  # test_subjects <- sort(sample(subjects, 100)) #100 test subjects
  # template_subjects <- sort(setdiff(subjects, test_subjects)) #362 template subjects
  # save(template_subjects, test_subjects, file='subjects3.RData')
  load(file='subjects3.RData')

  data_dir <- '/Volumes/Lab_Data_Drive/data/HCP_Resting_State'

  #file path to all training subjects (test and retest)
  cifti_fullnames1 <- file.path(data_dir, paste0(template_subjects,'/',cifti_fnames[1]))
  cifti_fullnames2 <- file.path(data_dir, paste0(template_subjects,'/',cifti_fnames[2]))

  #8-12 hours for 25-100 ICs
  system.time(template <- estimate_template(BOLD = cifti_fullnames1,
                                BOLD2 = cifti_fullnames2,
                                GICA = GICA_fname,
                                inds = IC_order,
                                scale = "local",
                                scale_sm_FWHM = 2, #default
                                FC = TRUE,
                                #hpf=0, #don't need because data is already filtered
                                keep_FC = TRUE,
                                brainstructures = 'all'))

  # #for separate exploratory analysis for Cholesky-based FC template
  # tempFC <- template$template
  # tempFC$mean <- tempFC$varUB <- tempFC$varNN <- NULL
  # FC_vals0 <- template$FC
  # FC_vals <- template$FC_chol #this was just the Cholesky decompositions of the FC matrices
  # save(tempFC, FC_vals, FC_vals0, file='template_chol_TEST.rds')

  saveRDS(template, file.path(main_dir,'templates',paste0('template_',nQ,'.rds')))

  ###############################################################
  ## VISUALIZE TEMPLATES
  ###############################################################

  template0 <- template
  #template0$template$FC_Chol <- NULL
  template_exp <- export_template(template0)
  empirical_FC_mean <- template_exp$FC$mean_empirical
  empirical_FC_var <- template_exp$FC$var_empirical
  template_FC_mean <- template_exp$FC$mean
  template_FC_var <- template_exp$FC$var
  template_FC_Chol_mean <- template0$template$FC_Chol$FC_samp_mean
  template_FC_Chol_var <- template0$template$FC_Chol$FC_samp_var


  #VISUALIZE IW FC TEMPLATE MEAN AND SD
  pdf(file.path(main_dir,'plots',paste0('GICA',nQ),'FC_template_IW.pdf'), width=5, height=5)
  diag(template_FC_var) <- NA
  plot_FC(template_FC_mean, zlim=c(-0.6, 0.6), break_by = 0.2, title = "IW FC Template Mean", cor=TRUE, lines = change, labels = labels, lwd_lines=1.5, cex=cex)
  plot_FC(sqrt(template_FC_var), zlim=c(0.1, 0.4), break_by = 0.05, digits_legend = 2, cols = viridis(12, option='B'), title = "IW FC Template SD", cor=FALSE, lines = change, labels=labels, col_lines = 'white', lwd_lines=1.5, cex=cex)
  dev.off()

  #VISUALIZE Cholesky FC TEMPLATE MEAN AND SD
  pdf(file.path(main_dir,'plots',paste0('GICA',nQ),'FC_template_Chol.pdf'), width=5, height=5)
  diag(template_FC_Chol_var) <- NA
  plot_FC(template_FC_Chol_mean, zlim=c(-0.6, 0.6), break_by = 0.2, title = "Cholesky FC Template Mean", cor=TRUE, lines = change, labels = labels, lwd_lines=1.5, cex=cex)
  plot_FC(sqrt(template_FC_Chol_var), zlim=c(0.1, 0.4), break_by = 0.05, digits_legend = 2, cols = viridis(12, option='B'), title = "Cholesky FC Template SD", cor=FALSE, lines = change, labels=labels, col_lines = 'white', lwd_lines=1.5, cex=cex)
  dev.off()

  #VISUALIZE EMPIRICAL MEAN AND SD OF FC
  pdf(file.path(main_dir,'plots',paste0('GICA',nQ),'FC_empirical.pdf'), width=5, height=5)
  diag(empirical_FC_var) <- NA
  plot_FC(empirical_FC_mean, zlim=c(-0.6, 0.6), break_by = 0.2, title = "FC Empirical Mean", cor=TRUE, lines = change, labels = labels, lwd_lines=1.5)
  plot_FC(sqrt(empirical_FC_var), zlim=c(0.1, 0.4), break_by = 0.05, digits_legend = 2, cols = viridis(12, option='B'), title = "FC Empirical SD", cor=FALSE, lines = change, labels = labels, col_lines = 'white', lwd_lines=1.5)
  dev.off()

  #VISUALIZE TEMPLATE MEAN AND SD FOR SPATIAL ICS
  template_exp$mean$meta$cifti$names <- paste0('IC',IC_order)
  template_exp$var$meta$cifti$names <- paste0('IC',IC_order)
  if(nQ==25) lim <- 0.2
  #if(nQ==100) lim <- 0.1
  if(make_images){
    view_xifti_surface(template_exp$mean, idx=1:nL, zlim=c(-lim, lim), fname = file.path('templates',paste0('template',nQ),'template_mean'))
    view_xifti_surface(template_exp$var, idx=1:nL, zlim=c(0,lim^2), fname = file.path('templates',paste0('template',nQ),'template_var'), color_mode='sequential')
  }
}
load(file='subjects3.RData')
template <- readRDS(file.path(main_dir,'templates',paste0('template_',nQ,'.rds')))

###############################################################
## ANALYZE TEST SUBJECTS
###############################################################

# Analysis plan:
# Run 100 subjects and analyze performance of FC, spatial maps, and activation maps.
# Visualize results for first 5 subjects

data_dir <- '/Volumes/Lab_Data_Drive/data/HCP_Resting_State'
result_dir <- file.path('results',paste0('GICA',nQ))
template_exp <- export_template(template)
#determine thresholds for activations
sd_Q <- sqrt(colVars(as.matrix(template_exp$mean))) #thresholds to use in activations (x2, x3, ...)
thr <- sd_Q*2

# Visualize Activation Thresholds
template_exp$mean$meta$cifti$names <- paste0('IC',IC_order)
if(make_images){
  thr_mat <- matrix(thr, nrow=nrow(template_exp$mean), ncol=ncol(template_exp$mean), byrow=TRUE)
  template_mean <- as.matrix(template_exp$mean)
  template_mean_thr <- ( abs(template_mean) > thr_mat ) * template_mean
  lim <- 0.2
  view_xifti_volume(newdata_xifti(template_exp$mean, template_mean_thr), idx=1:nL, zlim=c(-lim,lim),
                    fname = file.path('templates',paste0('template',nQ),paste0('template_mean_thr_IC',IC_order,'_sub')))
  template_mean_thr[template_mean_thr == 0] <- NA
  view_xifti_surface(newdata_xifti(template_exp$mean, template_mean_thr), idx=1:nL, zlim=c(-lim,lim),
                     fname = file.path('templates',paste0('template',nQ),'template_mean_thr_IC'))
}

nV <- 18792 #actual resolution after resampling
nI <- 100
nJ <- 2

algos <- c('VB1','VB2') #also DR and tICA but saved within other algorithms

comptime <- array(NA, dim=c(nI, nJ, 4)) # 4 algorithms (VB1, VB2, tICA, DR)
FC <- array(NA, dim=c(nQ, nQ, nI, nJ, 5)) # 5 algorithms including DR & DR2 (using GICA vs using S)
S <- array(NA, dim=c(nV, nQ, nI, nJ, 4)) # 4 algorithms including DR
SE <- array(NA, dim=c(nV, nQ, nI, nJ, 3)) # 3 algorithms

#comptime <- readRDS(file=file.path(result_dir,'comptime_VB2.rds'))
#FC <- readRDS(file=file.path(result_dir,'FC_VB2.rds'))
#S <- readRDS(file=file.path(result_dir,'S_VB2.rds'))
#SE <- readRDS(file=file.path(result_dir,'SE_VB2.rds'))

for(i in 1:5){
  print(paste0('~~~~~~~~~~~~~~~~ SUBJECT ',i,' ~~~~~~~~~~~~~~~~'))

  for(j in 1:nJ){
    print(paste0('~~~~~~~~~~~~~~~~ SESSION ',j,' ~~~~~~~~~~~~~~~~'))
    cifti_ij <- file.path(data_dir, paste0(test_subjects[i],'/',cifti_fnames[j]))

    for(aa in 1:length(algos)){
      algo <- algos[aa]
      print(paste0('~~~~~~~~~~~~~~~~ ALGORITHM: ',algo,' ~~~~~~~~~~~~~~~~'))

      #20-25 min for VB1, 4 hours for VB2
      time_ijq <- system.time(
        result_ijq <- suppressMessages(templateICA(cifti_ij,
                                template,
                                method_FC = algo,
                                scale='local',
                                scale_sm_FWHM = 2, #default
                                Q2 = 200, #PESEL tends to estimate a large number of nuisance components, set to 200 to save time
                                hpf = 0, #data already temporally filtered
                                resamp_res=10000,
                                brainstructures=c('left','right'),
                                reduce_dim=FALSE,
                                usePar=TRUE,
                                verbose=FALSE)))
      print(time_ijq)

      if(i <= 5) save(result_ijq, time_ijq, file=file.path(result_dir,paste0('result_subj',i,'_sess',j,'_',algo,'.RData')))
      #if(i <= 5) load(file=file.path(result_dir,paste0('result_subj',i,'_sess',j,'_',algo,'.RData')))
      saveRDS(result_ijq$FC, file=file.path(result_dir,paste0('FC_subj',i,'_sess',j,'_',algo,'.rds')))

      xii <- result_ijq$subjICmean
      temp_mean <- result_ijq$template_mean

      #save computation time in seconds
      comptime[i,j,aa] <- result_ijq$comptime["FC-tICA"]
      if(aa==1){
        comptime[i,j,3] <- result_ijq$comptime["tICA"]
        comptime[i,j,4] <- result_ijq$comptime["DR"]
      }

      #save FC estimates
      FC[,,i,j,aa] <- cor(result_ijq$A) #FC-tICA
      if(aa==1){
        FC[,,i,j,3] <- cor(result_ijq$result_tICA$A) #tICA
        FC[,,i,j,4] <- cor(result_ijq$result_DR$A) #DR, based on GICA
        FC[,,i,j,5] <- cor(result_ijq$result_DR$A2) #DR, based on S ("DR+")
      }

      #save IC estimates
      S[,,i,j,aa] <- as.matrix(result_ijq$subjICmean) #FC-tICA
      if(aa==1){
        S[,,i,j,3] <- result_ijq$result_tICA$subjICmean #tICA
        S[,,i,j,4] <- t(result_ijq$result_DR$S) #DR
      }

      #save IC standard errors
      SE[,,i,j,aa] <- as.matrix(result_ijq$subjICse) #FC-tICA
      if(aa==1){
        SE[,,i,j,3] <- result_ijq$result_tICA$subjICse #tICA
      }

      saveRDS(comptime, file=file.path(result_dir,'comptime.rds'))
      saveRDS(FC, file=file.path(result_dir,'FC.rds'))
      saveRDS(S, file=file.path(result_dir,'S.rds'))
      saveRDS(SE, file=file.path(result_dir,'SE.rds'))

   } #end loop over algorithms

  } #end loop over visits

} #end loop over subjects

comptime <- readRDS(file=file.path(result_dir,'comptime.rds'))
FC <- readRDS(file=file.path(result_dir,'FC.rds'))
S <- readRDS(file=file.path(result_dir,'S.rds'))
SE <- readRDS(file=file.path(result_dir,'SE.rds'))


#--------------------------
# COMPUTATION TIME


comptime_VB2 <- readRDS(file=file.path(result_dir,'comptime_VB2.rds'))
comptime[,,2] <- comptime_VB2[,,2]
# saveRDS(comptime, file=file.path(result_dir,'comptime.rds'))

apply(comptime, 3, mean)/60
# VB1, VB2, tICA, DR:
# 23.45712  278.23319  16.64218  9.36696 #note: VB2 based on incomplete results


#--------------------------
# FC MATRICES

FC_VB2 <- readRDS(file=file.path(result_dir,'FC_VB2.rds'))
FC[,,,,2] <- FC_VB2[,,,,2]
# saveRDS(FC, file=file.path(result_dir,'FC.rds'))

### Visualize 3 subjects
for(i in 1:3){

  for(j in 1:2){

    pdf(file.path(main_dir,'plots',paste0('GICA',nQ),paste0('FC_subj',i,'_sess',j,'.pdf')), width=5, height=5)
    plot_FC(FC[,,i,j,1], zlim=c(-0.8, 0.8), break_by = 0.2, title = "FC Template ICA (VB1)", cor=TRUE, lines = change, labels = labels, lwd_lines=1.5, cex=cex)
    plot_FC(FC[,,i,j,2], zlim=c(-0.8, 0.8), break_by = 0.2, title = "FC Template ICA (VB2)", cor=TRUE, lines = change, labels = labels, lwd_lines=1.5, cex=cex)
    plot_FC(FC[,,i,j,3], zlim=c(-0.8, 0.8), break_by = 0.2, title = "Template ICA", cor=TRUE, lines = change, labels = labels, lwd_lines=1.5, cex=cex)
    plot_FC(FC[,,i,j,4], zlim=c(-0.8, 0.8), break_by = 0.2, title = "DR", cor=TRUE, lines = change, labels = labels, lwd_lines=1.5, cex=cex)
    plot_FC(FC[,,i,j,5], zlim=c(-0.8, 0.8), break_by = 0.2, title = "DR+", cor=TRUE, lines = change, labels = labels, lwd_lines=1.5, cex=cex)
    dev.off()

  }
}

### Distance from population mean (less distance reflects shrinkage)

FC_pop_mean <- template_exp$FC$mean_empirical
FC_dist <- apply(FC, 3:5, function(x) (x - FC_pop_mean)^2)
FC_dist_mean <- apply(FC_dist, c(1,4), mean)
FC_dist_mean <- array(FC_dist_mean, dim = c(nQ, nQ, 5))

pdf(file.path(main_dir,'plots','GICA25',paste0('FC_dist.pdf')), width=5, height=5)
cols <- viridis_pal(option = 'A')(7)
plot_FC(sqrt(FC_dist_mean[,,1]), zlim=c(0, 0.4), break_by = 0.2, title = "FC Template ICA (VB1)", cor=TRUE, lines = change, labels = labels, lwd_lines=1.5, cols=cols)#, col_lines = 'white')
plot_FC(sqrt(FC_dist_mean[,,2]), zlim=c(0, 0.4), break_by = 0.2, title = "FC Template ICA (VB2)", cor=TRUE, lines = change, labels = labels, lwd_lines=1.5, cols=cols)#, col_lines = 'white')
plot_FC(sqrt(FC_dist_mean[,,3]), zlim=c(0, 0.4), break_by = 0.2, title = "Template ICA", cor=TRUE, lines = change, labels = labels, lwd_lines=1.5, cols=cols)
dev.off()

### ICC

FC1 <- FC[,,,1,]
FC2 <- FC[,,,2,]
FC_var_tot1 <- apply(FC1, c(1,2,4), var)
FC_var_tot2 <- apply(FC2, c(1,2,4), var)
FC_var_tot <- (FC_var_tot1 + FC_var_tot2)/2
FC_var_noise <- 1/2 * apply(FC1 - FC2, c(1,2,4), var)
FC_var_sig <- FC_var_tot - FC_var_noise
FC_var_sig[FC_var_sig < 0] <- 0
FC_ICC <- FC_var_sig/FC_var_tot

pdf(file.path(main_dir,'plots','GICA25',paste0('FC_ICC.pdf')), width=5, height=5)
cols <- viridis_pal(option = 'D')(7)
plot_FC(FC_ICC[,,1], zlim=c(0, 0.6), break_by = 0.3, title = "FC Template ICA (VB1)", cor=TRUE, lines = change, labels = labels, lwd_lines=1.5, cols=cols)#, col_lines = 'white')
plot_FC(FC_ICC[,,2], zlim=c(0, 0.6), break_by = 0.3, title = "FC Template ICA (VB2)", cor=TRUE, lines = change, labels = labels, lwd_lines=1.5, cols=cols)#, col_lines = 'white')
plot_FC(FC_ICC[,,3], zlim=c(0, 0.6), break_by = 0.3, title = "Template ICA", cor=TRUE, lines = change, labels = labels, lwd_lines=1.5, cols=cols)
dev.off()

pdf(file.path(main_dir,'plots','GICA25',paste0('FC_ICC_change.pdf')), width=5, height=5)
plot_FC(FC_ICC[,,1] - FC_ICC[,,3], zlim=c(-0.3, 0.3), break_by = 0.3, title = "FC-tICA (VB1) vs. tICA", cor=TRUE, lines = change, labels = labels, lwd_lines=1.5, cols_rev = TRUE)#, col_lines = 'white')
plot_FC(FC_ICC[,,2] - FC_ICC[,,3], zlim=c(-0.3, 0.3), break_by = 0.3, title = "FC-tICA (VB2) vs. tICA", cor=TRUE, lines = change, labels = labels, lwd_lines=1.5, cols_rev = TRUE)#, col_lines = 'white')
dev.off()

#-------------------------------------------
# IC Spatial Maps and Standard Errors

S_VB2 <- readRDS(file=file.path(result_dir,'S_VB2.rds'))
S[,,,,2] <- S_VB2[,,,,2]
# saveRDS(S, file=file.path(result_dir,'S.rds'))

SE_VB2 <- readRDS(file=file.path(result_dir,'SE_VB2.rds'))
SE[,,,,2] <- SE_VB2[,,,,2]
# saveRDS(SE, file=file.path(result_dir,'SE.rds'))


### Visualize 3 subjects

idx <- c(11,19)

for(i in 1:3){

  for(j in 1:2){

      ### IC MAPS
      zlim <- c(-0.2,0.2)
      plot(newdata_xifti(xii, S[,,i,j,1]), idx = idx, zlim=zlim, fname=file.path(main_dir,'images',paste0('GICA',nQ),paste0('S_subj',i,'_sess',j,'_FCtICA1')))
      plot(newdata_xifti(xii, S[,,i,j,2]), idx = idx, zlim=zlim, fname=file.path(main_dir,'images',paste0('GICA',nQ),paste0('S_subj',i,'_sess',j,'_FCtICA2')))
      plot(newdata_xifti(xii, S[,,i,j,3]), idx = idx, zlim=zlim, fname=file.path(main_dir,'images',paste0('GICA',nQ),paste0('S_subj',i,'_sess',j,'_tICA')))
      plot(newdata_xifti(xii, S[,,i,j,4]), idx = idx, zlim=zlim, fname=file.path(main_dir,'images',paste0('GICA',nQ),paste0('S_subj',i,'_sess',j,'_DR')))

      ### IC Standard Errors
      zlim <- c(0.025,0.035)
      plot(newdata_xifti(xii, SE[,,i,j,1]), idx = idx, zlim=zlim, fname=file.path(main_dir,'images',paste0('GICA',nQ),paste0('SE_subj',i,'_sess',j,'_FCtICA1')))
      plot(newdata_xifti(xii, SE[,,i,j,2]), idx = idx, zlim=zlim, fname=file.path(main_dir,'images',paste0('GICA',nQ),paste0('SE_subj',i,'_sess',j,'_FCtICA2')))
      plot(newdata_xifti(xii, SE[,,i,j,2]), idx = idx, zlim=zlim, fname=file.path(main_dir,'images',paste0('GICA',nQ),paste0('SE_subj',i,'_sess',j,'_tICA')))
  }
}

### ICC of Spatial IC Maps

#compute & summarize
var_tot1 <- apply(S[,,,1,], c(1,2,4), var) #var across subjects for each location, IC, and method (session 1)
var_tot2 <- apply(S[,,,2,], c(1,2,4), var) #var across subjects for each location, IC, and method (session 2)
var_tot <- (var_tot1 + var_tot2)/2
var_noise <- 1/2 * apply(S[,,,1,] - S[,,,2,], c(1,2,4), var)
var_sig <- var_tot - var_noise
var_sig[var_sig < 0] <- 0
ICC <- var_sig/var_tot
saveRDS(ICC, file = file.path(result_dir, 'S_ICC.rds'))

#ICC
plot(newdata_xifti(xii, ICC[,,1]), idx=1:nQ, zlim=c(0,0.6), colors='viridis', fname=file.path(main_dir,'images','GICA25','S_ICC_FCtICA1'))
plot(newdata_xifti(xii, ICC[,,2]), idx=1:nQ, zlim=c(0,0.6), colors='viridis', fname=file.path(main_dir,'images','GICA25','S_ICC_FCtICA2'))
plot(newdata_xifti(xii, ICC[,,3]), idx=1:nQ, zlim=c(0,0.6), colors='viridis', fname=file.path(main_dir,'images','GICA25','S_ICC_tICA'))
plot(newdata_xifti(xii, ICC[,,4]), idx=1:nQ, zlim=c(0,0.6), colors='viridis', fname=file.path(main_dir,'images','GICA25','S_ICC_DR'))

#ICC change vs. tICA
ICC_df <- NULL
for(qq in 1:nQ){
  ICC_df_q <- data.frame(IC = IC_order[qq], 
                         temp_mean = temp_mean[,qq], 
                         ICC = c(ICC[,qq,]), 
                         ICC_tICA = ICC[,qq,3],
                         method = rep(c('FCtICA1','FCtICA2','tICA','DR'), each=18792))
  ICC_df <- rbind(ICC_df, ICC_df_q)
}
ICC_df$ICC_diff_tICA <- ICC_df$ICC - ICC_df$ICC_tICA 

quantile(abs(ICC_df$temp_mean), 0.99) #0.2453841

ICC_df$method <- factor(ICC_df$method, 
                         levels = c('FCtICA1','FCtICA2','tICA','DR'), 
                         labels=c('FC-tICA (VB1)', 'FC-tICA (VB2)', 'tICA', 'DR'))
ggplot(ICC_df, aes(x=abs(temp_mean), y=ICC, group=method, color=method)) + 
  geom_smooth() + xlab('IC Magnitude') + xlim(0, 0.2453841) + 
  theme_few()
ICC_df2 <- subset(ICC_df, method != 'tICA')

pdf(file.path(main_dir, 'plots', 'S_ICC_vs_tICA.pdf'), width=5, height=5)
ggplot(ICC_df2, aes(x=abs(temp_mean), y=ICC_diff_tICA, group=method, color=method)) + 
  geom_smooth() + geom_hline(yintercept=0, linetype=2) +
  xlab('IC Magnitude') + ylab('ICC change versus tICA') +
  xlim(0, 0.2453841) +
  scale_color_manual('Algorithm', values=c('blue','turquoise','grey')) +
  theme_few() + theme(legend.position='bottom')
dev.off()

### Explore residual autocorrelation

load(file.path(result_dir,paste0('result_subj5_sess1_VB1.RData')))

Y <- result_ijq$BOLD
A <- result_ijq$A
S <- as.matrix(result_ijq$subjICmean)
E <- Y - S %*% t(A)

#just motor networks
A_motor <- A[,8,drop=FALSE]
S_motor <- t(solve(crossprod(A_motor)) %*% crossprod(A_motor, t(Y)))
E_motor <- Y - S_motor %*% t(A_motor)

#randomly select 100 voxels
acf_df <- NULL
set.seed(1234)
rvox <- sample(1:18792, 100)
for(k in 1:100){
  v <- rvox[k]
  acf_v <- acf(E[v,], plot=FALSE, lag.max=10)$acf
  acf_motor_v <- acf(E_motor[v,], plot=FALSE, lag.max=10)$acf
  acf_df_v <- data.frame(acf = c(acf_v, acf_motor_v), 
                         model = rep(c('ICA','task'), each=11),
                         lag = 0:10,
                         voxel = v)
  acf_df <- rbind(acf_df, acf_df_v)
}

pdf('plots/ACFs.pdf', height=3, width=6)
ggplot(acf_df, aes(x=lag, y=acf, color=model)) +
  geom_line(alpha = 0.1, aes(group=interaction(model, voxel))) +
  stat_summary(fun.y=mean, aes(group=model), geom='line') +
  scale_x_continuous(breaks=0:10) +
  theme_few()
dev.off()


