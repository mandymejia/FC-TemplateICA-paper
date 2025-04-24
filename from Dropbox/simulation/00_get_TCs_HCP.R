### This code generates dual regression-based time courses for FC template ICA simulation study

library(ciftiTools)
ciftiTools.setOption('wb_path', '/Applications/workbench/')
#library(devtools)
#install_github('mandymejia/templateICAr')
library(templateICAr)

main_dir <- '~/Dropbox/RESEARCH/FCTemplateICA/'
setwd(main_dir)

data_dir <- '/Volumes/Lab_Data_Drive/data/HCP_Resting_State'

cifti_fname <- 'rfMRI_REST1_LR/rfMRI_REST1_LR_Atlas.dtseries.nii'

#get list of subjects
numfiles <- table(substr(list.files(data_dir), 1, 6))
subjects <- names(numfiles)
N <- length(subjects)

#pick ICs to generate TCs
GICA_fname <- 'melodic_IC_25.dscalar.nii'
GICA <- as.matrix(read_cifti(GICA_fname))
ICs <- c(1, 3, 4, 2, 16) #three visual ICs, 1 DMN, 1 motor
L <- length(ICs)

#file path to all subjects
cifti_fullnames <- file.path(data_dir, subjects, cifti_fname)

Amats <- array(dim=c(1200, L, N))
for(ii in 1:N){
  
  if(!file.exists(cifti_fullnames[ii])) next()
  print(ii)
  BOLD_ii <- as.matrix(read_cifti(cifti_fullnames[ii]))
  if(ncol(BOLD_ii) != 1200) next() 
  DR_ii <- dual_reg(BOLD_ii, GICA, detrend_DCT=10)
  Amats[,,ii] <- DR_ii$A[,ICs]
  
  if(ii/100 == round(ii/100)) {
    print('saving')
    saveRDS(Amats, 'simulation/TCs.RDS')
  }
}


