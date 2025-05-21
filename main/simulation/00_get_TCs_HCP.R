### This code generates dual regression-based time courses for FC template ICA simulation study

library(ciftiTools)
ciftiTools.setOption('wb_path','~/') #re-define wb_path because this code is run locally or in RED

#library(devtools)
#install_github('mandymejia/templateICAr')

### Use templateICAr version 7.0 from  https://github.com/mandymejia/templateICAr/tree/7.0 just for this file and rest of the
### files uses templateICAr latest version(10.0)
library(templateICAr)

main_dir <- '~/Documents/Github/FC-TemplateICA-paper'
result_dir <- '	/N/project/FCTemplateICA' #Slate directory to store results
github_results <- file.path(main_dir, 'results')
plot_dir <- file.path(main_dir, 'plots') 
setwd(main_dir)

data_directory <- "/N/project/hcp_dcwan"
#list of subjects
subfolders <- list.dirs(data_directory, full.names = FALSE, recursive = FALSE)
subjects <- subfolders[grepl("^\\d{6}$", subfolders)]
N <- length(subjects)


#pick ICs to generate TCs
GICA_fname <- '~/Documents/Github/FC-TemplateICA-paper/data/melodic_IC_25.dscalar.nii'
GICA <- as.matrix(read_cifti(GICA_fname))
ICs <- c(1, 3, 4, 2, 16) #three visual ICs, 1 DMN, 1 motor
L <- length(ICs)

#file path to all subjects
cifti_fullnames <- file.path(data_directory, subjects, "MNINonLinear", "Results", "rfMRI_REST1_LR", "rfMRI_REST1_LR_Atlas.dtseries.nii")


Amats <- array(dim=c(1200, L, N))
for(ii in 1:N){
  
  if(!file.exists(cifti_fullnames[ii])) next()
  print(ii)
  BOLD_ii <- as.matrix(read_cifti(cifti_fullnames[ii]))
  if(ncol(BOLD_ii) != 1200) next() 
  DR_ii <- templateICAr:::dual_reg(BOLD_ii, GICA, detrend_DCT=10)
  Amats[,,ii] <- DR_ii$A[,ICs]
  
  if(ii/100 == round(ii/100)) {
    print('saving')
    output_file <- "~/Documents/Github/FC-TemplateICA-paper/results/00_get_TCs_HCP/TCs.RDS"
    
    # Create directory if it doesn't exist
    output_dir <- dirname(output_file)
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    
    # Save the object
    saveRDS(Amats, output_file)
  }
}



