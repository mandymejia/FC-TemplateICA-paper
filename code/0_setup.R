# If you're on a HPC computing system you will need to 'module load mesa' and 'module load pandoc'
# and will need to set your local library path to install specific packages, and their required versions
# .libPaths(c(normalizePath("~/R/libs"), .libPaths())) #  Add and prioritize a local library path for HPC systems
library(templateICAr) #0.10.0
library(fMRItools) #0.5.0
library(ciftiTools) #0.16.0
library(dplyr) #1.1.4
library(viridis) #0.6.5
library(RColorBrewer) #1.1-3
library(matrixStats) #1.3.0
library(ggplot2) #3.5.1
library(ggthemes) #5.1.0
library(piggyback) #0.1.5
# library(rgl) #1.3.18 May need for HPC if pdfs and pngs do not render.
# library(manipulateWidget) #0.11.1 May need for HPC if pdfs and pngs do not render.
# library(rmarkdown) #2.29 May need for HPC if pdfs and pngs do not render.


# Note:
# Some packages require specific versions and this how you can install them:
# install.packages("remotes") 
# install.packages("devtools")
# remotes::install_version("templateICAr", version = "0.10.0", lib = "~/R/libs", repos = "https://cloud.r-project.org")
# devtools::install_github('mandymejia/templateICAr', ref='10.0') # For development versions

# Set your ciftiTools wb_path to the location of your wb_command
# ciftiTools.setOption('wb_path', '/Applications') # MacOS, Set path
ciftiTools.setOption('wb_path', '~/workbench/bin_rh_linux64/wb_command') # Linux HPC, Set path

source('code/functions.R') 

# Variable definitions for file paths:
##
##
session_names <- c('rfMRI_REST1_LR','rfMRI_REST2_LR')
data_dir <- '/N/project/hcp_dcwan' # HCP data directory
cifti_fnames <- file.path('MNINonLinear/Results', session_names, paste0(session_names, '_Atlas_hp2000_clean.dtseries.nii'))

# Note:
# The code assumes the existence of certain directories (e.g., "plots", "results") for saving outputs.
# If running the scripts from scratch, be sure to create these directories manually or add checks like:
# if (!dir.exists("my_directory")) {
#   dir.create("my_directory")
# }