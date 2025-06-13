# If you're on a HPC computing system you will need to 'module load mesa' and 'module load pandoc'
# as well as the other necessary modules to run R.
# Also you will need to set your local library path to install specific packages, and their required versions
# .libPaths(c(normalizePath("~/R/libs"), .libPaths()))  # Add and prioritize a local library path for HPC systems
library(rgl)
library(reshape2) #1.4.4
library(ggplot2) #3.5.1
library(ggthemes) #5.1.0
library(scales) #alpha() #1.3.0
library(RColorBrewer) #1.1-3
library(abind) #1.4-5
library(viridis) #0.6.5
library(dplyr) #1.1.4
library(fMRItools) #0.5.0
library(templateICAr) #0.10.0
library(ciftiTools) #0.16.0
library(matrixStats) #1.3.0
library(piggyback) #0.1.5
library(xtable) #1.8.4
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
ciftiTools.setOption('wb_path', '~/workbench/bin_rh_linux64/wb_command') # Linux HPC,  Set path

source('sim_funs.R')

# Variable definitions for file paths:
##
##
n <- 500 + 50 # 500 training subjects (use to re-estimate template) + 50 test subjects
BOLD_fnames_user_setup <- paste0('/N/project/FCTemplateICA/data/1_make_templates/subj', 1:n, '.dtseries.nii') # Set data path to simulated data
BOLD_fnames_user_setup_longsim <- paste0('/N/project/FCTemplateICA/data/5_longsim/subj', 1:500, '.dtseries.nii') # Set data path to simulated data (long)


if (exists("calling_script") && calling_script == '5_longsim.R') {
  # For long simulation data
  BOLD_fnames <- BOLD_fnames_user_setup_longsim
} else {
  BOLD_fnames <- BOLD_fnames_user_setup
}


# Note:
# The code assumes the existence of certain directories (e.g., "plots", "results") for saving outputs.
# If running the scripts from scratch, be sure to create these directories manually or add checks like:
# if (!dir.exists("my_directory")) {
#   dir.create("my_directory")
# }