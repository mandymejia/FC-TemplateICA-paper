# If you're on a HPC computing system you will need to 'module load mesa' and 'module load pandoc'
# and will need to set your local library path to install specific packages, and their required versions
.libPaths(c(normalizePath("~/R/libs"), .libPaths()))  # Add and prioritize a local library path
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
# library(matrixStats) #colVars
# library(manipulateWidget) # May need for HPC if pdfs and pngs do not render.
# library(rmarkdown) # May need for HPC if pdfs and pngs do not render.

# Note:
# Some packages require specific versions and this how you can install them:
# install.packages("remotes") 
# install.packages("devtools")
# remotes::install_version("templateICAr", version = "0.10.0", lib = "~/R/libs", repos = "https://cloud.r-project.org")
# devtools::install_github('mandymejia/templateICAr', ref='10.0') # For development versions

# Set your ciftiTools wb_path to the location of your wb_command
# ciftiTools.setOption('wb_path', '/Applications') 
# ciftiTools.setOption('wb_path', '/Applications/workbench/bin_macosxub/') # MacOS
ciftiTools.setOption('wb_path', '~/workbench/bin_rh_linux64/wb_command') # Linux HPC

source('sim_funs.R')

# Variable definitions for file paths:
##
##
if (calling_script == '5_longsim.R') {
  # For long simulation data
  BOLD_fnames <- paste0('/N/project/FCTemplateICA/data/5_longsim/subj', 1:500, '.dtseries.nii')
} else {
  BOLD_fnames <- paste0('/N/project/FCTemplateICA/data/1_make_templates/subj', 1:n, '.dtseries.nii')
}
# Note:
# The code creates a lot of directories for organizational purposes, so if you will be running the scripts
# from scratch, makes sure to create the directors or insert lines such as these:
# if (!dir.exists("my_directory")) {
#   dir.create("my_directory")
# }