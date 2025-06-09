.libPaths(c(normalizePath("~/R/libs"), .libPaths())) # for HPC, must install library(rgl) and module load mesa <add more robust comments if needed>
library(rgl)
# install.packages("remotes")
# remotes::install_version("templateICAr", version = "0.10.0") # To install specific versions.

library(templateICAr) #0.10.0
library(ciftiTools) #0.16.0
#### REINSTALL THIS ONE #library(fMRItools) #0.5.0

#### ANI CHANGE next two lines
# ciftiTools.setOption('wb_path', '/Applications')
#ciftiTools.setOption('wb_path', '/Applications/workbench/bin_macosxub/')
ciftiTools.setOption('wb_path', '~/workbench/bin_rh_linux64/wb_command')

#library(matrixStats) #colVars
library(reshape2) #1.4.4
library(ggplot2) #3.5.1
library(ggthemes) #5.1.0
library(scales) #alpha() #1.3.0
library(RColorBrewer) #1.1-3
library(abind) #1.4-5
library(viridis) #0.6.5
library(dplyr) #1.1.4
library(fMRItools)


source('sim_funs.R')

# library(manipulateWidget) # for HPC
# library(rmarkdown) # for HPC
# install.packages("rmarkdown") # or module load pandoc (HPC)
# install.packages("manipulateWidget")

