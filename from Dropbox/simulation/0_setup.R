#devtools::install_github('mandymejia/templateICAr', ref='8.0')
library(templateICAr) #0.8.0
library(ciftiTools) #0.16.0
#### REINSTALL THIS ONE #library(fMRItools) #0.5.0

#### ANI CHANGE next two lines
ciftiTools.setOption('wb_path', '/Applications')
#ciftiTools.setOption('wb_path', '/Applications/workbench/bin_macosxub/')

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