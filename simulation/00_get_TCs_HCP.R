### This code generates dual regression-based time courses for FC template ICA simulation study
### The output file TCs.rds has been supplied in this repository so it is NOT necessary to run this 
### script unless the user choses to

library(ciftiTools)
ciftiTools.setOption('wb_path', '/Applications/workbench/') # Set path
 
main_dir <- '~/Documents/Github/FC-TemplateICA-paper/simulation/'
setwd(main_dir)

data_dir <- '/N/project/hcp_dcwan' # HCP data directory, set path

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

###########################################################################################
### Hardcoded dual_reg function from an older version of templateICAr package used      ###
### exclusively for this script.                                                        ###
###########################################################################################

dual_reg <- function(
  BOLD, GICA,
  scale=c("global", "local", "none"), scale_sm_xifti=NULL, scale_sm_FWHM=2,
  detrend_DCT=0, center_Bcols=FALSE, normA=FALSE){

  stopifnot(is.matrix(BOLD))
  stopifnot(is.matrix(GICA))
  if (is.null(scale) || isFALSE(scale)) { scale <- "none" }
  if (isTRUE(scale)) {
    warning(
      "Setting `scale='global'`. Use `'global'` or `'local'` ",
      "instead of `TRUE`, which has been deprecated."
    )
    scale <- "global"
  }
  scale <- match.arg(scale, c("global", "local", "none"))
  if (!is.null(scale_sm_xifti)) { stopifnot(ciftiTools::is.xifti(scale_sm_xifti)) }
  stopifnot(is.numeric(scale_sm_FWHM) && length(scale_sm_FWHM)==1)
  stopifnot(is.logical(normA) && length(normA)==1)

  if (any(is.na(BOLD))) { stop("`NA` values in `BOLD` not supported with DR.") }
  if (any(is.na(GICA))) { stop("`NA` values in `GICA` not supported with DR.") }

  nV <- nrow(BOLD) #number of data locations
  nT <- ncol(BOLD) #length of timeseries
  if(nV < nT) warning('More time points than voxels. Are you sure?')
  if(nV != nrow(GICA)) {
    stop('The number of voxels in dat (', nV, ') and GICA (', nrow(GICA), ') must match')
  }

  nQ <- ncol(GICA) #number of ICs
  if(nQ > nV) warning('More ICs than voxels. Are you sure?')
  if(nQ > nT) warning('More ICs than time points. Are you sure?')

  # Center each voxel timecourse. Do not center the image at each timepoint.
  # Standardize scale if `scale`, and detrend if `detrend_DCT`.
  # Transpose it: now `BOLD` is TxV.
  BOLD <- t(norm_BOLD(
    BOLD, center_rows=TRUE, center_cols=center_Bcols,
    scale=scale, scale_sm_xifti=scale_sm_xifti, scale_sm_FWHM=scale_sm_FWHM,
    detrend_DCT=detrend_DCT
  ))

  # Center each group IC across space. (Used to be a function argument.)
  center_Gcols <- TRUE
  if (center_Gcols) { GICA <- fMRItools::colCenter(GICA) }

  # Estimate A (IC timeseries).
  # We need to center `BOLD` across space because the linear model has no intercept.
  A <- ((BOLD - rowMeans(BOLD, na.rm=TRUE)) %*% GICA) %*% chol2inv(chol(crossprod(GICA)))

  # Center each subject IC timecourse across time.
  # (Redundant. Since BOLD is column-centered, A is already column-centered.)
  # A <- fMRItools::colCenter(A)

  # Normalize each subject IC timecourse if `normA`.
  if (normA) { A <- scale(A) }

  # Estimate S (IC maps).
  # Don't worry about the intercept: `BOLD` and `A` are centered across time.
  S <- solve(a=crossprod(A), b=crossprod(A, BOLD))

  #return result
  list(S = S, A = A)
}