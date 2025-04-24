main_dir <- '~/Dropbox/RESEARCH/FCTemplateICA/simulation/'
#main_dir <- '~/Dropbox (Brown)/FCTemplateICA/simulation/'
setwd(main_dir)
source('0_setup.R')

roxygen2::roxygenize('~/Documents/Github/templateICAr/')

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
template_mean <- newdata_xifti(GICA, template$template$mean)
template_noFC <- template; template_noFC$template$FC <- NULL
template_FC_mean <- template$template$FC$psi/(template$template$FC$nu - Q - 1)
FC_mean <- cov2cor(template_FC_mean)

### ANALYZE FULL-DURATION DATA

algos <- c('tICA','VB1','VB2')
num_alg <- length(algos)
results_list <- vector('list', num_alg)
names(results_list) <- algos

### GET SURFACE MESH FOR SPATIAL ESS

xii <- read_cifti(BOLD_fnames_test[1],
                  brainstructures='left',
                  surfL_fname = ciftiTools.files()$surf[1])
FV <- as.matrix(xii$surf$cortex_left$faces)
P <- as.matrix(xii$surf$cortex_left$vertices)
inds <- which(xii$meta$cortex$medial_wall_mask$left)
mesh <- list(P, FV, inds)

result_dir <- 'results'

for(ii in 1:50){
  print(paste0('~~~~~~~~~~~~~~~~ SUBJECT ',ii,' ~~~~~~~~~~~~~~~~'))
  #load(file=paste0('results/testsubj',ii,'.RData'))

  results_ii <- results_list
  #comptime_ii <- rep(NA, num_alg)

  for(aa in 2:num_alg){ #skip tICA since the results are available within FC-tICA
    print(alg <- algos[aa])

    result_aa <- templateICA(BOLD = BOLD_fnames_test[ii],
                             template = template,
                             scale = 'none',
                             Q2 = 0,
                             #eps_inter = eps_inter,
                             brainstructures = 'left', # For tICA at least, it looks like the areas of engagement stay very similar with additional convergence (0.001 to 1e-6), but the background areas get a little less accurate, maybe due to over-fitting noise?
                             verbose=TRUE,
                             method_FC = alg,
                             reduce_dim = FALSE,
                             time_inds = 1:600) # use half the data for model estimation, the other half for ground truth FC

    #save subject i results
    results_ii[[aa]] <- result_aa
    save(results_ii, file=file.path(result_dir, paste0('testsubj',ii,'.RData')))

  } #end loop over algos
} #end loop over subjects




