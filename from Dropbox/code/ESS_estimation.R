library(Matrix)
library(INLA)
library(fmesher)
lik <- function(theta, Y, G,C,ind = NULL) {
  kappa <- exp(theta[1])
  tau <- exp(theta[2])
  Q <- tau^2*(kappa^2*C + G)

  if(!is.null(ind)) {
    ind2 <- setdiff(1:dim(Q)[1], ind)
    Q <- Q[ind,ind] - Q[ind,ind2]%*%solve(Q[ind2,ind2],Q[ind2,ind])
  }
  R <- chol(Q)
  ld <- c(determinant(R, logarithm = TRUE, sqrt = TRUE)$modulus)
  lik <- as.double(ld - 0.5*t(Y)%*%Q%*%Y)

  return(-lik)
}

#' Estimation of effective sample size
#'
#' @param mesh fmesher mesh
#' @param Y data
#' @param ind index of the data locations in the mesh
#' @param use_inla should INLA be used to compute standard deviations? If FALSE, the
#' standard deviations is approximated as being constant over the locations
#' @param trace If true, an alternative formula based on the trace is used
#' @details
#' The functions computes the effective sample size as 1^T Q.corr 1 where Q.corr is
#' the inverse of the correlation matrix of a fitted SPDE model with alpha = 1.
#' If `trace=TRUE`, the alternative formula trace(Q^{-1})/trace(Q*Q) is used.
#'
#' @return Estimate of the effective sample size
estimate.ESS <- function(mesh, Y, ind = NULL, use_inla = TRUE, trace = FALSE) {
  fem <- fmesher::fm_fem(mesh)
  res <- optim(c(log(1),log(1)), lik, Y = Y, C = fem$c0, G = fem$g1, ind = ind)
  Q <- exp(res$par[2])^2*(exp(res$par[1])^2*fem$c0 + fem$g1)
  if(!is.null(ind)) {
    ind2 <- setdiff(1:dim(Q)[1], ind)
    Q <- Q[ind,ind] - Q[ind,ind2]%*%solve(Q[ind2,ind2],Q[ind2,ind])
  }

  if(trace == FALSE){
    if(use_inla) {
      std <- sqrt(diag(inla.qinv(Q)))

    } else {
      if(mesh$manifold == "R1") {
        if(is.null(ind)) {
          loc <- mesh$interval[1] + diff(mesh$interval)/2
          A <- fm_basis(mesh, loc)
        } else {
          A <- rep(0,dim(Q)[1])
          A[1] <- 1
        }

      } else {
        if(is.null(ind)) {
          loc <- colMeans(mesh$loc)[1:2]
          A <- fm_basis(mesh, matrix(loc,1,2))
        } else {
          A <- matrix(0,1,dim(Q)[1])
          A[1] <- 1
        }
      }
      std <- rep(sqrt(as.double(A%*%solve(Q, t(A)))), dim(Q)[1])
    }
    return(as.double(t(std)%*%Q%*%std))
  } else {
    Sigma <- solve(Q)
    return(sum(diag(Sigma))^2/sum(Sigma^2))
  }
}

# Example 1d
n <- 1001
s <- seq(from = 0, to = 1, length.out = n)

kappa <- 10
tau <- 1

mesh1 <- fmesher::fm_mesh_1d(loc=s)
fem <- fmesher::fm_fem(mesh1)

Q <- kappa^2*fem$c0 + fem$g1
R <- chol(Q)
Y <- solve(R,rnorm(n))

ESS <- estimate.ESS(mesh1,Y)
ESS2 <- estimate.ESS(mesh1,Y, use_inla = FALSE)
ESS3 <- estimate.ESS(mesh1,Y, use_inla = FALSE, trace = TRUE)
cat("estimated ESS:", ESS, ESS2, ESS3)

#example 2d

#For simplicitly, we require the mesh to contain only the observation locations
n <- 200
loc_2d_mesh <- matrix(runif(n * 2), n, 2)
mesh2 <- fm_mesh_2d(
  loc = loc_2d_mesh,
  cutoff = 0,
  max.edge = c(0.1, 0.5)
)

plot(mesh2)
points(mesh2$loc[mesh2$idx$loc,1:2],col=2)
fem <- fm_fem(mesh2)

Q <- kappa^2*fem$c0 + fem$g1
R <- chol(Q)
u <- solve(R,rnorm(mesh2$n))
ind <- mesh2$idx$loc
Y <- u[ind]

ESS <- estimate.ESS(mesh2,Y,ind)
ESS2 <- estimate.ESS(mesh2,Y,ind, use_inla = FALSE)
cat("estimated ESS:", ESS, ESS2)

# For a general cortical surface, we do:

# step 1: get the mesh

# mesh = ...

# Step 2: get the indices in the mesh where the data is

# ind = ...

# Step 3: Run estimate.ESS:

# ESS <- estimate.ESS(mesh,Y,ind)
