
library(fastICA)
library(MASS)
library(steadyICA)

rdirichlet = function(n, alpha){
	l <- length(alpha)
	x <- matrix(rgamma(l*n, alpha),ncol = l, byrow=TRUE)
	sm <- x %*% rep(1, l)
	return(x/as.vector(sm))
}

rinvgamma = function(n, shape, scale = 1){
	return(1/rgamma(n = n, shape = shape, rate = scale))
}

source("/Users/ani/Dropbox/Bayesian\ ICA/R\ Code/StVal/binning.R")
source("/Users/ani/Dropbox/Bayesian\ ICA/R\ Code/StVal/weights.R")
source("/Users/ani/Dropbox/Bayesian\ ICA/R\ Code/transformation.R")
source("/Users/ani/Dropbox/Bayesian\ ICA/R\ Code/amari.R")
source("/Users/ani/Dropbox/Bayesian\ ICA/R\ Code/bicaRSGS.R")

########### Generate the data ##############

TT = 5
Q = 2
V = 10000
N.update = 2000
N.burn = 1000

	N = c(15, 20)
	n.b=500
	sigmae.true=0.05

error.am = c()
error.a = c()
sige = c()
run.t = c()

N.sim = 1
for(n.sim in 1:N.sim){

##### Generate S.true

	shape0=4
	scale0=0.25
	mu0=shape0*scale0
	sigma0=sqrt(mu0*scale0)
	lambw3=3
	kw3=1
	muw3=kw3*gamma(1+1/lambw3)
	varw3=(kw3^2*gamma(1+2/lambw3)-muw3^2)


	u1 = runif(V,0,1)
	u2 = runif(V,0,1)
	u3 = runif(V,0,1)
	u4 = runif(V,0,1)	
	S1 = (rweibull(V,shape=lambw3,scale=kw3)-muw3)/sqrt(varw3)
 	S2 = rnorm(V,-0.8,0.3)*(u1<0.6)+rnorm(V,1.2,0.3)*(u1>0.6)	
 	S3 = (rgamma(V,shape=shape0,scale=scale0)-mu0)/sigma0
 	
	S.true = rbind(S1, S2)
##### Generate the error

    E=matrix(rnorm(TT*V,0,sigmae.true),TT,V)

##### Define A.true
	A.true = matrix(c(0.75, -0.25, -0.5, 0.5, 0.3, 0.25, 0.52, 1, 0.12, 0.35), TT, Q)
	temp.true = transform.sa(S.true, A.true, rows = TRUE)
	A.true = temp.true[[2]]
	S.true = temp.true[[1]]
		
##### Transform S.true and A.true	

	X=A.true%*%S.true+E

## Parameters of the priors
	alphae = 100
	betae = 0.01
	mua = 0
	sigmaa = 0.5
	alpha = c()
	for(q in 1:Q){
		alpha = c(alpha, list(rep(1,N[q])))
	}

## Initial values for everything

	sigmae = 0.05
	fa = fastICA(t(X), n.comp = Q)
	A = t(fa$A)
	S = t(fa$S)
	temp1 = transform.sa(S, A, rows=TRUE)
	A = temp1[[2]]
	S = temp1[[1]]
	
## Identify good initial values for theta
	counts = matrix(0, Q, n.b)
	Tj = matrix(0, Q, n.b)	

	fout = fbins(S, Q, n.b, rows = TRUE)
	bins = fout[[1]]
	Tj = fout[[2]]	
	for(q in 1:Q){
	counts[q,] = table(cut(S[q,], breaks=bins[q,],include.lowest=TRUE))}
	counts[counts == 0] = 1
	Co = counts/V

	mu = c()
	sigmaN = c()
	for(q in 1:Q){

		mu = c(mu, list(seq(min(S[q,]), max(S[q,]), l = N[q])))
		sigmaN = c(sigmaN, 2/3*(mu[[q]][N[q]]-mu[[q]][1])/(N[q]-1))
	}

	theta = c()
	for(q in 1:Q){
	
		theta = c(theta, list(weight.seq(Tj[q,], Co[q,], mu[[q]], sigmaN[q], 0)))
	}

### Initial values for Z
	
	Z=matrix(0,Q,V)
	for(q in 1:Q){
		Z[q,] = sample(N[q], V, replace=TRUE, prob=theta[[q]])
	}
	
folds = 8  # , folds = folds

st.time=proc.time()

res = ICA.gibbs(X = X, TT = TT, Q = Q, V = V, N = N, alpha = alpha, sigmaa = sigmaa, mua = mua, betae = betae, alphae = alphae, sigmae = sigmae, A = A, theta = theta, Z = Z, S = S, mu = mu, sigmaN = sigmaN, N.update = N.update, folds = folds)

end.time=proc.time()-st.time

A.matr = res[[1]]
S.matr = res[[2]]
sigmae.vec = res[[3]]

A.m = apply(A.matr[,N.burn:N.update], 1, mean)
A.m = matrix(A.m, TT, Q, byrow=TRUE)

S.m = S.matr[, N.burn:N.update]
S.m = apply(S.m, 1, mean)
S.m = matrix(S.m, V, Q)

temp = transform.sa(t(S.m), A.m, rows=TRUE)
A.m = temp[[2]]
S.m = temp[[1]]

error.am = c(error.am, amari(solve(A.m), solve(A.true)))
error.a = c(error.a, amari(solve(A), solve(A.true)))
sige = c(sige, mean(sigmae.vec[N.burn:N.update]))
run.t = c(run.t, end.time[3]/60)

}

