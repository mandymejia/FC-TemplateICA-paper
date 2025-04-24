
ICA.gibbs = function(X, V, TT, Q, N=20, alpha, sigmaa, mua, betae, alphae, sigmae, A, theta, Z, S, mu, sigmaN, N.update, folds = 4, sigma0, mu0){

sigN = diag(sigmaN^2)
sigE = diag(rep(sigmae^2, TT))
	
## Sets of updates
sigmae.vec=rep(1, N.update)
S.matr=matrix(0, V*Q, N.update)
A.matr=matrix(0, TT*Q, N.update)

sigmae.vec[1]=sigmae
S.matr[,1] = as.vector(t(S))
A.matr[,1]=as.vector(t(A))

# Variables for updates
#n.v = rep(0, N)
n.v = c()
for(q in 1:Q){
	n.v = c(n.v, list(rep(0, N[q])))
}
Z.temp = matrix(0, Q, V)
S.temp = matrix(0, Q, V)
A.temp = matrix(0, TT,Q)

z.ind = matrix(1:V, folds, V/folds)

n=2
## The Gibbs Sampler
while(n<=N.update){

# Update the posterior value of S
mu.mat = matrix(0, Q, V)
for(q in 1:Q){
	mu.mat[q,]=mu[[q]][Z[q,]]
}
	
sigma.s = solve(t(A) %*% solve(sigE) %*% A + solve(sigN))
mu.s = sigma.s  %*%  (t(A) %*% solve(sigE) %*% X + solve(sigN) %*% mu.mat)
sq=t(chol(sigma.s))
st.n=matrix(rnorm(Q*V, 0, 1),Q,V)
S.temp=mu.s+sq%*%st.n

# Update the posterior value of A
sigA = diag(rep(sigmaa^2,Q))
muA = matrix(mua,Q, TT)
sigma.a = solve(S%*%t(S)/sigmae^2+solve(sigA))
mu.a = sigma.a%*%(S%*%t(X)/sigmae^2 + solve(sigA)%*%muA)
sq.a = t(chol(sigma.a))
st.na=matrix(rnorm(Q*TT),Q,TT)
A.temp=t(mu.a+sq.a %*% st.na)

zz = n %% folds
if(zz == 0){zz = folds}

Z.temp = Z
# Update the posterior value of Z
for(q in 1:Q){
	prob.val1 = matrix(0, dim(z.ind)[2], N[q])
	for(k in 1:N[q]){
		prob.val1[,k] = exp(log(theta[[q]][k])-(S[q,z.ind[zz,]]-mu[[q]][k])^2/(2*sigmaN[q]^2)) 
	}
	prob.sum1 = apply(prob.val1, 1, sum)
	prob.matr1 = exp(log(prob.val1)-log(prob.sum1))
	Z.temp[q,z.ind[zz,]] = apply(prob.matr1, 1, function(x){sample(1:N[q], 1, prob = x)})
}
Z=Z.temp
#Z[Z == 0] = 1

temp=transform.sa(S.temp, A.temp,rows=TRUE)
S=temp[[1]]
A=temp[[2]]

mua = rnorm(1, (mean(A.temp)*sigma0^2 + mu0*sigmaa^2/TT^2)/(sigma0^2+sigmaa^2/TT^2), sqrt((sigmaa^2*sigma0^2)/(sigma0^2*TT^2+sigmaa^2)))
sigmaa = sqrt(rinvgamma(1, alphae + TT^2/2, betae + 1/2*sum((A.temp-mua)^2)))

theta.temp = c() ##
# Update the posterior value of theta
for(q in 1:Q){
	for(k in 1:N[q]){
		n.v[[q]][k]=sum(Z[q,]==k)
		}
	theta.temp = c(theta.temp, list(rdirichlet(n=1,alpha[[q]]+n.v[[q]])))
}
theta=theta.temp

# Update the posterior value of sigmae
sigmae = sqrt(rinvgamma(1,alphae+V*TT/2,betae+1/2*sum((X-A%*%S)^2)))
sigE=diag(rep(sigmae^2,TT))

sigmae.vec[n]=sigmae
S.matr[,n]=as.vector(t(S))
A.matr[,n]=as.vector(t(A))

if(n==(N.update/2)){print(n)}
n=n+1
}

return(c(list(A.matr),list(S.matr),list(sigmae.vec)))
}



