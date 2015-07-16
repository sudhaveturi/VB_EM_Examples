######### MCMC and EM algorithm for BRR
library(BGLR)
library(MASS)
library(corpcor)
data(wheat)
y<-wheat.Y[,1]

X<- wheat.X
p=ncol(X)
n=nrow(X)
ETA<-list(X=list(X=X,model="BRR"))

#MCMC
fm<-BGLR(y=y,ETA=ETA,nIter=12000,burnIn=2000)

# EM
# initial values

sigma2_old = runif(1)
sigma2_b_old = runif(1)
reps = 100
sigma2 = sigma2_b = matrix(,reps,1)
sigma2[1,] = sigma2_old
sigma2_b[1,] = sigma2_b_old
beta_hat = matrix(,reps,p)

for(i in 2:reps){
  lambda = sigma2[i-1,]/sigma2_b[i-1,]
  D = diag(p)
  D[1,1] <- 0
  beta_hat[i-1,] = chol2inv(chol(make.positive.definite(crossprod(X)+lambda*D)))%*%crossprod(X,y)
  sigma2_b[i,] = (beta_hat[i-1,]%*%D%*%beta_hat[i-1,])/p  
  sigma2[i,] = crossprod((y-(X%*%beta_hat[i-1,])))/n
  cat(i,"\n")
}

beta_EM = beta_hat[99,]
sigma2_EM = sigma2
sigma2_b_EM = sigma2_b

beta_MCMC = fm$ETA[[1]]$b
sigma2_MCMC = fm$varE
sigma2_b_MCMC = fm$ETA[[1]]$varB

plot(beta_EM)
points(beta_MCMC, col="red")

######### VB algorithm for the GBLUP

X<-scale(wheat.X,center=TRUE,scale=TRUE)
G<-tcrossprod(X)/ncol(X)
library(geoR)
set.seed(100)
M_alpha0 = S2_alpha0 = matrix(,8,1)
for(i in 1:8){
  M_alpha0[i,] = rnorm(1)
  S2_alpha0[i,] = runif(1)
}
M_mu0 = rnorm(1)
S2_mu0 = (runif(1))^2
M2_tau0 = 0
for(i in 1:8)
  M2_tau0 = M2_tau0+((M_alpha0[i,]-M_mu0)^2+S2_alpha0[i,]+S2_mu0)
M2_tau0 = M2_tau0/7

reps=100
alpha = matrix(,reps+1,8)
mu = tau2 =  KL = matrix(,reps+1,1)

M_alpha = S2_alpha = matrix(,reps+1,8)
M_mu = S2_mu = M2_tau = matrix(,reps+1,1)

M_alpha[1,] = as.vector(M_alpha0)
M_mu[1,] = M_mu0
S2_alpha[1,] = as.vector(S2_alpha0)
S2_mu[1,] = S2_mu0
M2_tau[1,] = M2_tau0

for(j in 1:reps){
  for(i in 1:8)
    alpha[j+1,i] = rnorm(1,M_alpha[j,i],sqrt(S2_alpha[j,i]))
  mu[j+1,] = rnorm(1,M_mu[j,],sqrt(S2_mu[j,]))
  tau2[j+1,] = rinvchisq(1,7,M2_tau[j,])
  
  KL[j,] = sum((((y-M_alpha[j,])^2)+S2_alpha[j,])/(2*sigma^2))+(8*log(sqrt(M2_tau[j,])))+
    sum((((M_alpha[j,]-M_mu[j,])^2)+S2_alpha[j,]+S2_mu[j,])/(2*M2_tau[j,]))
  -sum(log(sqrt(S2_alpha[j,])))-log(sqrt(S2_mu[j,]))-(10*log(sqrt(M2_tau[j,])))
  
  M_alpha[j+1,] = (((y/sigma^2)+((1/M2_tau[j,])*M_mu[j,]))/((1/sigma^2)+(1/M2_tau[j,])))
  S2_alpha[j+1,] = 1/((1/sigma^2)+(1/M2_tau[j,]))
  M_mu[j+1,] = mean(M_alpha[j+1,])
  S2_mu[j+1,] = (1/8)*(1/(1/M2_tau[j,]))
  tmp = 0
  for(k in 1:8)
    tmp = tmp+((M_alpha[j+1,k]-M_mu[j+1,])^2+S2_alpha[j+1,k]+S2_mu[j+1,]) 
  M2_tau[j+1,] = tmp/7
}

alpha_VB = alpha
mu_VB = mu
tau_VB = sqrt(tau2)











