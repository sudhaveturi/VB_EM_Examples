rm(list=ls())
# Code to run a simple example using Variational Bayes, EM and MCMC algorithms 

### VB for Univariate Gaussian

N = 1000
set.seed(100)
x <- rnorm(N)

mu00 = rnorm(1)
b00 = a00 = runif(1)
lambda00  = runif(1)

mu0 = rnorm(1)
b0 = runif(1)
a0 = runif(1)
lambda0  = runif(1)


### grid search

nGrid<-300

mu<-seq(from=-1,to=1,length=nGrid)
tau<-seq(from=0.1,to=2,length=nGrid)

logp<- logq1<- logq2 <- logq3 <- logq4 <- matrix(nrow=length(mu),ncol=length(tau),NA)


reps=1000

lambda_N =  b_N = E_tau = E_mu = E_mu2 = q_mu = q_tau = KL = matrix(,reps,1)
E_tau[1,] = a0/b0

mu_N = ((lambda0*mu0)+(N*mean(x)))/(lambda0+N)
a_N = a0+((N+1)/2)

lambda_N[1,] = (lambda0+N)*E_tau[1,]
q_mu[1,] = rnorm(1,mu_N,lambda_N[1,])
E_mu[1,] = mu_N
E_mu2[1,] = (1/lambda_N[1,])+mu_N^2
b_N[1,] = b0+0.5*((E_mu2[1,])*(N+lambda0)-(2*E_mu[1,]*((N*mean(x))+(lambda0*mu0)))+((lambda0*(mu0^2))+sum(x^2)))

q_tau[1,] = rgamma(1,a_N,b_N[1,])


for(i in 2:reps){
  
  lambda_N[i,] = (lambda0+N)*E_tau[i-1,]
  q_mu[i,] = rnorm(1,mu_N,lambda_N[i,])
  E_mu[i,] = mu_N
  E_mu2[i,] = (1/lambda_N[i,])+mu_N^2
  b_N[i,] = b0+0.5*((E_mu2[i,])*(N+lambda0)-(2*E_mu[i,]*((N*mean(x))+(lambda0*mu0)))+((lambda0*(mu0^2))+sum(x^2)))
  q_tau[i,] = rgamma(1,a_N,b_N[i,])
  E_tau[i,] = a_N/b_N[i,]
  
}

for(i in 1:nrow(logp)){
  
  for(j in 1:ncol(logp)){
    
    logp[i,j]<-sum(log(dnorm(x,mu[i],(1/sqrt(tau[j])))))+log(dnorm(mu[i],mu00,(1/sqrt(lambda00*tau[j]))))+log(dgamma(tau[j],a00,b00))
    logq1[i,j]<- log(dnorm(mu[i],mu_N,(1/sqrt(lambda_N[1,]))))+log(dgamma(tau[j],a_N,b_N[1,])) 
    logq2[i,j]<- log(dnorm(mu[i],mu_N,(1/sqrt(lambda_N[2,]))))+log(dgamma(tau[j],a_N,b_N[1,])) 
    logq3[i,j]<- log(dnorm(mu[i],mu_N,(1/sqrt(lambda_N[2,]))))+log(dgamma(tau[j],a_N,b_N[2,])) 
    logq4[i,j]<- log(dnorm(mu[i],mu_N,(1/sqrt(lambda_N[100,]))))+log(dgamma(tau[j],a_N,b_N[100,])) 
    
    #  logq1[i,j]<--(E_tau[1,]/2)*((lambda0*(mu[i]-mu0)^2)+sum((x-mu[i])^2))+ ((((N-1)/2)+a0)*log(tau[j]))-
    #     (tau[j]/2)*(b0+0.5*((E_mu2[1,])*(N+lambda0)-(2*E_mu[1,]*((N*mean(x))+(lambda0*mu0)))+((lambda0*(mu0^2))+sum(x^2))))
    print(paste(i,j,sep='-'))
    
  }
  
}

contour(x=mu,y=tau,z=logp,levels=c(seq(from=-5000,to=5,by=20)),xlab='mu',ylab='tau')
contour(x=mu,y=tau,z=logq4,levels=c(seq(from=-5000,to=5,by=20)),xlab='mu',ylab='tau',col="red",add=T)
###############################################################

# SAT-V scores
y = c(28,8,-3,7,-1,1,18,12)
sigma = c(15,10,16,11,9,11,10,18)

# VB
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


### Gibbs sampler

J=8
alpha.update <- function (){
  alpha.hat <- (mu/tau^2 + y/sigma^2)/(1/tau^2 + 1/sigma^2)
  V.alpha <- 1/(1/tau^2 + 1/sigma^2)
  rnorm (J, alpha.hat, sqrt(V.alpha))
}
mu.update <- function (){
  rnorm (1, mean(alpha), tau/sqrt(J))
}
tau.update <- function (){
  sqrt(sum((alpha-mu)^2)/rchisq(1,J-1))
}

n.chains <- 5
n.iter <- 1000
sims <- array (NA, c(n.iter, n.chains, J+2))
dimnames (sims) <- list (NULL, NULL,
                         c (paste ("alpha[", 1:8, "]", sep=""), "mu", "tau"))
for (m in 1:n.chains){
  mu <- rnorm (1, mean(y), sd(y))
  tau <- runif (1, 0, sd(y))
  for (t in 1:n.iter){
    alpha <- alpha.update ()
    mu <- mu.update ()
    tau <- tau.update ()
    sims[t,m,] <- c (alpha, mu, tau)
  }
}


## EM algorithm

# initial values

mu_old = rnorm(1)
tau_old = runif(1)
reps = 100
mu = tau = matrix(,reps,1)
mu[1,] = mu_old
tau[1,] = tau_old
alpha_hat = V_alpha_hat = matrix(,reps,8)

for(i in 2:reps){
  alpha_hat[i-1,] = ((y/sigma^2)+(mu[i-1,]/tau[i-1,]^2))/((1/sigma^2)+(1/tau[i-1,]^2))
  V_alpha_hat[i-1,] = 1/((1/sigma^2)+(1/tau[i-1,]^2))
  mu[i,] = (1/8)* sum(alpha_hat[i-1,])   
  tau[i,] = sqrt((1/7)*sum(((alpha_hat[i-1,]-mu[i,])^2)+V_alpha_hat[i-1,]))
}
i=101
alpha_hat[i-1,] = ((y/sigma^2)+(mu[i-1,]/tau[i-1,]^2))/((1/sigma^2)+(1/tau[i-1,]^2))
V_alpha_hat[i-1,] = 1/((1/sigma^2)+(1/tau[i-1,]^2))

alpha_EM = alpha_hat
mu_EM = mu
tau_EM = tau


#############################
####Summary

# predicts + interval
plotfn <- function(parameter,mean,S2){
  df = data.frame(x=c(1:100),y=parameter)
  newx <- seq(min(df$x), max(df$x), length.out=100)
  preds50 <- data.frame(parameter,mean-(0.67*sqrt(S2)),mean+(0.67*sqrt(S2)))
  preds90 <- data.frame(parameter,mean-(1.645*sqrt(S2)),mean+(1.645*sqrt(S2)))
  
  colnames(preds) = c("fit","lower","upper")
  # plot
  plot(preds[,1], type= "n", ylim=range(preds90),xlab="Iterations",ylab="",xlim=c(0,105))
  # add fill
  polygon(c(rev(newx), newx), c(rev(preds90[ ,3]), preds90[ ,2]), col = 'grey80', border = NA)
  polygon(c(rev(newx), newx), c(rev(preds50[ ,3]), preds50[ ,2]), col = 'grey20', border = NA)
  
  # model
  lines(mean)
  # intervals
  lines(newx, preds90[ ,3], lty = 'dashed', col = 'red')
  lines(newx, preds90[ ,2], lty = 'dashed', col = 'red')
  lines(newx, preds50[ ,3], lty = 'dashed', col = 'red')
  lines(newx, preds50[ ,2], lty = 'dashed', col = 'red')
  
}

plotfn(alpha_VB[-1,1],M_alpha[-1,1],S2_alpha[-1,1])
points(rep(102,3),quantile(sims[,,1])[2:4],pch=c(1,19,1))
points(104,alpha_EM[100,1],pch=19,col="red")
title(bquote("variational inference for"~alpha[1]))

plotfn(mu_VB[-1,],M_mu[-1,],S2_mu[-1,])
points(rep(102,3),quantile(sims[,,9])[2:4],pch=c(1,19,1))
points(104,mu_EM[100,1],pch=19,col="red")
title(bquote("variational inference for"~mu))

plotfn(tau_VB[-1,],sqrt((7/5)*M2_tau[-1,]),sqrt(((2*(7^2))/(((7-2)^2)*(7-4)))*M2_tau[-1,]^2))
points(rep(102,3),quantile(sims[,,10])[2:4],pch=c(1,19,1))
points(104,tau_EM[100,1],pch=19,col="red")
title(bquote("variational inference for"~tau))

###############################################################



