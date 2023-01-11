#############################################################
# Code to implement the scenario A of the FLCM paper
#############################################################

# Parameters
n=60 # sample size (n=60,100)
dd=0  # d value (dd=0,3,7)
p=1
m=81
T<-seq(0,1,l=m)
sigma<-.9
M=1000

#Progress bar
pb <- txtProgressBar(max = M, style = 3) #barra progreso
cat("\n")
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)


pvalues=c()
for(m_ind in 1:M){
set.seed(m_ind)
Arand<-matrix(0,n,p)
Brand<-matrix(0,n,p)
Crand<-matrix(0,n,p)
for(i in 1:n)
{for(j in 1:p)
{
  Arand[i,j]<-rnorm(1,0,1) 
}
}
for(i in 1:n)
{for(j in 1:p)
{
  Brand[i,j]<-rnorm(1,0,.85) 
}
}
for(i in 1:n)
{for(j in 1:p)
{
  Crand[i,j]<-rnorm(1,0,.7) 
}
}
X<-function(i,j,t){Arand[i,j]+(Brand[i,j]*sqrt(2)*sin(pi*t))+(Crand[i,j]*sqrt(2)*cos(pi*t))}
W<-matrix(0,n,m)
for(i in 1:n)
{for(l in 1:m)
{
  W[i,l]<-rnorm(1,0,.6)+ X(i,j,T[l])
}
}
Beta0<-function(t){1+(2*t)+(t^2)}
Beta1<-function(t){(t/8)}
c<-rnorm(n,0,sqrt(2))
d<-rnorm(n,0,.75)
E<-matrix(0,m,n)
for(l in 1:m)
{for(i in 1:n)
{
  E[l,i]<-sqrt(2)*c[i]*cos(pi*T[l])+sqrt(2)*d[i]*sin(pi*T[l])+rnorm(1,0,sigma) 
}
}
#simulatin Y under NULL
Y<-matrix(0,m,n)
for(l in 1:m)
{for(i in 1: n){
  Y[l,i]=Beta0(T[l])+X(i,1,T[l])*(dd)*Beta1(T[l])+E[l,i] 
}
}
Ymat<-t(Y)
Xmat<-W  ##observed
library(MASS)
library(mgcv)
library(refund)
library(fda)
library(Rcpp)
library(RcppArmadillo)
library(Matrix)
source("test.R")
pvalue<-FLCM.test1(Ymat,Xmat,nbas=7)
pvalues=c(pvalues,pvalue)
setTxtProgressBar(pb,m_ind)
}


cat("pvalues for levels 0.01, 0.05 and 0.1 are:\n",sep="")
c(mean(pvalues<=0.01),mean(pvalues<=0.05),mean(pvalues<=0.1))

