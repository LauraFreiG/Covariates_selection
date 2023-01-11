#####################################################################
# FLCM competitor in Scenario A: linear model
#####################################################################


#Simulation parameters
M=2000 # Monte-Carlo replicates
p=2
t=0:24
n=20 #sample size (n=20,60,100)
set.seed(2020)
f1=5*sin(t*pi/12)
f2=-(t-20)^2/50-4
beta1=-(((t-15)/10)^2+0.8)
beta2=0.01*( (t-12)^2 - 12^2 +100 )
B_boot_n=1000 #Bootstrap resamples
hyp="H0" #simulation hypothesis (hyp="H0","Ha","Ha2")

#Libraries
library(Rcpp)
library(RcppArmadillo)
library(MASS)
library(fda.usc)
library(refund)
library(Matrix)
source("test.R")


#Progress bar
pb <- txtProgressBar(max = M, style = 3) #barra progreso
cat("\n")
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)



pvalue_vec1=pvalue_vec2=rep(0,M)
for(m in 1:M){

    # Data generation----
    set.seed(2020+m)
    gp_func<-function(a, f){ return(t(apply(gp[a,], 1, function(z){z+f})))  }
    gp= rproc2fdata((p+1)*n, t, sigma = "vexponential",
                    par.list = list(scale = 0.1, theta=10))$data
    x1= gp_func(1:n,f1)
    x2= gp_func((n+1):(2*n),f2)
    if(hyp=="H0"){
      y= gp[(p*n+1):((p+1)*n),]
    } else if(hyp=="Ha"){
      y=  t(apply(x2,1,function(z){z=z*beta2; return(z)})) +
          gp[(p*n+1):((p+1)*n),]
    } else if(hyp=="Ha2"){
      y=  t(apply(x1,1,function(z){z=z*beta1; return(z)})) +
          t(apply(x2,1,function(z){z=z*beta2; return(z)})) +
          gp[(p*n+1):((p+1)*n),]
    }else{
      stop("Wrong hypothesis provided")
    }
    
    
    # Partial tests----
    pvalue_vec1[m]=FLCM.test1(y,x1,nbas=7)
    pvalue_vec2[m]=FLCM.test1(y,x2,nbas=7)

    setTxtProgressBar(pb, m)
    
  }
close(pb)



cat("pvalues for X1 at levels 0.01, 0.05 and 0.1 are:\n",sep="")
cat(mean(pvalue_vec1<=0.01),", ",mean(pvalue_vec1<=0.05)," and ",mean(pvalue_vec1<=0.1),sep="")
cat("\n")
cat("pvalues for X2 at levels 0.01, 0.05 and 0.1 are:\n",sep="")
cat(mean(pvalue_vec2<=0.01),", ",mean(pvalue_vec2<=0.05)," and ",mean(pvalue_vec2<=0.1),sep="")

















