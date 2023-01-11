#####################################################################
# ANFCM competitor in Scenario B: nonlinear model
#####################################################################


#Simulation parameters
M=2000
p=2
t=seq(0,1,len=25)
n=20 #sample size (n=20,60,100)
set.seed(2020)
fx1=5*sin(24*t*pi/12)
fx2=-(24*t-20)^2/50-4
f1_fun<-function(t,x){exp((24*t+1)*x/20)-2}
f2_fun<-function(t,x){-1.2*log(x^2)*sin(2*pi*24*t/24)}
B_boot_n=1000 #Bootstrap resamples
hyp="H0" #simulation hypothesis (hyp="H0","Ha","Ha2")

#Libraries
library(Rcpp)
library(RcppArmadillo)
library(MASS)
library(fda.usc)
library(refund)
library(mgcv)
source("test_ANFCM.R")


#Progress bar
pb <- txtProgressBar(max = M, style = 3) #barra progreso
cat("\n")
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)



pvalue_vec1=pvalue_vec2=rep(0,M)
cont_err=0
for(m in 1:M){

  
    # Data generation----
    gen_data<-function(m){
      set.seed(2020+m)
      # independent gaussian processes
      gp_func<-function(a, f){ return(t(apply(gp[a,], 1, function(z){z+f})))  }
      gp= rproc2fdata((p+1)*n, t, sigma = "vexponential",
                      par.list = list(scale = 0.02, theta=10))$data
      x1= gp_func(1:n,fx1)
      x2= gp_func((n+1):(2*n),fx2)
      f1= t(apply(x1,1,function(z){z=f1_fun(t,z); return(z)}))
      f2= t(apply(x2,1,function(z){z=f2_fun(t,z); return(z)}) )
      if(hyp=="H0"){
        y= gp[(p*n+1):((p+1)*n),]
      } else if(hyp=="Ha"){
        y=  f2 + gp[(p*n+1):((p+1)*n),]
      } else if(hyp=="Ha2"){
        y=  f1 + f2 + gp[(p*n+1):((p+1)*n),]
      }else{
        stop("Wrong hypothesis provided")
      }
      rm(list=c("gp","f1","f2"))
      return(list(y,x1,x2))
    }
    
    # If error, simulate other data sample
    sal1=sal2=F
    dat_aux=gen_data(m)
    y=dat_aux[[1]]
    x1=dat_aux[[2]]
    x2=dat_aux[[3]]
    
    sal1=tryCatch(test.anova(
      y=y,	# n-by-m
      w=list(x1),	# n-by-m
      test.type = 1,
      nbasis.null=7,		# number of null model basis 
      nbasis.full=c(7,7),		# number of full model basis 
      pve=0.99,		# percent of variance explained is 99% 
      B=B_boot_n,			# number of bootstrap replication 
      seed.B=c(1:B_boot_n),
      l=m
    ),error=function(e){TRUE})
    
    
    sal2=tryCatch(test.anova(
      y=y,	# n-by-m
      w=list(x2),	# n-by-m
      test.type = 1,
      nbasis.null=7,		# number of null model basis 
      nbasis.full=c(7,7),		# number of full model basis 
      pve=0.99,		# percent of variance explained is 99% 
      B=B_boot_n,			# number of bootstrap replication 
      seed.B=c(1:B_boot_n),
      l=m
    ),error=function(e){TRUE})
    
    
    while((is.logical(sal1) | is.logical(sal2)) && (sal1==T | sal2==T)){
      cont_err=cont_err+1
      dat_aux=gen_data(M+cont_err)
      y=dat_aux[[1]]
      x1=dat_aux[[2]]
      x2=dat_aux[[3]]
      
      sal1=tryCatch(test.anova(
        y=y,	# n-by-m
        w=list(x1),	# n-by-m
        test.type = 1,
        nbasis.null=7,		# number of null model basis 
        nbasis.full=c(7,7),		# number of full model basis 
        pve=0.99,		# percent of variance explained is 99% 
        B=B_boot_n,			# number of bootstrap replication 
        seed.B=c(1:B_boot_n),
        l=m
      ),error=function(e){TRUE})
      
      
      sal2=tryCatch(test.anova(
        y=y,	# n-by-m
        w=list(x2),	# n-by-m
        test.type = 1,
        nbasis.null=7,		# number of null model basis 
        nbasis.full=c(7,7),		# number of full model basis 
        pve=0.99,		# percent of variance explained is 99% 
        B=B_boot_n,			# number of bootstrap replication 
        seed.B=c(1:B_boot_n),
        l=m
      ),error=function(e){TRUE})
    }
    
    # Partial tests----
    pvalue_vec1[m]=sal1$pvalue
    pvalue_vec2[m]=sal2$pvalue
    
    setTxtProgressBar(pb, m)
  
}
close(pb)



cat("pvalues for X1 at levels 0.01, 0.05 and 0.1 are:\n",sep="")
cat(mean(pvalue_vec1<=0.01),", ",mean(pvalue_vec1<=0.05)," and ",mean(pvalue_vec1<=0.1),sep="")
cat("\n")
cat("pvalues for X2 at levels 0.01, 0.05 and 0.1 are:\n",sep="")
cat(mean(pvalue_vec2<=0.01),", ",mean(pvalue_vec2<=0.05)," and ",mean(pvalue_vec2<=0.1),sep="")


















