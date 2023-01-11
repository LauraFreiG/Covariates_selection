#############################################################
# Code to implement the scenario B of the ANFCM paper
#############################################################

# Parameters
n=60 # sample size (n=60,100)
d=0  # d value (d=0,3,7)
m=81
B_boot=200
M=1000

#Libraries
library(MASS)
library(refund)
library(mgcv)
source("datagenALL_bis.r")
source("test_ANFCM.R")


#Progress bar
pb <- txtProgressBar(max = M, style = 3) #barra progreso
cat("\n")
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)


pvalue_vec=rep(0,M)
for(i in 1:M){
  
  # Data generation
  null.data.dn = anova.datagen(n=n, m=m, mi_low=NULL, mi_up=NULL, F=F.anova2, Q=2, d=d, 
                               error_type=3, tau=0.6^2, sig2noise=0.9^2, seed=10+i)

  
  # Test
  result = test.anova(
    y=null.data.dn$Yanova,	# n-by-m
    w=list(null.data.dn$Weval[[2]]),	# n-by-m
    test.type = 1,
    nbasis.null=7,		# number of null model basis 
    nbasis.full=c(7,7),		# number of full model basis 
    pve=0.99,		# percent of variance explained is 99% 
    B=B_boot,			# number of bootstrap replication 
    seed.B=c(1:B_boot),
    l=i
  )
  
  pvalue_vec[i]=result$pvalue
  
  setTxtProgressBar(pb, i)
}
close(pb)



cat("\n Pvalues for 1%, 5% and 10%: \n")
c(mean(pvalue_vec<=0.01),mean(pvalue_vec<=0.05),mean(pvalue_vec<=0.1))



