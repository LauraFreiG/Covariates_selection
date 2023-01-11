##########################################################################
# Code to implement the scenario B of the ANFCM paper using the MDD test
#########################################################################

# Parameters
n=60 # sample size (n=60,100)
d=0  # d value (d=0,3,7)
m=81
B_boot_n=200
M=1000

#Libraries
library(MASS)
library(refund)
library(mgcv)
library(Rcpp)
library(RcppArmadillo)
library(MASS)
library(fda.usc)
library(refund)
library(EDMeasure)
sourceCpp("signif_test_MDD_FCM.cpp")
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
  
  # Partial test----
  Xmat=null.data.dn$Weval[[2]]
  Ymat=null.data.dn$Yanova
  p=1
  Tn_vec=c()
  Tn_boot_mat=c()
  for(t_ind in 1:m){
    Xn=as.matrix(Xmat[,t_ind])
    Y=Ymat[,t_ind]
    
    # T_n calculation----
    MDD=apply(Xn,2,function(z){mdd(z,Y,compute="C",center="U")})
    cn= ((n-3)^4/(n-1)^4) + (2*(n-3)^4/((n-1)^4*(n-2)^3)) + (2*(n-3)/((n-1)^4*(n-2)^3))
    
    
    A=list()
    for(j in 1:p){
      X_aux=as.vector(Xn[,j])
      XX_aux=rep(0,n^2)
      A_aux=UCenter_X(n,1,X_aux,XX_aux)
      A[[j]]=matrix(A_aux,n,n) 
    }
    Y_aux=as.vector(Y)
    YY_aux=rep(0,n^2)
    B=UCenter_Y(n,1,Y_aux,YY_aux)
    B=matrix(B,n,n)
    
    suma=0
    for(k in 1:(n-1)){
      for(l in (k+1):n){
        for(j1 in 1:p){
          A1=A[[j1]]
          for(j2 in 1:p){
            A2=A[[j2]]
            suma = suma + (A1[k,l]*A2[k,l]*(B[k,l]^2))
          }
        }
      }
    }
    S= (2/(n*(n-1)*cn))*suma
    Tn=sqrt(choose(n,2))*sum(MDD)/sqrt(S)
    Tn_vec=c(Tn_vec,Tn)
    
    
    
    # Bootstrap---
    Tn_boot=c()
    for(b in 1:B_boot_n){
      set.seed(b)
      e=rnorm(n)
      E1=matrix(e,n,n)
      E2=t(E1)
      B_boot=B*E1*E2
      MDD_boot=rep(0,p)
      for(j in 1:p){
        Aj=A[[j]]
        suma_boot=0
        suma_boot= sum(Aj*B_boot) - sum(diag(Aj*B_boot))
        MDD_boot[j]= suma_boot/(n*(n-1))
      }
      suma=suma_S_boot(e, n, p, A, B)
      S_boot=suma/choose(n,2)
      Tn_boot=c(Tn_boot, sqrt(choose(n,2))*sum(MDD_boot)/sqrt(S_boot) )
    }
    
    Tn_boot_mat=cbind(Tn_boot_mat,Tn_boot)
    
  }
  
  T=seq(0,1,l=m)
  Tn=int.simpson2(x=T,y=Tn_vec,method="CSR")
  Tn_boot=apply(Tn_boot_mat,1,
                function(z){int.simpson2(x=T,y=z,method="CSR")} )
  
  pvalue_vec[i]=mean(Tn_boot>=Tn)
  
  setTxtProgressBar(pb, i)
}
close(pb)



cat("\n Pvalues for 1%, 5% and 10%: \n")
c(mean(pvalue_vec<=0.01, na.rm=T),mean(pvalue_vec<=0.05, na.rm=T),mean(pvalue_vec<=0.1, na.rm=T))



