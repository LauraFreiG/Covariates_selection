##########################################################
# Global MDD significance test for the synchronous FCM
# Scenario A: linear model
##########################################################


#Simulation parameters
M=2000 # Monte-Carlo replicates
p=2
t=0:24
n=20 #sample size (n=20,40,60,80,100)
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
library(EDMeasure)
sourceCpp("signif_test_MDD_FCM.cpp")


#Progress bar
pb <- txtProgressBar(max = M, style = 3) #barra progreso
cat("\n")
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)



pvalue_vec=rep(0,M)
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
    
    
    # Global test----
    Tn_vec=c()
    Tn_boot_mat=c()
    for(t_ind in 1:length(t)){
      X=cbind(x1[,t_ind],x2[,t_ind])
      Y=y[,t_ind]
      
      # Tn calculation----
      MDD=apply(X,2,function(z){mdd(z,Y,compute="C",center="U")})
      cn= ((n-3)^4/(n-1)^4) + (2*(n-3)^4/((n-1)^4*(n-2)^3)) + (2*(n-3)/((n-1)^4*(n-2)^3))
      
      
      A=list()
      for(j in 1:p){
        X_aux=as.vector(X[,j])
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
    
    
    Tn=int.simpson2(x=t,y=Tn_vec,method="CSR")
    Tn_boot=apply(Tn_boot_mat,1,
                   function(z){int.simpson2(x=t,y=z,method="CSR")} )
    
    pvalue_vec[m]=mean(Tn_boot>=Tn)

    
    setTxtProgressBar(pb, m)
    
}
  


close(pb)



cat("pvalues for levels 0.01, 0.05 and 0.1 are:\n",sep="")
cat(mean(pvalue_vec<=0.01),", ",mean(pvalue_vec<=0.05)," and ",mean(pvalue_vec<=0.1),sep="")















