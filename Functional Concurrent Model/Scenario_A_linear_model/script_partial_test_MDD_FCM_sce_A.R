##########################################################
# Partial MDD significance tests for the synchronous FCM
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
        y=t(apply(x1,1,function(z){z=z*beta1; return(z)})) +
          t(apply(x2,1,function(z){z=z*beta2; return(z)})) +
        gp[(p*n+1):((p+1)*n),]
    }else{
      stop("Wrong hypothesis provided")
    }
    
    
    # Partial tests----
    Tn_vec_1=c()
    Tn_vec_2=c()
    Tn_boot_mat_1=c()
    Tn_boot_mat_2=c()
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
      
      suma_1=suma_2=0
      A1=A[[1]]
      A2=A[[2]]
      for(k in 1:(n-1)){
        for(l in (k+1):n){
          suma_1 = suma_1 + (A1[k,l]*A1[k,l]*(B[k,l]^2))
          suma_2 = suma_2 + (A2[k,l]*A2[k,l]*(B[k,l]^2))
        }
        
      }
      
      
      S1= (2/(n*(n-1)*cn))*suma_1
      S2= (2/(n*(n-1)*cn))*suma_2
      
      Tn_1=sqrt(choose(n,2))*MDD[1]/sqrt(S1)
      Tn_2=sqrt(choose(n,2))*MDD[2]/sqrt(S2)
      
      Tn_vec_1=c(Tn_vec_1,Tn_1)
      Tn_vec_2=c(Tn_vec_2,Tn_2)
      
      
      
      # Bootstrap---
      Tn_boot_1=c()
      Tn_boot_2=c()
      for(b in 1:B_boot_n){
        set.seed(b)
        e=rnorm(n)
        
        E1=matrix(e,n,n)
        E2=t(E1)
        B_boot=B*E1*E2
        MDD_boot_1=MDD_boot_2=0
        suma_boot_1= sum(A1*B_boot) - sum(diag(A1*B_boot))
        suma_boot_2= sum(A2*B_boot) - sum(diag(A2*B_boot))
        MDD_boot_1= suma_boot_1/(n*(n-1))
        MDD_boot_2= suma_boot_2/(n*(n-1))
        
        suma_1=suma_S_boot_parcial(e, n, 1, A1, B)
        suma_2=suma_S_boot_parcial(e, n, 1, A2, B)
        
        S_boot_1=suma_1/choose(n,2)
        S_boot_2=suma_2/choose(n,2)
        Tn_boot_1=c(Tn_boot_1, sqrt(choose(n,2))*MDD_boot_1/sqrt(S_boot_1) )
        Tn_boot_2=c(Tn_boot_2, sqrt(choose(n,2))*MDD_boot_2/sqrt(S_boot_2) )
        
      }
      
      
      Tn_boot_mat_1=cbind(Tn_boot_mat_1,Tn_boot_1)
      Tn_boot_mat_2=cbind(Tn_boot_mat_2,Tn_boot_2)
      
      
    }

    Tn_1=int.simpson2(x=t,y=Tn_vec_1,method="CSR")
    Tn_2=int.simpson2(x=t,y=Tn_vec_2,method="CSR")
    
    Tn_boot_1=apply(Tn_boot_mat_1,1,
                  function(z){int.simpson2(x=t,y=z,method="CSR")} )
    Tn_boot_2=apply(Tn_boot_mat_2,1,
                    function(z){int.simpson2(x=t,y=z,method="CSR")} )
    
    
    pvalue_vec1[m]=mean(Tn_boot_1>=Tn_1)
    pvalue_vec2[m]=mean(Tn_boot_2>=Tn_2)
    
    
    setTxtProgressBar(pb, m)
  }

close(pb)




cat("pvalues for X1 at levels 0.01, 0.05 and 0.1 are:\n",sep="")
cat(mean(pvalue_vec1<=0.01),", ",mean(pvalue_vec1<=0.05)," and ",mean(pvalue_vec1<=0.1),sep="")
cat("\n")
cat("pvalues for X2 at levels 0.01, 0.05 and 0.1 are:\n",sep="")
cat(mean(pvalue_vec2<=0.01),", ",mean(pvalue_vec2<=0.05)," and ",mean(pvalue_vec2<=0.1),sep="")











