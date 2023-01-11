######################################################
# Gait data example
######################################################

library(fda)
library(fda.usc)
Hip=t(gait[,,"Hip Angle"])
Knee=t(gait[,,"Knee Angle"])
n=dim(Hip)[1]
p=1
t=(1:dim(Hip)[2])/dim(Hip)[2]
B_boot_n=1000



# Progress bar
pb <- txtProgressBar(max = length(t), style = 3) #
cat("\n")
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)


# libraries
library(Rcpp)
library(RcppArmadillo)
library(MASS)
sourceCpp("signif_test_MDD_FCM.cpp")
  
sal=c()
for(t_ind in 1:length(t)){  
  
  
  Y=as.matrix(Knee[,t_ind])
  X=as.matrix(Hip[,t_ind])
  
  
  # Tn calculation----
  library(EDMeasure)
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
  
  
  
  # Bootstrap---
  Tn_boot=c()
  S_boot=c()
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
  
  sal=rbind(sal,c(Tn,Tn_boot))
  
  setTxtProgressBar(pb, t_ind)
  
}

cat("\n")
close(pb)

Tn_vec=sal[,1]
Tn_boot_mat=sal[,-1]
Tn=int.simpson2(x=t,y=Tn_vec,method="CSR")
Tn_boot=apply(Tn_boot_mat,2,
              function(z){int.simpson2(x=t,y=z,method="CSR")} )

pvalue=mean(Tn_boot>=Tn)

cat("pvalue=", pvalue, "\n", sep="")



