########################################################################
# Code to implement the scenario A of the FLCM paper using the MDD test
########################################################################


# Parameters
n=60 # sample size (n=60,100)
dd=0  # d value (dd=0,3,7)
p=1
m=81
T<-seq(0,1,l=m)
sigma<-.9
M=1000
B_boot_n=200 


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
library(Rcpp)
library(RcppArmadillo)
library(MASS)
library(fda.usc)
sourceCpp("signif_test_MDD_FCM.cpp")


# MDD test----
Tn_vec=c()
Tn_boot_mat=c()
for(t_ind in 1:length(T)){
  X=as.matrix(Xmat[,t_ind])
  Y=Ymat[,t_ind]
  
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


Tn=int.simpson2(x=T,y=Tn_vec,method="CSR")
Tn_boot=apply(Tn_boot_mat,1,
              function(z){int.simpson2(x=T,y=z,method="CSR")} )

pvalue=mean(Tn_boot>=Tn)

pvalues=c(pvalues,pvalue)
setTxtProgressBar(pb,m_ind)
}
close(pb)

cat("pvalues for levels 0.01, 0.05 and 0.1 are:\n",sep="")
c(mean(pvalues<=0.01),mean(pvalues<=0.05),mean(pvalues<=0.1))
