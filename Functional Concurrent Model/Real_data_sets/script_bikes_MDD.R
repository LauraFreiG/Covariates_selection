
#Script for bikes rent

#Libraries
library("fda.usc")
library(Rcpp)
library(RcppArmadillo)
library(MASS)
sourceCpp("signif_test_MDD_FCM.cpp")


##load bike data
data<-read.csv("bikenp.csv")
#####################################
# @ n=2512, p=4                     #
# @ casual = response
# @ temp = covariate
# @ atemp = covariate
# @ hum = covariate
# @ windspeed = covariate
#####################################
attach(data)
data$casual<-log(1+casual) 
n=length(unique(id))
m=length(unique(hour))
dim(data)[1]-m*n 


#Recovering missing data
t_seq=as.numeric(levels(as.factor(data$hour)))
i_ind=c()
y=X1=X2=X3=X4=matrix(NA,n,m)                 
for(i in 1:n){
  id_aux=i
  ind = data$id==id_aux
  
  if(sum(ind)==m){
    y[i,] = data$casual[ind]
    X1[i,] = data$temp[ind]
    X2[i, ] = data$atemp[ind]
    X3[i, ] = data$hum[ind]
    X4[i, ] = data$windspeed[ind]
  } else {
    #Interpolation
    i_ind=c(i_ind,i)
    t_seq_aux= data$hour[ind]
    ind_t = t_seq %in% t_seq_aux
    y[i,ind_t] = data$casual[ind]
    y[i,] = spline(t_seq,y[i,],xout=t_seq)$y
    X1[i,ind_t] = data$temp[ind]
    X1[i,] = spline(t_seq,X1[i,],xout=t_seq)$y
    X2[i,ind_t] = data$atemp[ind]
    X2[i,] = spline(t_seq,X2[i,],xout=t_seq)$y
    X3[i,ind_t] = data$hum[ind]
    X3[i,] = spline(t_seq,X3[i,],xout=t_seq)$y
    X4[i,ind_t] = data$windspeed[ind]
    X4[i,] = spline(t_seq,X4[i,],xout=t_seq)$y
  }
}




# Specification test MDD
#-------------------------
B_boot_n=1000


#Progress bar
pb <- txtProgressBar(max = length(t), style = 3) 
cat("\n")
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)


partial=TRUE # if partial (TRUE) or global (FALSE) test
var=3 #variable por partial test (var=1,2,3,4)


t=t_seq
sal=c()
for(t_ind in 1:length(t)){
    
  
  Y=as.matrix(y[,t_ind])
  if(partial){
    p=1
    if(var==1){
      X=as.matrix(X1[,t_ind])
    }else if(var==2){
      X=as.matrix(X1[,t_ind])
    }else if(var==3){
      X=as.matrix(X3[,t_ind])
    }else if(var==4){
      X=as.matrix(X4[,t_ind])
    }else{
      stop("Wrong variable choice")
    }
  }else{
    X=X=cbind(X1[,t_ind],X2[,t_ind],X3[,t_ind],X4[,t_ind])
    p=dim(X)[2]
  }

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
  # Tn_vec=c(Tn_vec,Tn)
  
  
  
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

