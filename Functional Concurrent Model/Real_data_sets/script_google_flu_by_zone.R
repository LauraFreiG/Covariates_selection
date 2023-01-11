######################################################
# Google flu data example
######################################################

#Read data
activity=read.csv(file = 'google_flu_activity.csv', header=T, row.names=1)
MDTV=read.csv(file = 'MDTV.csv', header=T, row.names=1)


midwest=c("ND","SD","KS","MN","IA","MO","IL","MI","IN","OH")
north_east=c("ME","NH","VT","MA","CT","NY","NJ","PA")
south=c("DE","MD","WV","VA","NC","KY","TN","AL","GA","MS","FL","LA","AR","TX")
west=c("WA","OR","ID","WY","CO","NM","AZ","UT","NV","CA")
midwest_loc=which(row.names(MDTV) %in% midwest)
north_east_loc=which(row.names(MDTV) %in% north_east)
south_loc=which(row.names(MDTV) %in% south)
west_loc=which(row.names(MDTV) %in% west)
n=dim(MDTV)[1]

#Create an adjacency matrix for the states in the US
states=as.vector(t(read.table("states.txt")))
states_use=states[-which(states %in% c("Montana","Nebraska","Oklahoma",
                                      "Wisconsin","Rhode Island","District of Columbia",
                                      "South Carolina","Oklahoma"))]

#Libraries
library(fda.usc)
library(Rcpp)
library(RcppArmadillo)
library(MASS)
sourceCpp("signif_test_MDD_FCM.cpp")


#Parameters
B_boot_n=1000
t=1:dim(MDTV)[2]

#Progress bar
pb <- txtProgressBar(max = length(t), style = 3) 
cat("\n")
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)


#Global dependence test
sal=vector("list", length(t))
for(t_ind in 1:length(t)){
  
  MDTV_aux=MDTV[,t_ind]
  naindex=which(is.na(MDTV[t_ind]))
  #by zone
  MDTV_aux_mid=MDTV_aux[midwest_loc[!midwest_loc %in% naindex]]
  Y_mid=as.matrix(activity[midwest_loc[!midwest_loc %in% naindex],t_ind])
  MDTV_aux_ne=MDTV_aux[north_east_loc[!north_east_loc %in% naindex]]
  Y_ne=as.matrix(activity[north_east_loc[!north_east_loc %in% naindex],t_ind])
  MDTV_aux_s=MDTV_aux[south_loc[!south_loc %in% naindex]]
  Y_s=as.matrix(activity[south_loc[!south_loc %in% naindex],t_ind])
  MDTV_aux_w=MDTV_aux[west_loc[!west_loc %in% naindex]]
  Y_w=as.matrix(activity[west_loc[!west_loc %in% naindex],t_ind])
  
  MDTV_aux_mid=scale(MDTV_aux_mid,center=FALSE,scale=TRUE)
  Y_mid=scale(Y_mid,center=FALSE,scale=TRUE)
  MDTV_aux_ne=scale(MDTV_aux_ne,center=FALSE,scale=TRUE)
  Y_ne=scale(Y_ne,center=FALSE,scale=TRUE)
  MDTV_aux_s=scale(MDTV_aux_s,center=FALSE,scale=TRUE)
  Y_s=scale(Y_s,center=FALSE,scale=TRUE)
  MDTV_aux_w=scale(MDTV_aux_w,center=FALSE,scale=TRUE)
  Y_w=scale(Y_w,center=FALSE,scale=TRUE)
  
  X_list=list(MDTV_aux_mid,MDTV_aux_ne,MDTV_aux_s,MDTV_aux_w)
  Y_list=list(Y_mid,Y_ne,Y_s,Y_w)
  names(X_list)=names(Y_list)=c("midwest","north_east","south","west")


  n_vec=c( length(Y_mid), length(Y_ne), length(Y_s), length(Y_w)  )
  p_vec=rep(1,4)
  library(EDMeasure)

  
  Tn_aux=c()
  for(z in 1:4){
    X=X_list[[z]]
    Y=Y_list[[z]]
    n=n_vec[z]
    p=p_vec[z]
    
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
    
    Tn_aux=rbind( Tn_aux, c(Tn,Tn_boot) )
  }
  
  
  
  sal[[t_ind]]=Tn_aux

  
  setTxtProgressBar(pb, t_ind)
  
}

cat("\n")
close(pb)


Tn_vec_mid=Tn_vec_ne=Tn_vec_s=Tn_vec_w=c()
Tn_boot_mat_mid=Tn_boot_mat_ne=Tn_boot_mat_s=Tn_boot_mat_w=c()
for(i in 1:length(t)){
  sal_aux=sal[[i]]
  
  Tn_vec_mid=c(Tn_vec_mid,sal_aux[1,1])
  Tn_vec_ne=c(Tn_vec_ne,sal_aux[2,1])
  Tn_vec_s=c(Tn_vec_s,sal_aux[3,1])
  Tn_vec_w=c(Tn_vec_w,sal_aux[4,1])
  
  Tn_boot_mat_mid=rbind( Tn_boot_mat_mid, sal_aux[1,-1] )
  Tn_boot_mat_ne=rbind( Tn_boot_mat_ne, sal_aux[2,-1] )
  Tn_boot_mat_s=rbind( Tn_boot_mat_s, sal_aux[3,-1] )
  Tn_boot_mat_w=rbind( Tn_boot_mat_w, sal_aux[4,-1] )

}
Tn_mid=int.simpson2(x=t,y=Tn_vec_mid,method="CSR")
Tn_ne=int.simpson2(x=t,y=Tn_vec_ne,method="CSR")
Tn_s=int.simpson2(x=t,y=Tn_vec_s,method="CSR")
Tn_w=int.simpson2(x=t,y=Tn_vec_w,method="CSR")


Tn_boot_mid=apply(Tn_boot_mat_mid,2,
              function(z){int.simpson2(x=t,y=z,method="CSR")} )
Tn_boot_ne=apply(Tn_boot_mat_ne,2,
                  function(z){int.simpson2(x=t,y=z,method="CSR")} )
Tn_boot_s=apply(Tn_boot_mat_s,2,
                  function(z){int.simpson2(x=t,y=z,method="CSR")} )
Tn_boot_w=apply(Tn_boot_mat_w,2,
                  function(z){int.simpson2(x=t,y=z,method="CSR")} )


pvalues=c( mean(Tn_boot_mid>=Tn_mid), mean(Tn_boot_ne>=Tn_ne),
           mean(Tn_boot_s>=Tn_s), mean(Tn_boot_w>=Tn_w))

cat("pvalues=", pvalues, "\n", sep=" ")









