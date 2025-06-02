

Sig=diag(p)
aux=eigen(Sig)
aux$values


# scenario 2 (UTOEP)
#--------------------------------------------------------
library(MASS)
p=100
rho=0.9 # rho=0.5,0.9


Sig=matrix(0,p,p)
values=rho^{0:(p-1)}
for(i in 1:p){
  Sig[i,1:p>=i]=Sig[1:p>=i,i]=values[1:sum(1:p>=i)]
}

# UTOEP-B
Sig_aux=Sig[1:15,1:15]
aux=eigen(Sig_aux)$values
cumsum(aux)/sum(aux)*100

# UTOEP-S
Sig_aux=Sig[seq(3,30,by=3),seq(3,30,by=3)]
aux=eigen(Sig_aux)$values
cumsum(aux)/sum(aux)*100



# scenario 3 (TOEP)
#--------------------------------------------------------
library(MASS)
p=100
rho=0.9 # rho=0.5,0.9

Sig=matrix(0,p,p)
values=rho^{0:(p-1)}
for(i in 1:p){
  Sig[i,1:p>=i]=Sig[1:p>=i,i]=values[1:sum(1:p>=i)]
}
Sig_UTOEP=Sig

# RC.TOEP-S
D=diag(p)
diag(D)[seq(3,30,by=3)]=sqrt( c(0.5,0.5,1,1,3,3,10,10,25,25) )
Sig=D%*%Sig_UTOEP%*%t(D)

Sig_aux=Sig[seq(3,30,by=3),seq(3,30,by=3)]
aux=eigen(Sig_aux)$values
cumsum(aux)/sum(aux)*100


# RNC.TOEP-S
D=diag(p)
diag(D)[seq(3,30,by=3)]=sqrt( c(0.5,0.5,1,1,3,3,10,10,25,25) )
diag(D)[c(2,5,8,11,14,17,20,23,26,29,32,35)]=sqrt( c(0.5,0.5,1.5,1.5,3,3,10,10,25,25,50,50) )
Sig=D%*%Sig_UTOEP%*%t(D)

Sig_aux=Sig[seq(3,30,by=3),seq(3,30,by=3)]
aux=eigen(Sig_aux)$values
cumsum(aux)/sum(aux)*100



