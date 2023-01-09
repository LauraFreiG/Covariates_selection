############################################################
# Code for scenario 2: Toeplitz covariance with unit scales
############################################################

data_gener<-function(p,n,rho,scen){

# p=100 # nº covariates
# n=400 # sample size (n=25,50,100,200,400)
# rho=0.5 # rho value (rho=0.5,0.9)
# scen="a" # scenario (scen="a","b")
  
  

s=10
beta=rep(0,p)
per=0.9 #explanation of the 90% of deviance
Sig=diag(1,p,p) #Sigma matrix
if(scen=="a"){
  beta[seq(3,30,by=3)]=0.5
  var_vec=c(0.5, 0.5, 1, 1, 3, 3, 10, 10, 25, 25)
  seq_vec=seq(3,30,by=3)
  sum_total=0
  for(jjj in 1:(s-1)){  
    jjj=seq_vec[jjj]
    kkk=seq_vec[which(seq_vec>jjj)]
    sj=sqrt(var_vec[jjj/3])
    sk=sqrt(var_vec[kkk/3])
    sum_total=sum_total+sum(sj*sk*rho^abs(jjj-kkk))
  }
  dev=sqrt( ((1-per)/per)*(19.75+0.5*sum_total) )
  var_vec_aux=rep(1,p)
  var_vec_aux[seq_vec]=var_vec
  for(i in 1:p){
    for(j in i:p){
      Sig[i,j]=Sig[j,i]=sqrt(var_vec_aux[i])*sqrt(var_vec_aux[j])*(rho^{abs(i-j)})
    }
  }
}else if(scen=="b"){
  beta[seq(3,30,by=3)]=0.5
  var_vec=c(0.5, 0.5, 1, 1, 3, 3, 10, 10, 25, 25,
            0.5, 0.5, 1.5, 1.5, 3, 3, 10, 10, 25, 25, 50, 50)
  seq_vec=seq(3,30,by=3)
  sum_total=0
  for(jjj in 1:(s-1)){  
    jjj=seq_vec[jjj]
    kkk=seq_vec[which(seq_vec>jjj)]
    sj=sqrt(var_vec[jjj/3])
    sk=sqrt(var_vec[kkk/3])
    sum_total=sum_total+sum(sj*sk*rho^abs(jjj-kkk))
  }
  dev=sqrt( ((1-per)/per)*(19.75+0.5*sum_total) )
  ind=seq(3,30,by=3)
  ind2=seq(2,35,by=3)
  Sig=matrix(0,p,p)
  var_vec_aux=rep(1,p)
  var_vec_aux[ind]=var_vec[1:s] # relevant with different scales
  var_vec_aux[ind2]=var_vec[(s+1):22] # unrelevant with different scales
  for(i in 1:p){ 
    for(j in i:p){
      Sig[i,j]=Sig[j,i]=sqrt(var_vec_aux[i])*sqrt(var_vec_aux[j])*(rho^{abs(i-j)})
    }
  }
}else{
  stop("Wrong scenario")
}  


#data generation
x=mvrnorm(n=n,mu=rep(0,p),Sigma=Sig)
x=x-matrix(colMeans(x),n,p,byrow=TRUE)
eps=rnorm(n=n, sd=dev)
y=x%*%beta+eps; y=as.numeric(y)
y=y-mean(y)


return(list("x"=x,"y"=y))}