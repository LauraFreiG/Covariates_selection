############################################
# Code for scenario 1: independence
############################################

data_gener<-function(p,n,s,scen){
  
# p=100 # nº covariates
# n=25 # sample size (n=25,50,100,150,300)
# s=10 # relevant terms 
# scen="a" # scenario (scen="a","b","c")
  

#beta vector
nbeta=s #first s non-zero coefficients
beta=c(rep(1.25,nbeta),rep(0,p-nbeta))


per=0.9 #explanation of the 90% of deviance
Sig=diag(1,p,p) #Sigma matrix
if(scen=="a"){
  dev=sqrt((1-per)/per*sum(beta^2))
}else if(scen=="b"){
  var_vec=c(0.5,0.5,1,1,3,3,10,10,25,25)
  diag(Sig)[1:nbeta]=var_vec
  var_sum= sum( (beta[1:s]^2)*var_vec )
  dev=sqrt((1-per)/per*var_sum)
}else if(scen=="c"){
  var_vec=c(0.5,0.5,1,1,3,3,10,10,25,25, 
            0.5,0.5,1.5,1.5,3,3,10,10,25,25,50,50)
  diag(Sig)[1:length(var_vec)]=var_vec
  var_sum= sum( (beta[1:s]^2)*var_vec[1:s] )
  dev=sqrt((1-per)/per*var_sum)
}else{
  stop("Wrong scenario")
}  


#data generation
library(MASS)
x=mvrnorm(n=n,mu=rep(0,p),Sigma=Sig)
x=x-matrix(colMeans(x),n,p,byrow=TRUE)
eps=rnorm(n=n, sd=dev)
y=x%*%beta+eps; y=as.numeric(y)
y=y-mean(y)


return(list("x"=x,"y"=y))}