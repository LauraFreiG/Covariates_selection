#n=200
#aa=c(1,1,1,1)
#Simulation model
#dat=data.frame(z1=runif(n),z2=rnorm(n),z3=rnorm(n),z4=runif(n),z5=rnorm(n))

#ap3=3*dat$z1
#ap4=dat$z2^2
#ap5=rnorm(n,sd=0.25)
#dat$y=10+aa[3]*ap3+aa[4]*ap4+ap5

# ldata contains in ldata$df the scalar variables and the rest of elements in the list are funcional data. df must be the first object in the list ldata.

gam.vs=function(resp,dat,alpha=0.05,limdcor=0.01){
library("mgcv")
library("fda.usc")
n=length(dat[[resp]])
nvar=length(names(dat))-1 #Global number of variates
namesc=names(dat)[which(names(dat)!=resp)] #covariates names
dcor=matrix(0,nrow=nvar,ncol=nvar) #nvar x nvar == p x p
ind=numeric(nvar) #vector indices
colnames(dcor)=namesc
names(ind)=colnames(dcor)
rownames(dcor)=1:nvar
form=as.formula(paste0(resp,"~1"))

j=1
Mset=c() #variates included in the model
Sset=namesc #variates which can be included
modelo=gam(formula=form,data=dat)
res=dat[,resp]-mean(dat[,resp]) #residue
modfinish=FALSE
sal=0 #si el mod gam da error que devuelva lo que tiene

while(!modfinish){
  for (i in 1:length(namesc)){
	  if (namesc[i] %in% Sset) {
	    if (is.factor(dat[[namesc[i]]])) {
	      y=model.matrix(as.formula(paste0("~-1+",namesc[i])),dat)
	    } else {
	      y=dat[[namesc[i]]]
	    }
	  tt=dcor.xy(res,y,n=length(res))
	  dcor[j,i]=tt$estimate*(tt$p.value<alpha)} else {dcor[j,i]=NA}
    #print(paste0("Variable:",names(ldata$df)[namesc[i]]))
  }

  #dcor[j,colnames(dcor) %in% Mset]=-dcor[j,colnames(dcor) %in% Mset]
  jj=which(dcor[j,]==max(dcor[j,],na.rm=TRUE))
  if (is.null(jj)) modfinish=TRUE
  #if(sum(dcor[j,jj])==0) {cat("Dcor<0.01","\n"); break}
  if(sum(dcor[j,jj])==0) { break}
    #if(length(jj)>1) cat(paste0("Iter:",j," dcor(",colnames(dcor)[jj[1]],")= ",round(dcor[j,jj[1]],4),"\n"))
    #if(length(jj)==1) cat(paste0("Iter:",j," dcor(",colnames(dcor)[jj],")= ",round(dcor[j,jj],4),"\n"))
  if (max(dcor[j,],na.rm=TRUE)>limdcor) {
		modant=modelo
		formant=form

		Mset=c(Mset,namesc[jj])
		Sset=setdiff(Sset,namesc[jj])
		if (is.factor(dat[,namesc[jj]])){
		  form=update.formula(form, paste0(".~.+",namesc[jj]))
		}else{
		  form=update.formula(form, paste0(".~.+s(",namesc[jj],")"))
		}
		
    j=j+1
  } else {modfinish=TRUE}
  
  sal=tryCatch(gam(formula=form,data=dat), error=function(e){return(TRUE)})
  if(is.logical(sal)) {break}
  #if(is.logical(sal)) {cat("problemas de convergencia ", "\n" ); break}
  #if(sal==TRUE) break
	modelo=sal
  #gam(formula=form,data=dat)
	aa=summary(modant)$sp.criterion
	ab=summary(modelo)$sp.criterion
  if (aa>0 & (ab/aa)>0.995){
	  form=formant
	  modelo=modant
	  Mset=setdiff(Mset,namesc[jj])
    #	modfinish=TRUE
  } else {res=residuals(modelo)}
	
}


i.predictor=rep(0,length(namesc))
if(length(Mset)==0){
  Mset=0
} else { Mset=setdiff(Mset,namesc[jj])}
i.predictor[which(namesc %in% Mset)]=1
#i.predictor[which(namesc %in% attr(modelo$terms,"term.labels"))]=1
names(i.predictor)=namesc
return( list(form=modelo$formula, data=dat, model=modelo,dcor=dcor,i.predictor=i.predictor, 
            indices=Mset) )
}

#final=gam.vs("y",dat)