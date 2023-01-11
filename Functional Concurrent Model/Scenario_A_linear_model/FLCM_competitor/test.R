##########test functions for dense data, testing for global effect of x#######
#the GetMatrices R function to get the stacked matrices from the input data ymat,xmat,nbas
GetMatrices <- function(ymat, xmat, nbas)
{
  GrandY <- as.vector(t(ymat))
  n = nrow(ymat)
  p = 1
  m = ncol(xmat)
  T <- seq(0, 1, l = m)
  #Using FPCA to preprocess noisy data
  w.sm = fpca.sc(xmat, pve = .99, var = TRUE)
  xhat = w.sm$Yhat
  #equispaced knots used, one can used alternate knots on quantiles or other config based on domain knowldge
  knots0 = seq(min(T),max(T),l=5)     #determine based on T
  knots1 = seq(min(T),max(T),l=nbas-2)
  norder = 4
  nbasis0 = length(knots0) + norder - 2
  nbasis1 = length(knots1) + norder - 2
  dayrng = c(0,1)
  bbasis0 = create.bspline.basis(dayrng,nbasis0,norder,knots0)
  bbasis1 = create.bspline.basis(dayrng,nbasis1,norder,knots1)
  bbasisMat0 = eval.basis(T,bbasis0)
  bbasisMat1 = eval.basis(T,bbasis1)
  in.mat0 = inprod(bbasis0,bbasis0) #penalty matrix for function
  in.mat1 = inprod(bbasis1,bbasis1)
  Kphi0<-in.mat0
  basismat0<-bbasisMat0
  colnames(basismat0)<-NULL
  Kphi1<-in.mat1
  basismat1<-bbasisMat1
  colnames(basismat1)<-NULL
  XMatstar <- array(0, dim = c(n, p, m, nbasis1))
  {
    for (i in 1:n)
    {
      for (j in 1:p)
        for (l in 1:m) {
          XMatstar[i, j, l, ] = xhat[i, l] * basismat1[l, ]
        }
    }
  }
  B<-basismat0
  R0<-chol(Kphi0)
  L0<-t(R0)
  M0<-solve(R0)
  R1<-chol(Kphi1)
  L1<-t(R1)
  M1<-solve(R1)
  ZMatstar <- array(0, dim = c(n, p, m, nbasis1))
  {
    for (i in 1:n)
    {
      for (j in 1:p)
        for (l in 1:m) {
          ZMatstar[i, j, l, ] = XMatstar[i, j, l, ] %*% M1
        }
    }
  }
  Zimatlist1 <- list()
  for (i in 1:n)
  {
    Zimatlist1[[i]] <- matrix(0, m, nbasis1)
    for (l in 1:m)
    {
      Zimatlist1[[i]][l, 1:nbasis1] <- ZMatstar[i, 1, l, ]
    }
  }
  B_repara <- B %*% M0
  Grandz1 <- NULL
  for (i in 1:n)
  {
    Grandz1 <- rbind(Grandz1, Zimatlist1[[i]])
  }
  GrandB <- NULL
  for (i in 1:n)
  {
    GrandB <- rbind(GrandB, B_repara)
  }
  matlist <- list(GrandY, GrandB, Grandz1, T, n, m)
  names(matlist) <- c("GrandY", "GrandB", "Grandz1", "T", "n", "m")
  return(matlist)
}

###scorefuncpval R Function for testing: Requires Matrices in stacked form
scorefuncpval <- function(GrandY, GrandB, Grandz1, T, n, m) {
  Rcpp::sourceCpp('test.cpp')
  nbasis0<-ncol(GrandB)
  nbasis1<-ncol(Grandz1)
  PP <- list(GrandB = list(diag(nbasis0)), Grandz1 = list(diag(nbasis1)))
  #Fitting the full model to get the residuals
  fit <- bam(GrandY ~ -1 + GrandB + Grandz1, paraPen = PP, method = "REML")
  res <- fit$residuals
  resmat <- matrix(res, n, m, byrow = TRUE)
  #Using FPCA to get estimate of the covariance matrix
  fpca1 <- fpca.sc(resmat,center = TRUE,nbasis = 15,argvals = T,pve = 0.999,var = TRUE)
  # fpca1 <- fpca.sc(resmat,center = TRUE,argvals = T,pve = 0.999,var = TRUE)
  if (length(fpca1$evalues) > 1)
  {
    sigmat <- fpca1$efunctions %*% diag(fpca1$evalues) %*% t(fpca1$efunctions) + fpca1$sigma2 *diag(m)
  }
  if (length(fpca1$evalues) == 1)
  {
    sigmat <-fpca1$evalues * fpca1$efunctions %*% t(fpca1$efunctions) + fpca1$sigma2 *diag(m)
  }
  R <- chol(sigmat)
  L <- t(R)
  Linv <- forwardsolve(L, diag(dim(L)[1]))
  rm(list = setdiff(
    ls(),
    list("GrandY", "GrandB", "Grandz1", "Linv", "n", "m", "p")
  ))
  llist <- replicate(n, Linv, simplify = FALSE)
  Linvhalf <- bdiag(llist)
  NewB <- as.matrix(Linvhalf %*% GrandB)
  Newz1 <- as.matrix(Linvhalf %*% Grandz1)
  NewY <- as.matrix(Linvhalf %*% GrandY)
  Rcpp::sourceCpp('test.cpp')
  nbasis0<-ncol(NewB)
  nbasis1<-ncol(Newz1)
  PPreduced <- list(NewB = list(diag(nbasis0)))
  #Fitting the reduced model under null
  fitfinalreduced <- bam(NewY ~ -1 + NewB, paraPen = PPreduced, method = "REML")
  mvar <- (1 / fitfinalreduced$sp)
  #Getting the MLE of the variance components under null
  estt0 <- mvar[[1]]
  #Gettng the appropiate statistics for performing the test
  Wc <- sqrt(estt0) * NewB
  meat <- diag(ncol(Wc)) + eigenMapMatMult(t(Wc), Wc)
  meatinv <- chol2inv(chol(meat))
  Yvi <- NewY - eigenMapMatMult(eigenMapMatMult(eigenMapMatMult(Wc, meatinv), t(Wc)), NewY)
  Eigenmat <- eigenMapMatMult(t(Newz1), Newz1) - eigenMapMatMult(eigenMapMatMult(t(Newz1), eigenMapMatMult(eigenMapMatMult(Wc, meatinv), t(Wc))), Newz1)
  lambda <- eigen(Eigenmat, symmetric = TRUE)$values
  s1l <- sum(lambda)
  s2l <- .5 * sum((lambda) ^ 2)
  term1 <- eigenMapMatMult(t(Newz1), Yvi)
  Score <- -.5 * (s1l - eigenMapMatMult(t(term1), term1))
  EigenmatB <- eigenMapMatMult(t(NewB), NewB) - eigenMapMatMult(eigenMapMatMult(t(NewB), eigenMapMatMult(eigenMapMatMult(Wc, meatinv), t(Wc))), NewB)
  lambdaB <- eigen(EigenmatB, symmetric = TRUE)$values
  s2lB <- .5 * sum((lambdaB) ^ 2)
  EigenmatBz1 <- eigenMapMatMult(t(Newz1), NewB) - eigenMapMatMult(eigenMapMatMult(t(Newz1), eigenMapMatMult(eigenMapMatMult(Wc, meatinv), t(Wc))), NewB)
  trac <- function(A) {
    sum(diag(A))
  }
  s2lBz1 <- .5 * trac(eigenMapMatMult(EigenmatBz1, t(EigenmatBz1)))
  #The information matix
  Inform <- matrix(0, 2, 2)
  Inform[1, 1] <- s2lB
  Inform[2, 2] <- s2l
  Inform[1, 2] <- s2lBz1
  Inform[2, 1] <- s2lBz1
  I <- Inform[2, 2] - t(Inform[1, 2]) %*% solve(Inform[1, 1]) %*% Inform[1, 2]
  if (Score > 0)
  {
    Tsobs <- (t(Score) %*% Score) / I
  }
  if (Score <= 0)
  {
    Tsobs = 0
  }
  #The observed test statistic
  Tobs <- as.numeric(Tsobs)
  #Calculating the P-Value
  simuT <- function()
  {
    w <- rnorm(ncol(Newz1))
    s <- sum(lambda * w ^ 2) - s1l
    T1 <- (1/4) * (s) ^ 2
    if (s > 0)
    {
      return((T1/I))
    }
    if (s <= 0)
    {
      return(0)
    }
  }
  Tsnull <- replicate(100000, simuT())
  pval <- (sum(Tsnull >= Tobs)) / 100000
  pvallist <- list(Tobs, Score, pval)
  names(pvallist) <- c("Tobs", "Score", "pval")
  return(pvallist)
}

### The wrapper function testpval which uses the above two functions 
### Reprots P-Value of the test
FLCM.test1 <- function(ymat, xmat, nbas = 7)
{
  Matlist <- GetMatrices(ymat, xmat, nbas) # Getting the stacked matrices
  GrandY <- Matlist$GrandY
  GrandB <- Matlist$GrandB
  Grandz1 <- Matlist$Grandz1
  T <- Matlist$T
  n <- Matlist$n
  m <- Matlist$m
  scoremethodresult <- scorefuncpval(GrandY, GrandB, Grandz1, T, n, m) # Performing the test
  scorepval <- scoremethodresult$pval
  scorepval
}
########################################################################################
##############test functions for dense data,testing for x2 in presence of x1 ###########
GetMatrices2 <- function(ymat, x1mat, x2mat, nbas)
{
  
  GrandY <- as.vector(t(ymat))
  n = nrow(ymat)
  p = 2
  m = ncol(x1mat)
  T <- seq(0, 1, l = m)
  #Using FPCA to preprocess noisy data
  w.sm = fpca.sc(x1mat, pve = .99, var = TRUE)
  xhat = w.sm$Yhat
  w.sm2 = fpca.sc(x2mat, pve = .99, var = TRUE)
  xhat2 = w.sm2$Yhat
  #equispaced knots used, one can used alternate knots on quantiles or other config based on domain knowldge
  library(fda)
  #int
  knots0 = seq(min(T), max(T), l = 5)
  knots1 = seq(min(T), max(T), l = 5)
  knots2 = seq(min(T), max(T), l = nbas - 2)
  # We'll use fourth-order B-splines
  norder = 4
  # this implies the number of basis functions
  nbasis0 = length(knots0) + norder - 2
  nbasis1 = length(knots1) + norder - 2
  nbasis2 = length(knots2) + norder - 2
  dayrng = c(0, 1)
  bbasis0 = create.bspline.basis(dayrng, nbasis0, norder, knots0)
  bbasis1 = create.bspline.basis(dayrng, nbasis1, norder, knots1)
  bbasis2 = create.bspline.basis(dayrng, nbasis2, norder, knots2)
  #bbasis$nbasis    # number of basis functions
  #bbasis$rangeval   # basis range
  bbasisMat0 = eval.basis(T, bbasis0)
  bbasisMat1 = eval.basis(T, bbasis1)
  bbasisMat2 = eval.basis(T, bbasis2)
  #dim(bbasisMat)
  #100 102
  in.mat0 = inprod(bbasis0, bbasis0) #penalty matrix for function
  in.mat1 = inprod(bbasis1, bbasis1)
  in.mat2 = inprod(bbasis2, bbasis2)
  
  #par(mfrow=c(1,1))
  #image(in.mat)
  Kphi0 <- in.mat0
  basismat0 <- bbasisMat0
  colnames(basismat0) <- NULL
  Kphi1 <- in.mat1
  basismat1 <- bbasisMat1
  colnames(basismat1) <- NULL
  Kphi2 <- in.mat2
  basismat2 <- bbasisMat2
  colnames(basismat2) <- NULL
  XMatstar <- array(0, dim = c(n, 1, m, nbasis1))
  {
    for (i in 1:n)
    {
      for (j in 1:1)
        for (l in 1:m) {
          XMatstar[i, j, l, ] = xhat[i, l] * basismat1[l, ]
        }
    }
  }
  XMatstar2 <- array(0, dim = c(n, 1, m, nbasis2))
  {
    for (i in 1:n)
    {
      for (j in 1:1)
        for (l in 1:m) {
          XMatstar2[i, j, l, ] = xhat2[i, l] * basismat2[l, ]
        }
    }
  }
  R0 <- chol(Kphi0)
  L0 <- t(R0)
  M0 <- solve(R0)
  R1 <- chol(Kphi1)
  L1 <- t(R1)
  M1 <- solve(R1)
  R2 <- chol(Kphi2)
  L2 <- t(R2)
  M2 <- solve(R2)
  
  ZMatstar <- array(0, dim = c(n, 1, m, nbasis1))
  {
    for (i in 1:n)
    {
      for (j in 1:1)
        for (l in 1:m) {
          ZMatstar[i, j, l, ] = XMatstar[i, 1, l, ] %*% M1
        }
    }
  }
  Zimatlist1 <- list()
  for (i in 1:n)
  {
    Zimatlist1[[i]] <- matrix(0, m, nbasis1)
    for (l in 1:m)
    {
      Zimatlist1[[i]][l, 1:nbasis1] <- ZMatstar[i, 1, l, ]
    }
  }
  ZMatstar2 <- array(0, dim = c(n, 1, m, nbasis2))
  {
    for (i in 1:n)
    {
      for (j in 1:1)
        for (l in 1:m) {
          ZMatstar2[i, j, l, ] = XMatstar2[i, 1, l, ] %*% M2
        }
    }
  }
  Zimatlist2 <- list()
  for (i in 1:n)
  {
    Zimatlist2[[i]] <- matrix(0, m, nbasis2)
    for (l in 1:m)
    {
      Zimatlist2[[i]][l, 1:nbasis2] <- ZMatstar2[i, 1, l, ]
    }
  }
  B <- basismat0
  B_repara <- B %*% M0
  Grandz1 <- NULL
  for (i in 1:n)
  {
    Grandz1 <- rbind(Grandz1, Zimatlist1[[i]])
  }
  Grandz2 <- NULL
  for (i in 1:n)
  {
    Grandz2 <- rbind(Grandz2, Zimatlist2[[i]])
  }
  GrandB <- NULL
  for (i in 1:n)
  {
    GrandB <- rbind(GrandB, B_repara)
  }
  matlist <- list(GrandY, GrandB, Grandz1, Grandz2, T, n, m)
  names(matlist) <- c("GrandY", "GrandB", "Grandz1", "Grandz2", "T", "n", "m")
  return(matlist)
}

###scorefuncpval2 R Function for testing: Requires Matrices in stacked form
scorefuncpval2 <- function(GrandY, GrandB, Grandz1, Grandz2, T, n, m) {
  Rcpp::sourceCpp('test.cpp')
  nbasis0 <- ncol(GrandB)
  nbasis1 <- ncol(Grandz1)
  nbasis2 <- ncol(Grandz2)
  library(mgcv)
  PP <-
    list(
      GrandB = list(diag(nbasis0)),
      Grandz1 = list(diag(nbasis1)),
      Grandz2 = list(diag(nbasis2))
    )
  fit <-
    bam(GrandY ~ -1 + GrandB + Grandz1 + Grandz2,
        paraPen = PP,
        method = "REML")
  res <- fit$residuals
  library(refund)
  resmat <- matrix(res, n, m, byrow = TRUE)
  fpca1 <-
    fpca.sc(
      resmat,
      center = TRUE,
      nbasis = 15,
      argvals = T,
      pve = 0.999,
      var = TRUE
    )
  if (length(fpca1$evalues) > 1)
  {
    sigmat <-
      fpca1$efunctions %*% diag(fpca1$evalues) %*% t(fpca1$efunctions) + fpca1$sigma2 *
      diag(m)
  }
  if (length(fpca1$evalues) == 1)
  {
    sigmat <-
      fpca1$evalues * fpca1$efunctions %*% t(fpca1$efunctions) + fpca1$sigma2 *
      diag(m)
  }
  R <- chol(sigmat)
  L <- t(R)
  Linv <- forwardsolve(L, diag(dim(L)[1]))
  rm(list = setdiff(
    ls(),
    list("GrandY", "GrandB", "Grandz1", "Grandz2", "Linv", "n", "m", "p")
  ))
  llist <- replicate(n, Linv, simplify = FALSE)
  Linvhalf <- bdiag(llist)
  #see<- eigenMapMatMult(eigenMapMatMult((Linvhalf),Grandsigmat),t(Linvhalf))
  NewB <- as.matrix(Linvhalf %*% GrandB)
  Newz1 <- as.matrix(Linvhalf %*% Grandz1)
  Newz2 <- as.matrix(Linvhalf %*% Grandz2)
  NewY <- as.matrix(Linvhalf %*% GrandY)
  Rcpp::sourceCpp('test.cpp')
  #PPfull <-list(NewB=list(diag(102)),Newz1=list(diag(102)),Newz2=list(diag(102)))
  #fitfinalfull <-gam(NewY~ -1+NewB+Newz1+Newz2,paraPen=PP2,method="REML")
  nbasis0 <- ncol(NewB)
  nbasis1 <- ncol(Newz1)
  nbasis2 <- ncol(Newz2)
  PPreduced <-
    list(NewB = list(diag(nbasis0)), Newz1 = list(diag(nbasis1)))
  fitfinalreduced <-
    bam(NewY ~ -1 + NewB + Newz1, paraPen = PPreduced, method = "REML")
  mvar <-
    (1 / fitfinalreduced$sp)      #gam.vcomp(fitfinalreduced) #use 1/sp
  #w<-m^2
  estt0 <- mvar[[1]]
  estt1 <- mvar[[2]]
  Wc <- cbind(sqrt(estt0) * NewB, sqrt(estt1) * Newz1)
  meat <- diag(ncol(Wc)) + eigenMapMatMult(t(Wc), Wc)
  meatinv <- chol2inv(chol(meat))
  Yvi <-
    NewY - eigenMapMatMult(eigenMapMatMult(eigenMapMatMult(Wc, meatinv), t(Wc)), NewY)
  Eigenmat <-
    eigenMapMatMult(t(Newz2), Newz2) - eigenMapMatMult(eigenMapMatMult(t(Newz2), eigenMapMatMult(eigenMapMatMult(Wc, meatinv), t(Wc))), Newz2)
  lambda <- eigen(Eigenmat)$values
  s1l <- sum(lambda)
  s2l <- .5 * sum((lambda) ^ 2)
  term1 <- eigenMapMatMult(t(Newz2), Yvi)
  Score <- -.5 * (s1l - eigenMapMatMult(t(term1), term1))
  EigenmatB <-
    eigenMapMatMult(t(NewB), NewB) - eigenMapMatMult(eigenMapMatMult(t(NewB), eigenMapMatMult(eigenMapMatMult(Wc, meatinv), t(Wc))), NewB)
  lambdaB <- eigen(EigenmatB)$values
  s2lB <- .5 * sum((lambdaB) ^ 2)
  Eigenmatz1 <-
    eigenMapMatMult(t(Newz1), Newz1) - eigenMapMatMult(eigenMapMatMult(t(Newz1), eigenMapMatMult(eigenMapMatMult(Wc, meatinv), t(Wc))), Newz1)
  lambdaz1 <- eigen(Eigenmatz1)$values
  s2lz1 <- .5 * sum((lambdaz1) ^ 2)
  EigenmatBz1 <-
    eigenMapMatMult(t(Newz1), NewB) - eigenMapMatMult(eigenMapMatMult(t(Newz1), eigenMapMatMult(eigenMapMatMult(Wc, meatinv), t(Wc))), NewB)
  EigenmatBz2 <-
    eigenMapMatMult(t(Newz2), NewB) - eigenMapMatMult(eigenMapMatMult(t(Newz2), eigenMapMatMult(eigenMapMatMult(Wc, meatinv), t(Wc))), NewB)
  Eigenmatz2z1 <-
    eigenMapMatMult(t(Newz1), Newz2) - eigenMapMatMult(eigenMapMatMult(t(Newz1), eigenMapMatMult(eigenMapMatMult(Wc, meatinv), t(Wc))), Newz2)
  trac <- function(A) {
    sum(diag(A))
  }
  s2lBz1 <- .5 * trac(eigenMapMatMult(EigenmatBz1, t(EigenmatBz1)))
  s2lBz2 <- .5 * trac(eigenMapMatMult(EigenmatBz2, t(EigenmatBz2)))
  s2lz2z1 <- .5 * trac(eigenMapMatMult(Eigenmatz2z1, t(Eigenmatz2z1)))
  Inform <- matrix(0, 3, 3)
  Inform[1, 1] <- s2lB
  Inform[2, 2] <- s2lz1
  Inform[3, 3] <- s2l
  Inform[1, 2] <- s2lBz1
  Inform[2, 1] <- s2lBz1
  Inform[1, 3] <- s2lBz2
  Inform[3, 1] <- s2lBz2
  Inform[2, 3] <- s2lz2z1
  Inform[3, 2] <- s2lz2z1
  I <-
    Inform[3, 3] - t(Inform[1:2, 3]) %*% solve(Inform[1:2, 1:2]) %*% Inform[1:2, 3]
  if (Score > 0)
  {
    Tsobs <- (t(Score) %*% Score) / I
  }
  if (Score <= 0)
  {
    Tsobs = 0
  }
  Tobs <- as.numeric(Tsobs)
  simuT <- function()
  {
    w <- rnorm(nbasis2)
    s <- sum(lambda * w ^ 2) - s1l
    T1 <- (1 / 4) * (s) ^ 2
    if (s > 0)
    {
      return((T1 / I))
    }
    if (s <= 0)
    {
      return(0)
    }
  }
  Tsnull <- replicate(100000, simuT())
  pval <- (sum(Tsnull >= Tobs)) / 100000
  pvallist <- list(Tobs, Score, pval)
  names(pvallist) <- c("Tobs", "Score", "pval")
  return(pvallist)
}

### The wrapper function testpval which uses the above two functions 
### Reprots P-Value of the test
FLCM.test2 <- function(ymat, x1mat, x2mat, nbas = 7)
{
  Matlist <- GetMatrices2(ymat, x1mat, x2mat, nbas) # Getting the stacked matrices
  GrandY <- Matlist$GrandY
  GrandB <- Matlist$GrandB
  Grandz1 <- Matlist$Grandz1
  Grandz2 <- Matlist$Grandz2
  T <- Matlist$T
  n <- Matlist$n
  m <- Matlist$m
  scoremethodresult <- scorefuncpval2(GrandY, GrandB, Grandz1, Grandz2, T, n, m) # Performing the test
  scorepval <- scoremethodresult$pval
  scorepval
}
#####################################################################################
##########test functions for sparse data,testing effect of single covariate x1 ###########
#Function for getting stacked matrices
GetMatricesSparse <- function(y, x1, T, n, ID, nbas)
{
  UT <- c()
  for (i in 1:n)
  {
    UT <- union(UT, T[[i]])
  }
  UT <- sort(UT)
  Y <- list()
  for (i in 1:n)
  {
    ind <- which(ID == i)
    Y[[i]] <- y[ind]
  }
  GrandY <- NULL
  for (i in 1:n)
  {
    GrandY <- c(GrandY, Y[[i]])
  }
  W <- list()
  for (i in 1:n)
  {
    ind <- which(ID == i)
    W[[i]] <- x1[ind]  
  }
  lenX <- c()
  for (i in 1:n)
  {
    lenX[i] <- length(T[[i]])
  }
  .id = rep(1:n, lenX)
  .index <- unlist(T)
  Wdata <- data.frame(.id , .index)
  Wdata$.value <- unlist(W)
  #Using FPCA to preprocess noisy data
  w.sm = fpca.sc(ydata = Wdata, pve = .99, var = TRUE)
  xhat = w.sm$Yhat
  timepoints <- as.numeric(attr(w.sm$Yhat, "dimnames")[[2]])
  Xhat_resp <- list()
  for (i in 1:n)
  {
    Xhat_resp[[i]] <- vector(mode = "logical", length = length(T[[i]]))
    for (l in 1:length(T[[i]])) {
      ind <- which(round(timepoints, 8) == round(T[[i]], 8)[l])
      Xhat_resp[[i]][l] = xhat[i, ind]
    }
  }
  #equispaced knots used, one can used alternate knots on quantiles or other config based on domain knowldge
  knots0 = seq(min(UT), max(UT), l = 5)
  knots1 = seq(min(UT), max(UT), l = nbas - 2)
  # We'll use fourth-order B-splines
  norder = 4
  # this implies the number of basis functions
  nbasis0 = length(knots0) + norder - 2
  nbasis1 = length(knots1) + norder - 2
  dayrng = c(0, 1)
  bbasis0 = create.bspline.basis(dayrng, nbasis0, norder, knots0)
  bbasis1 = create.bspline.basis(dayrng, nbasis1, norder, knots1)
  #bbasis$nbasis    # number of basis functions
  #bbasis$rangeval   # basis range
  bbasisMat0 <- list()
  for (i in 1:n)
  {
    bbasisMat0[[i]] = eval.basis(T[[i]], bbasis0)
  }
  bbasisMat1 <- list()
  for (i in 1:n)
  {
    bbasisMat1[[i]] = eval.basis(T[[i]], bbasis1)
  }
  in.mat0 = inprod(bbasis0, bbasis0) #penalty matrix for function
  in.mat1 = inprod(bbasis1, bbasis1) #penalty matrix for function
  Kphi0 <- in.mat0
  basismat0 <- bbasisMat0
  for (i in 1:n) {
    colnames(basismat0[[i]]) <- NULL
  }
  Kphi1 <- in.mat1
  basismat1 <- bbasisMat1
  for (i in 1:n) {
    colnames(basismat1[[i]]) <- NULL
  }
  
  XMatstar <- list()
  for (i in 1:n)
  {
    XMatstar[[i]] <- array(0, dim = c(1, length(T[[i]]), nbasis1))
    for (j in 1:1) {
      for (l in 1:length(T[[i]])) {
        XMatstar[[i]][j, l, ] = Xhat_resp[[i]][l] * basismat1[[i]][l, ]
      }
    }
  }
  B <- basismat0
  R0 <- chol(Kphi0)
  L0 <- t(R0)
  M0 <- solve(R0)
  R1 <- chol(Kphi1)
  L1 <- t(R1)
  M1 <- solve(R1)
  ZMatstar <- list()
  for (i in 1:n)
  {
    ZMatstar[[i]] <- array(0, dim = c(1, length(T[[i]]), nbasis1))
    for (j in 1:1) {
      for (l in 1:length(T[[i]])) {
        ZMatstar[[i]][j, l, ] = XMatstar[[i]][j, l, ] %*% M1
      }
    }
  }
  Zimatlist1 <- list()
  for (i in 1:n)
  {
    Zimatlist1[[i]] <- matrix(0, length(T[[i]]), nbasis1 * 1)
    for (l in 1:length(T[[i]]))
    {
      for (j in 1:1) {
        Zimatlist1[[i]][l, ((nbasis1 * (j - 1)) + 1):(nbasis1 * j)] <-
          ZMatstar[[i]][j, l, ]
      }
    }
  }
  #B<-basismat
  for (i in 1:n) {
    B[[i]] <- B[[i]] %*% M0
  }
  Grandz1 <- NULL
  for (i in 1:n)
  {
    Grandz1 <- rbind(Grandz1, Zimatlist1[[i]])
  }
  GrandB <- NULL
  for (i in 1:n)
  {
    GrandB <- rbind(GrandB, B[[i]])
  }
  matlist <- list(GrandY, GrandB, Grandz1, T, n)
  names(matlist) <- c("GrandY", "GrandB", "Grandz1", "T", "n")
  return(matlist)
}

#####Main Function For testing: Requires Matrices in stacked form#####
scorefuncpvalSparse <- function(GrandY, GrandB, Grandz1, T, n) {
  Rcpp::sourceCpp('test.cpp')
  nbasis0 <- ncol(GrandB)
  nbasis1 <- ncol(Grandz1)
  library(mgcv)
  PP <- list(GrandB = list(diag(nbasis0)), Grandz1 = list(diag(nbasis1)))
  fit <- bam(GrandY ~ -1 + GrandB + Grandz1,
             paraPen = PP,
             method = "REML")
  res <- fit$residuals
  lenT <- c()
  for (i in 1:n)
  {
    lenT[i] <- length(T[[i]])
  }
  .id = rep(1:n, lenT)
  .index <- unlist(T)
  resdata <- data.frame(.id , .index)
  resdata$.value <- res
  library(refund)
  fpca1 <-
    fpca.sc(
      ydata = resdata,
      center = TRUE,
      nbasis = 15,
      pve = 0.99,
      var = TRUE
    ) #tune pve if needed
  timeU <- as.numeric(attr(fpca1$Yhat, "dimnames")[[2]])
  m <- length(timeU)
  if (length(fpca1$evalues) > 1)
  {
    sigmat <-
      fpca1$efunctions %*% diag(fpca1$evalues) %*% t(fpca1$efunctions) + fpca1$sigma2 *
      diag(m)
  }
  if (length(fpca1$evalues) == 1)
  {
    sigmat <-
      fpca1$evalues * fpca1$efunctions %*% t(fpca1$efunctions) + fpca1$sigma2 *
      diag(m)
  }
  indtime <- list()
  for (i in 1:n)
  {
    indtime[[i]] <- vector(mode = "logical", length = length(T[[i]]))
    for (l in 1:length(T[[i]])) {
      indtime[[i]][l] = which(round(timeU, 8) == round(T[[i]][l], 8))
    }
  }
  subsigmat <- list()
  for (i in 1:n)
  {
    subsigmat[[i]] <- sigmat[indtime[[i]], indtime[[i]]]
  }
  
  Linv <- list()
  for (i in 1:n)
  {
    R <- chol(subsigmat[[i]])
    L <- t(R)
    Linv[[i]] <- forwardsolve(L, diag(dim(L)[1]))
  }
  Linvhalf <- bdiag(Linv)
  NewB <- as.matrix(Linvhalf %*% GrandB)
  Newz1 <- as.matrix(Linvhalf %*% Grandz1)
  NewY <- as.matrix(Linvhalf %*% GrandY)
  Rcpp::sourceCpp('test.cpp')
  nbasis0 <- ncol(NewB)
  nbasis1 <- ncol(Newz1)
  PPreduced <- list(NewB = list(diag(nbasis0)))
  fitfinalreduced <-
    bam(NewY ~ -1 + NewB, paraPen = PPreduced, method = "REML")
  mvar <-
    (1 / fitfinalreduced$sp)      #gam.vcomp(fitfinalreduced) #use 1/sp
  #w<-m^2
  estt0 <- mvar[[1]]
  Wc <- sqrt(estt0) * NewB
  meat <- diag(ncol(Wc)) + eigenMapMatMult(t(Wc), Wc)
  meatinv <- chol2inv(chol(meat))
  Yvi <-
    NewY - eigenMapMatMult(eigenMapMatMult(eigenMapMatMult(Wc, meatinv), t(Wc)), NewY)
  Eigenmat <-
    eigenMapMatMult(t(Newz1), Newz1) - eigenMapMatMult(eigenMapMatMult(t(Newz1), eigenMapMatMult(eigenMapMatMult(Wc, meatinv), t(Wc))), Newz1)
  lambda <- eigen(Eigenmat)$values
  s1l <- sum(lambda)
  s2l <- .5 * sum((lambda) ^ 2)
  term1 <- eigenMapMatMult(t(Newz1), Yvi)
  Score <- -.5 * (s1l - eigenMapMatMult(t(term1), term1))
  EigenmatB <-
    eigenMapMatMult(t(NewB), NewB) - eigenMapMatMult(eigenMapMatMult(t(NewB), eigenMapMatMult(eigenMapMatMult(Wc, meatinv), t(Wc))), NewB)
  lambdaB <- eigen(EigenmatB)$values
  s2lB <- .5 * sum((lambdaB) ^ 2)
  EigenmatBz1 <-
    eigenMapMatMult(t(Newz1), NewB) - eigenMapMatMult(eigenMapMatMult(t(Newz1), eigenMapMatMult(eigenMapMatMult(Wc, meatinv), t(Wc))), NewB)
  trac <- function(A) {
    sum(diag(A))
  }
  s2lBz1 <- .5 * trac(eigenMapMatMult(EigenmatBz1, t(EigenmatBz1)))
  Inform <- matrix(0, 2, 2)
  Inform[1, 1] <- s2lB
  Inform[2, 2] <- s2l
  Inform[1, 2] <- s2lBz1
  Inform[2, 1] <- s2lBz1
  I <- Inform[2, 2] - t(Inform[1, 2]) %*% solve(Inform[1, 1]) %*% Inform[1, 2]
  if (Score > 0)
  {
    Tsobs <- (t(Score) %*% Score) / I
  }
  if (Score <= 0)
  {
    Tsobs = 0
  }
  Tobs <- as.numeric(Tsobs)
  simuT <- function()
  {
    w <- rnorm(nbasis1)
    s <- sum(lambda * w ^ 2) - s1l
    T1 <- (1 / 4) * (s) ^ 2
    if (s > 0)
    {
      return((T1 / I))
    }
    if (s <= 0)
    {
      return(0)
    }
  }
  Tsnull <- replicate(100000, simuT())
  pval <- (sum(Tsnull >= Tobs)) / 100000
  pvallist <- list(Tobs, Score, pval)
  names(pvallist) <- c("Tobs", "Score", "pval")
  return(pvallist)
}

### The wrapper function testpvalSparse which uses the above two functions 
### Reprots P-Value of the test
FLCM.test.sparse1 <- function(y, x1, T, n, ID, nbas = 7)
{
  Matlist <- GetMatricesSparse(y, x1, T, n, ID, nbas)
  GrandY <- Matlist$GrandY
  GrandB <- Matlist$GrandB
  Grandz1 <- Matlist$Grandz1
  T <- Matlist$T
  n <- Matlist$n
  scoremethodresult <- scorefuncpvalSparse(GrandY, GrandB, Grandz1, T, n)
  scorepval <- scoremethodresult$pval
  scorepval
}

#####################################################################################
##########test functions for sparse data,testing for x2 in presence of x1 ###########
#Function for getting stacked matrices
GetMatricesSparse2 <- function(y, x1, x2, T, n, ID, nbas)
{
  UT <- c()
  for (i in 1:n)
  {
    UT <- union(UT, T[[i]])
  }
  UT <- sort(UT)
  Y <- list()
  for (i in 1:n)
  {
    ind <- which(ID == i)
    Y[[i]] <- y[ind]
  }
  GrandY <- NULL
  for (i in 1:n)
  {
    GrandY <- c(GrandY, Y[[i]])
  }
  W <- list()
  for (i in 1:n)
  {
    ind <- which(ID == i)
    W[[i]] <- x1[ind]  
  }
  lenX <- c()
  for (i in 1:n)
  {
    lenX[i] <- length(T[[i]])
  }
  .id = rep(1:n, lenX)
  .index <- unlist(T)
  Wdata <- data.frame(.id , .index)
  Wdata$.value <- unlist(W)
  #Using FPCA to preprocess noisy data
  w.sm = fpca.sc(ydata = Wdata, pve = .99, var = TRUE)
  xhat = w.sm$Yhat
  timepoints <- as.numeric(attr(w.sm$Yhat, "dimnames")[[2]])
  Xhat_resp <- list()
  for (i in 1:n)
  {
    Xhat_resp[[i]] <- vector(mode = "logical", length = length(T[[i]]))
    for (l in 1:length(T[[i]])) {
      ind <- which(round(timepoints, 8) == round(T[[i]], 8)[l])
      Xhat_resp[[i]][l] = xhat[i, ind]
    }
  }
  w2 <- list()
  for (i in 1:n)
  {
    ind <- which(ID == i)
    w2[[i]] <- x2[ind]
  }
  .id = rep(1:n, lenX)
  .index <- unlist(T)
  Wdata2 <- data.frame(.id , .index)
  Wdata2$.value <- unlist(w2)
  #Using FPCA to preprocess noisy data
  w.sm2 = fpca.sc(ydata = Wdata2, pve = .99, var = TRUE)
  xhat2 = w.sm2$Yhat
  timepoints <- as.numeric(attr(w.sm$Yhat, "dimnames")[[2]])
  Xhat2_resp <- list()
  for (i in 1:n)
  {
    Xhat2_resp[[i]] <- vector(mode = "logical", length = length(T[[i]]))
    for (l in 1:length(T[[i]])) {
      ind <- which(round(timepoints, 8) == round(T[[i]], 8)[l])
      Xhat2_resp[[i]][l] = xhat2[i, ind]
    }
  }
  #equispaced knots used, one can used alternate knots on quantiles or other config based on domain knowldge
  library(fda)
  knots0 = seq(min(UT), max(UT), l = 5)
  knots1 = seq(min(UT), max(UT), l = 5)
  knots2 = seq(min(UT), max(UT), l = nbas - 2)
  # We'll use fourth-order B-splines
  norder = 4
  # this implies the number of basis functions
  nbasis0 = length(knots0) + norder - 2
  nbasis1 = length(knots1) + norder - 2
  nbasis2 = length(knots2) + norder - 2
  dayrng = c(0, 1)
  bbasis0 = create.bspline.basis(dayrng, nbasis0, norder, knots0)
  bbasis1 = create.bspline.basis(dayrng, nbasis1, norder, knots1)
  bbasis2 = create.bspline.basis(dayrng, nbasis2, norder, knots2)
  #bbasis$nbasis    # number of basis functions
  #bbasis$rangeval   # basis range
  bbasisMat0 <- list()
  for (i in 1:n)
  {
    bbasisMat0[[i]] = eval.basis(T[[i]], bbasis0)
  }
  bbasisMat1 <- list()
  for (i in 1:n)
  {
    bbasisMat1[[i]] = eval.basis(T[[i]], bbasis1)
  }
  bbasisMat2 <- list()
  for (i in 1:n)
  {
    bbasisMat2[[i]] = eval.basis(T[[i]], bbasis2)
  }
  
  #dim(bbasisMat)
  #100 102
  in.mat0 = inprod(bbasis0, bbasis0) #penalty matrix for function
  in.mat1 = inprod(bbasis1, bbasis1) #penalty matrix for function
  in.mat2 = inprod(bbasis2, bbasis2) #penalty matrix for function
  
  #par(mfrow=c(1,1))
  #image(in.mat)
  Kphi0 <- in.mat0
  basismat0 <- bbasisMat0
  for (i in 1:n) {
    colnames(basismat0[[i]]) <- NULL
  }
  Kphi1 <- in.mat1
  basismat1 <- bbasisMat1
  for (i in 1:n) {
    colnames(basismat1[[i]]) <- NULL
  }
  Kphi2 <- in.mat2
  basismat2 <- bbasisMat2
  for (i in 1:n) {
    colnames(basismat2[[i]]) <- NULL
  }
  
  XMatstar <- list()
  for (i in 1:n)
  {
    XMatstar[[i]] <- array(0, dim = c(1, length(T[[i]]), nbasis1))
    for (j in 1:1) {
      for (l in 1:length(T[[i]])) {
        XMatstar[[i]][j, l, ] = Xhat_resp[[i]][l] * basismat1[[i]][l, ]
      }
    }
  }
  XMatstar2 <- list()
  for (i in 1:n)
  {
    XMatstar2[[i]] <- array(0, dim = c(1, length(T[[i]]), nbasis2))
    for (j in 1:1) {
      for (l in 1:length(T[[i]])) {
        XMatstar2[[i]][j, l, ] = Xhat2_resp[[i]][l] * basismat2[[i]][l, ]
      }
    }
  }
  B <- basismat0
  R0 <- chol(Kphi0)
  L0 <- t(R0)
  M0 <- solve(R0)
  R1 <- chol(Kphi1)
  L1 <- t(R1)
  M1 <- solve(R1)
  R2 <- chol(Kphi2)
  L2 <- t(R2)
  M2 <- solve(R2)
  ZMatstar <- list()
  for (i in 1:n)
  {
    ZMatstar[[i]] <- array(0, dim = c(1, length(T[[i]]), nbasis1))
    for (j in 1:1) {
      for (l in 1:length(T[[i]])) {
        ZMatstar[[i]][j, l, ] = XMatstar[[i]][j, l, ] %*% M1
      }
    }
  }
  
  ZMatstar2 <- list()
  for (i in 1:n)
  {
    ZMatstar2[[i]] <- array(0, dim = c(1, length(T[[i]]), nbasis2))
    for (j in 1:1) {
      for (l in 1:length(T[[i]])) {
        ZMatstar2[[i]][j, l, ] = XMatstar2[[i]][j, l, ] %*% M2
      }
    }
  }
  Zimatlist1 <- list()
  for (i in 1:n)
  {
    Zimatlist1[[i]] <- matrix(0, length(T[[i]]), nbasis1 * 1)
    for (l in 1:length(T[[i]]))
    {
      for (j in 1:1) {
        Zimatlist1[[i]][l, ((nbasis1 * (j - 1)) + 1):(nbasis1 * j)] <-
          ZMatstar[[i]][j, l, ]
      }
    }
  }
  
  Zimatlist2 <- list()
  for (i in 1:n)
  {
    Zimatlist2[[i]] <- matrix(0, length(T[[i]]), nbasis2 * 1)
    for (l in 1:length(T[[i]]))
    {
      for (j in 1:1) {
        Zimatlist2[[i]][l, ((nbasis2 * (j - 1)) + 1):(nbasis2 * j)] <-
          ZMatstar2[[i]][j, l, ]
      }
    }
  }
  #B<-basismat
  for (i in 1:n) {
    B[[i]] <- B[[i]] %*% M0
  }
  Grandz1 <- NULL
  for (i in 1:n)
  {
    Grandz1 <- rbind(Grandz1, Zimatlist1[[i]])
  }
  Grandz2 <- NULL
  for (i in 1:n)
  {
    Grandz2 <- rbind(Grandz2, Zimatlist2[[i]])
  }
  GrandB <- NULL
  for (i in 1:n)
  {
    GrandB <- rbind(GrandB, B[[i]])
  }
  matlist <- list(GrandY, GrandB, Grandz1, Grandz2, T, n)
  names(matlist) <- c("GrandY", "GrandB", "Grandz1", "Grandz2", "T", "n")
  return(matlist)
}
#####Main Function For testing: Requires Matrices in stacked form#####
scorefuncpvalSparse2 <- function(GrandY, GrandB, Grandz1, Grandz2, T, n) {
  Rcpp::sourceCpp('test.cpp')
  nbasis0<-ncol(GrandB)
  nbasis1<-ncol(Grandz1)
  nbasis2<-ncol(Grandz2)
  library(mgcv)
  PP <-list(GrandB=list(diag(nbasis0)),Grandz1=list(diag(nbasis1)),Grandz2=list(diag(nbasis2)))
  fit <-bam(GrandY~ -1+GrandB+Grandz1+Grandz2,paraPen=PP,method="REML")
  res<-fit$residuals
  lenT<-c()
  for (i in 1:n)
  {lenT[i]<-length(T[[i]])}
  .id=rep(1:n,lenT)
  .index<-unlist(T)
  resdata <- data.frame(.id ,.index)
  resdata$.value<-res
  library(refund)
  fpca1<-fpca.sc(ydata=resdata,center = TRUE, nbasis=15,pve=0.99,var=TRUE) #Tune pve if needed
  timeU<-as.numeric(attr(fpca1$Yhat,"dimnames")[[2]])
  m<-length(timeU)
  if(length(fpca1$evalues)>1)
  {sigmat<-fpca1$efunctions%*%diag(fpca1$evalues)%*%t(fpca1$efunctions) + fpca1$sigma2*diag(m)}
  if(length(fpca1$evalues)==1)
  {sigmat<-fpca1$evalues*fpca1$efunctions%*%t(fpca1$efunctions) + fpca1$sigma2*diag(m)}
  indtime<-list()
  for(i in 1:n)
  {indtime[[i]]<-vector(mode = "logical",length = length(T[[i]]))
  for(l in 1: length(T[[i]])){
    indtime[[i]][l]=which(round(timeU,8)==round(T[[i]][l],8))
  }
  }
  subsigmat<-list()
  for(i in 1:n)
  {subsigmat[[i]]<-sigmat[indtime[[i]],indtime[[i]]]
  }
  
  Linv<-list()
  for(i in 1:n)
  {R<-chol(subsigmat[[i]])
  L<-t(R)
  Linv[[i]]<-forwardsolve(L, diag(dim(L)[1]))}
  Linvhalf<-bdiag(Linv)
  NewB<-as.matrix(Linvhalf%*%GrandB)
  Newz1<-as.matrix(Linvhalf%*%Grandz1)
  Newz2<-as.matrix(Linvhalf%*%Grandz2)
  NewY<-as.matrix(Linvhalf%*%GrandY)
  Rcpp::sourceCpp('test.cpp')
  nbasis0<-ncol(NewB)
  nbasis1<-ncol(Newz1)
  nbasis2<-ncol(Newz2)
  PPreduced <-list(NewB=list(diag(nbasis0)),Newz1=list(diag(nbasis1)))
  fitfinalreduced <-bam(NewY~ -1+NewB+Newz1,paraPen=PPreduced,method="REML")
  mvar<-(1/fitfinalreduced$sp)      #gam.vcomp(fitfinalreduced) #use 1/sp
  #w<-m^2
  estt0<-mvar[[1]]
  estt1<-mvar[[2]]
  Wc<-cbind(sqrt(estt0)*NewB,sqrt(estt1)*Newz1)
  meat<-diag(ncol(Wc))+eigenMapMatMult(t(Wc),Wc)
  meatinv<-chol2inv(chol(meat))
  Yvi<-NewY-eigenMapMatMult(eigenMapMatMult(eigenMapMatMult(Wc,meatinv),t(Wc)),NewY)
  Eigenmat<-eigenMapMatMult(t(Newz2),Newz2)-eigenMapMatMult(eigenMapMatMult(t(Newz2),eigenMapMatMult(eigenMapMatMult(Wc,meatinv),t(Wc))),Newz2)
  lambda<-eigen(Eigenmat)$values
  s1l<-sum(lambda)
  s2l<-.5*sum((lambda)^2)
  term1<-eigenMapMatMult(t(Newz2),Yvi)
  Score<--.5*(s1l-eigenMapMatMult(t(term1),term1))
  EigenmatB<-eigenMapMatMult(t(NewB),NewB)-eigenMapMatMult(eigenMapMatMult(t(NewB),eigenMapMatMult(eigenMapMatMult(Wc,meatinv),t(Wc))),NewB)
  lambdaB<-eigen(EigenmatB)$values
  s2lB<-.5*sum((lambdaB)^2)
  Eigenmatz1<-eigenMapMatMult(t(Newz1),Newz1)-eigenMapMatMult(eigenMapMatMult(t(Newz1),eigenMapMatMult(eigenMapMatMult(Wc,meatinv),t(Wc))),Newz1)
  lambdaz1<-eigen(Eigenmatz1)$values
  s2lz1<-.5*sum((lambdaz1)^2)
  EigenmatBz1<-eigenMapMatMult(t(Newz1),NewB)-eigenMapMatMult(eigenMapMatMult(t(Newz1),eigenMapMatMult(eigenMapMatMult(Wc,meatinv),t(Wc))),NewB)
  EigenmatBz2<-eigenMapMatMult(t(Newz2),NewB)-eigenMapMatMult(eigenMapMatMult(t(Newz2),eigenMapMatMult(eigenMapMatMult(Wc,meatinv),t(Wc))),NewB)
  Eigenmatz2z1<-eigenMapMatMult(t(Newz1),Newz2)-eigenMapMatMult(eigenMapMatMult(t(Newz1),eigenMapMatMult(eigenMapMatMult(Wc,meatinv),t(Wc))),Newz2)
  trac<-function(A){sum(diag(A))}
  s2lBz1<-.5*trac(eigenMapMatMult(EigenmatBz1,t(EigenmatBz1)))
  s2lBz2<-.5*trac(eigenMapMatMult(EigenmatBz2,t(EigenmatBz2)))
  s2lz2z1<-.5*trac(eigenMapMatMult(Eigenmatz2z1,t(Eigenmatz2z1)))
  Inform<-matrix(0,3,3)
  Inform[1,1]<-s2lB
  Inform[2,2]<-s2lz1
  Inform[3,3]<-s2l
  Inform[1,2]<-s2lBz1
  Inform[2,1]<-s2lBz1
  Inform[1,3]<-s2lBz2
  Inform[3,1]<-s2lBz2
  Inform[2,3]<-s2lz2z1
  Inform[3,2]<-s2lz2z1
  I<-Inform[3,3]-t(Inform[1:2,3])%*%solve(Inform[1:2,1:2])%*%Inform[1:2,3]
  if(Score > 0)
  {
    Tsobs<-(t(Score)%*%Score)/I
  }
  if(Score<=0)
  {Tsobs=0
  }
  Tobs<-as.numeric(Tsobs)
  simuT<-function()
  {
    w<-rnorm(nbasis2)
    s<-sum(lambda*w^2) -s1l
    T1<-(1/4)*(s)^2
    if(s>0)
    {return((T1/I))}
    if(s<=0)
    {return(0)}
  }
  Tsnull<-replicate(100000,simuT())
  pval<-(sum(Tsnull>=Tobs))/100000
  pvallist<-list(Tobs,Score,pval)
  names(pvallist)<-c("Tobs","Score","pval")
  return(pvallist)
}
### The wrapper function testpvalSparse which uses the above two functions 
### Reprots P-Value of the test
FLCM.test.sparse2 <- function(y, x1, x2, T, n, ID, nbas = 7)
{
  Matlist <- GetMatricesSparse2(y, x1, x2, T, n, ID, nbas)
  GrandY <- Matlist$GrandY
  GrandB <- Matlist$GrandB
  Grandz1 <- Matlist$Grandz1
  Grandz2 <- Matlist$Grandz2
  T <- Matlist$T
  n <- Matlist$n
  scoremethodresult <- scorefuncpvalSparse2(GrandY, GrandB, Grandz1, Grandz2, T, n)
  scorepval <- scoremethodresult$pval
  scorepval
}
