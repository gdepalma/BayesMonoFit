

getLLKRJMCMC=function(yobs,y_mu,sigma2){
  
  yllk <- dnorm(yobs,y_mu,sqrt(sigma2),log=TRUE)
  
  llike <- sum(yllk)
  
  return(llike)
}

PriorLLKRJMCMC=function(coefs,sigma2){
  
  lam=4
  numInteriorKnots=length(coefs)-3
  
  lpriorNumKnots=dpois(numInteriorKnots,lam,log=TRUE)
  lpriorcoef=dnorm(coefs,0,100,log=TRUE)
  lpriorerr=log(1/dgamma(sigma2,.01,.01))
  
  return(lpriorNumKnots+sum(lpriorcoef)+lpriorerr)
}


updateCoefsRJMCMC=function(coefs,bases,knotseq,xobs,yobs,y_mu,sigma2){

  #get covariance
  designMatrix=getIsplineC(xobs,knotseq,bases)
  fit=lm(yobs~designMatrix-1)
  Sigma=vcov(fit)
  
  #### joint adaptive update
  propCoef=as.numeric(rmvnorm(n=1,coefs,sigma=Sigma))
  newY=as.numeric(propCoef%*%t(designMatrix))

  accept=0

  #log likelihood and priors
  oldLLK=getLLKRJMCMC(yobs,y_mu,sigma2)
  oldPRIORLLK=PriorLLKRJMCMC(coefs,sigma2)
  newLLK=getLLKRJMCMC(yobs,newY,sigma2)
  newPRIORLLK=PriorLLKRJMCMC(propCoef,sigma2)

  #accept/reject
  A=newLLK+newPRIORLLK-oldLLK-oldPRIORLLK
  if(is.nan(A)==FALSE & log(runif(1)) < A){
    y_mu=newY
    coefs=propCoef
    accept=1
  }

  return(list(acceptCoef=accept,coefs=coefs,y_mu=y_mu))

}


updateVarianceRJMCMC=function(yobs,y_mu,sigma2,coefs){
  
  accept=0
  propSigma2=rnorm(1,sigma2,.1)
  if(propSigma2<0)
    return(sigma2)
  #log likelihood and priors
  oldLLK=getLLKRJMCMC(yobs,y_mu,sigma2)
  oldPRIORLLK=PriorLLKRJMCMC(coefs,sigma2)
  newLLK=getLLKRJMCMC(yobs,y_mu,propSigma2)
  newPRIORLLK=PriorLLKRJMCMC(coefs,propSigma2)
                  
  #accept/reject
  A=newLLK+newPRIORLLK-oldLLK-oldPRIORLLK
  if(is.nan(A)==FALSE & log(runif(1)) < A){
    sigma2=propSigma2
  }
  
  return(sigma2)
}

updateKnots=function(knotseq,lowerept,upperept,bases,coefs,xobs,yobs,y_mu,sigma2,designMatrix,designMatrixGrid,xgrid){
  
  numInteriorKnots=length(coefs)-3
  
  choice=ifelse(numInteriorKnots==1,sample(2,1),sample(3,1))
  if(choice==1) parms=birthKnot(knotseq,lowerept,upperept,bases,coefs,xobs,yobs,y_mu,sigma2,designMatrix,designMatrixGrid,xgrid)
  if(choice==2) parms=moveKnot(knotseq,lowerept,upperept,bases,coefs,xobs,yobs,y_mu,sigma2,designMatrix,designMatrixGrid,xgrid)
  if(choice==3) parms=deathKnot(knotseq,lowerept,upperept,bases,coefs,xobs,yobs,y_mu,sigma2,designMatrix,designMatrixGrid,xgrid)

  
  # as.list( sys.call() )
}

birthKnot=function(knotseq,lowerept,upperept,bases,coefs,xobs,yobs,y_mu,sigma2,designMatrix,designMatrixGrid,xgrid){

  numInteriorKnots=length(coefs)-3
  interiorKnots=knotseq[5:(5+numInteriorKnots-1)]
  
  # interval for new knot
  interval=sample(numInteriorKnots+1,1)
  lowerKnot=knotseq[4+interval-1]
  upperKnot=knotseq[4+interval]
  newKnot=runif(1,lowerKnot,upperKnot)
  interval_diff=newKnot-lowerKnot

  # new knot sequence
  parms=Ispline(sort(c(lowerKnot,upperKnot,newKnot)),lowerept,upperept)
  propBases=parms$bases;
  propKnotseq=parms$knotseq

  
  # old LS coefficients
  fit=lm(yobs~designMatrix-1)
  coefOldLS=as.numeric(coef(fit))
  sigmaCoefOld=vcov(fit)
  oldY=as.numeric(coefs%*%t(designMatrix))
  
  
  # new LS coefficients
  designMatrix=getIsplineC(xobs,propKnotseq,propBases)
  fit=lm(yobs~designMatrix-1)
  coefNewLS=as.numeric(coef(fit))
  sigmaCoefNew=vcov(fit)

  # new proposal coefficients
  coefNew=tryCatch(as.numeric(rmvnorm(n=1,coefNewLS,sigma=sigmaCoefNew)),
                   error=function(e){print("Error in birth new coefficients")})
  newY=as.numeric(coefNew%*%t(designMatrix))
  
  # test monotonicity
  tempY=newY 
  DT=data.table(xobs,tempY)
  DT=DT[order(xobs)]
  test_mono_y=DT[,2]
  if(all(test_mono_y == cummin(test_mono_y))==FALSE) print("Birth not monotonic")  

  
  #likelihoods and priors
  newLLK=getLLKRJMCMC(newY,y_mu,sigma2)
  oldLLK=getLLKRJMCMC(yobs,y_mu,sigma2)
  newPRIORLLK=PriorLLKRJMCMC(coefNew,sigma2)
  oldPRIORLLK=PriorLLKRJMCMC(coefs,sigma2)
  
  
  # multivariate Normal likelihoods
  oldMV=dmvnorm(coefs,coefOldLS,sigmaCoefOld,log=TRUE)
  newMV=dmvnorm(coefNew,coefNewLS,sigmaCoefNew,log=TRUE)
  
  
  num=newLLK+newPRIORLLK+newMV
  dem=oldLLK+oldPRIORLLK+oldMV+log(1/interval_diff)
  
  A=num-dem
  
  if(is.na(A)==TRUE) A=-99
  
  #accept/reject
  if(log(runif(1)) < A){
    return(list(knotseq=propKnotseq,bases=propBases,accept=1,ytrue=newY,coef=coefNew))
  }else{
    return(list(knotseq=knotseq,bases=bases,accept=0,ytrue=ytrue1,coef=icoefs1))
  }
  
}

moveKnot=function(knotseq,lowerept,upperept,bases,coefs,xobs,yobs,y_mu,sigma2,designMatrix,designMatrixGrid,xgrid){
  
  
}

deathKnot=function(knotseq,lowerept,upperept,bases,coefs,xobs,yobs,y_mu,sigma2,designMatrix,designMatrixGrid,xgrid){
  
  
}
