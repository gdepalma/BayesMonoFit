

getLLKManyKnots=function(yobs,y_mu,sigma2){
  
  yllk <- dnorm(yobs,y_mu,sqrt(sigma2),log=TRUE)
  
  llike <- sum(yllk)
  
  return(llike)
}

PriorLLKManyKnots=function(coefs,smoothParam,sigma2){
  
  priorSmooth=dunif(smoothParam,0,2,log=T)
  

  lpriorcoef = dlnorm(coefs[1],0,100,log=T)
  for(i in 2:length(coefs))
    lpriorcoef=lpriorcoef+dnorm(log(coefs[i]),log(coefs[i-1]),smoothParam,log=T)
  
  
  lpriorerr=log(1/dgamma(sigma2,.01,.01))
  
  return(priorSmooth+lpriorcoef+lpriorerr)
}


updateCoefsManyKnots=function(smoothParam,coefs,iter,bases,y_mu,knotseq,xobs,yobs,sigma2,coefMat){

  if(iter<1000){
    Sigma=diag(.001,nrow=ncol(coefMat),ncol=ncol(coefMat))
  }else{
    Sigma=cov(coefMat[(iter-900):(iter-1),])
  }

  propCoef=tryCatch(as.numeric(exp(rmvnorm(n=1,as.numeric(log(coefs)),sigma=Sigma))),
          error=function(e){as.numeric(exp(rmvnorm(n=1,as.numeric(log(coefs)),
          sigma=diag(.001,nrow=nrow(bases),ncol=nrow(bases)))))})

  accept=0
  designMatrix=getIsplineC(xobs,knotseq,bases)
  newY=as.numeric(propCoef%*%t(designMatrix))

  #log likelihood and priors
  oldLLK=getLLKManyKnots(yobs,y_mu,sigma2)
  oldPriorLLKManyKnots=PriorLLKManyKnots(coefs,smoothParam,sigma2)
  newLLK=getLLKManyKnots(yobs,newY,sigma2)
  newPriorLLKManyKnots=PriorLLKManyKnots(propCoef,smoothParam,sigma2)

  #accept/reject
  A=newLLK+newPriorLLKManyKnots-oldLLK-oldPriorLLKManyKnots
  if(is.nan(A)==FALSE & log(runif(1)) < A){
    y_mu=newY
    coefs=propCoef
    accept=1
  }

  return(list(acceptCoef=accept,coefs=coefs,y_mu=y_mu))

}

updateSmoothParmManyKnots=function(smoothParam,coefs,sigma2){

  accept=0
  propSmooth=rnorm(1,smoothParam,.2)
  if(propSmooth<0 | propSmooth>2)
    return(list(smooth=smoothParam,accept=accept))
  oldPriorLLKManyKnots=PriorLLKManyKnots(coefs,smoothParam,sigma2)
  newPriorLLKManyKnots=PriorLLKManyKnots(coefs,propSmooth,sigma2)
  #accept/reject
  if(log(runif(1)) < newPriorLLKManyKnots-oldPriorLLKManyKnots){
    smoothParam=propSmooth
    accept=1
  }

  return(list(smooth=smoothParam,accept=accept))
}

updateVarianceManyKnots=function(yobs,y_mu,sigma2,smoothParam,coefs){
  
  accept=0
  propSigma2=rnorm(1,sigma2,.1)
  if(propSigma2<0)
    return(sigma2)
  #log likelihood and priors
  oldLLK=getLLKManyKnots(yobs,y_mu,sigma2)
  oldPriorLLKManyKnots=PriorLLKManyKnots(coefs,smoothParam,sigma2)
  newLLK=getLLKManyKnots(yobs,y_mu,propSigma2)
  newPriorLLKManyKnots=PriorLLKManyKnots(coefs,smoothParam,propSigma2)
                  
  #accept/reject
  A=newLLK+newPriorLLKManyKnots-oldLLK-oldPriorLLKManyKnots
  if(is.nan(A)==FALSE & log(runif(1)) < A){
    sigma2=propSigma2
  }
  
  return(sigma2)
}
