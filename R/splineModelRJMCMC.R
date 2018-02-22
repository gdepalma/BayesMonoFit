
runRJMCMC=function(initializeValues,numIter=9000,burnin=5000,thin=4){

  if(class(initializeValues)!="initialRJMCMC") stop("Wrong input data type")
  coefs=initializeValues$coefs; bases=initializeValues$bases; knotseq=initializeValues$knotseq
  y_mu=initializeValues$y_mu; xgrid=initializeValues$xgrid; designMatrix=initializeValues$designMatrix
  designMatrixGrid=initializeValues$designMatrixGrid; xobs=initializeValues$xobs; yobs=initializeValues$yobs
  
  lowerept=knotseq[1]
  upperept=knotseq[length(knotseq)]
  
  fitMat=matrix(nrow=numIter-burnin,ncol=length(xgrid))
  coefMat=matrix(nrow=numIter,ncol=length(coefs))
  sigma2_save=rep(NA,numIter)
  acceptCoef=rep(NA,numIter)
  
  numIter=9000; burnin=5000
  
  nobs=length(xobs)

  sigma2=5
  

  for(iter in 1:numIter){

    if(iter %% 500==0){
      cat('iteration: ',iter,'\n')
    }
    
    ### Update Coefs
    parms=updateCoefsRJMCMC(coefs,bases,knotseq,xobs,yobs,y_mu,sigma2)
    coefs=parms$coefs; y_mu=as.numeric(parms$y_mu); acceptCoef[iter]=parms$accept

    
    ### Update error variance
    sigma2=updateVarianceRJMCMC(yobs,y_mu,sigma2,coefs)
    
    ### RJMCMC interior knots
    

    coefMat[iter,]=coefs
    sigma2_save[iter]=sigma2
    if(iter>burnin){
      fitMat[iter-burnin,]=coefs%*%t(designMatrixGrid)
    }

  }
  
  plot(1:iter,sigma2_save,type='l')
  plot(1:iter,acceptCoef)
  mean(acceptCoef)
  
  resultData=list(fitMat=fitMat,coefMat=coefMat,sigma2_save=sigma2_save,
                        xobs=xobs,yobs=yobs,xgrid=xgrid,burnin=burnin)
  class(resultData)="RJMCMCResults"
  
  return(resultData)

}







