
runManyKnots=function(initializeValues,numIter=9000,burnin=5000,thin=4){

  if(class(initializeValues)!="initialManyKnots") stop("Wrong input data type")
  coefs=initializeValues$coefs; bases=initializeValues$bases; knotseq=initializeValues$knotseq
  y_mu=initializeValues$y_mu; xgrid=initializeValues$xgrid; designMatrix=initializeValues$designMatrix
  designMatrixGrid=initializeValues$designMatrixGrid; xobs=initializeValues$xobs; yobs=initializeValues$yobs
  
  
  fitMat=matrix(nrow=numIter-burnin,ncol=length(xgrid))
  coefMat=matrix(nrow=numIter,ncol=length(coefs))
  sigma2_save=rep(NA,numIter)
  smoothParam_save=rep(NA,numIter)
  acceptCoef=rep(NA,numIter)
  
  
  nobs=length(xobs)

  smoothParam=2
  sigma2=5
  

  for(iter in 1:numIter){

    if(iter %% 500==0){
      cat('iteration: ',iter,'\n')
    }
    
    ### Update Coefs
    parms=updateCoefsManyKnots(smoothParam,coefs,iter,bases,y_mu,knotseq,xobs,yobs,sigma2,coefMat)
    coefs=parms$coefs; y_mu=as.numeric(parms$y_mu); acceptCoef[iter]=parms$accept


    ### Update smoothness parameter
    parms=updateSmoothParmManyKnots(smoothParam,coefs,sigma2)
    smoothParam=parms$smooth
    
    ### Update error variance
    sigma2=updateVarianceManyKnots(yobs,y_mu,sigma2,smoothParam,coefs)

    coefMat[iter,]=log(coefs)
    sigma2_save[iter]=sigma2
    smoothParam_save[iter]=smoothParam
    if(iter>burnin){
      fitMat[iter-burnin,]=coefs%*%t(designMatrixGrid)
    }

  }
  
  resultData=list(fitMat=fitMat,coefMat=coefMat,sigma2_save=sigma2_save,smoothParam_save=smoothParam_save,
                        xobs=xobs,yobs=yobs,xgrid=xgrid,burnin=burnin)
  class(resultData)="manyKnotsResults"
  
  return(resultData)

}






