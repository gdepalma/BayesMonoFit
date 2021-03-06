findInitialCoefManyKnots=function(xobs,bases,knotseq,yobs,designMatrix){
  
  min=999999
  for(i in seq(.5,7,by=.1)){
    coefs=seq(.1,i,length=ncol(designMatrix))
    y_mu=coefs%*%t(designMatrix)
    if(sum(abs(yobs-y_mu))<min){ min=sum(abs(yobs-y_mu)); save=i}
  }
  coef=seq(.1,save,length=ncol(designMatrix))
  return(coef)
}


initializeManyKnots=function(xobs,yobs,xgrid,dist_knots=1){


  ### Initalize Isplines
  lowerept=min(xobs)-.5
  upperept=max(xobs)+.5
  parms=Ispline(seq(min(xobs),max(xobs),by=dist_knots),lowerept,upperept)
  bases=parms$bases;
  knotseq=parms$knotseq

  #### Initial spline coefficients
  designMatrix=getIsplineC(xobs,knotseq,bases)
  coefs=findInitialCoefManyKnots(xobs,bases,knotseq,yobs,designMatrix)
  y_mu=coefs%*%t(designMatrix)

  designMatrixGrid=getIsplineC(xgrid,knotseq,bases)

  initializeValues=list(coefs=coefs,bases=bases,knotseq=knotseq,y_mu=y_mu,xgrid=xgrid,
                        designMatrix=designMatrix,designMatrixGrid=designMatrixGrid,
                        xobs=xobs,yobs=yobs)
  class(initializeValues)='initialManyKnots'

  return(initializeValues)
}
