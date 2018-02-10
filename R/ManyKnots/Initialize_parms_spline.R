findInitialSpline=function(xobs,bases,knotseq,yobs,designMatrix){
  
  min=999999
  for(i in seq(.5,7,by=.1)){
    coefs=seq(.1,i,length=ncol(designMatrix))
    y_mu=coefs%*%t(designMatrix)
    if(sum(abs(yobs-y_mu))<min){ min=sum(abs(yobs-y_mu)); save=i}
  }
  coef=seq(.1,save,length=ncol(designMatrix))
  return(coef)
}


initialize_parms_spline=function(xobs,yobs,xgrid){


  ### Initalize Isplines
  lowept=min(xobs)-.5
  upperept=max(xobs)+.5
  dist=1
  parms=Ispline(seq(min(xobs),max(xobs),by=dist),lowept,upperept)
  bases=parms$bases;
  knotseq=parms$knotseq

  #### Initial spline coefficients
  designMatrix=getIsplineC(xobs,knotseq,bases)
  coefs=findInitialSpline(xobs,bases,knotseq,yobs,designMatrix)
  y_mu=coefs%*%t(designMatrix)

  designMatrixGrid=getIsplineC(xgrid,knotseq,bases)



  return(list(coefs=coefs,bases=bases,knotseq=knotseq,y_mu=y_mu,
              lowept=lowept,upperept=upperept,designMatrixGrid=designMatrixGrid))
}
