# initialize_parms_spline=function(xobs,yobs,xgrid){


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



#   return(list(coefs=coefs,xtrue=xtrue,bases=bases,knotseq=knotseq,ytrue=ytrue,
#               lowept=lowept,upperept=upperept,designMatrixGrid=designMatrixGrid))
# }
