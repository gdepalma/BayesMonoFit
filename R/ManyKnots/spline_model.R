# fitSplineMon=function(xobs,yobs,coefs,x_est,y_mu,xgrid,lowept,upperept,knotseq,bases,designMatrixGrid,
#                                     xsig=1,numIter=22000,burnin=14000,thin=4){

  numIter=6000
  burnin=2000
  fitMat=matrix(nrow=numIter-burnin,ncol=length(xgrid))
  coefMat=matrix(nrow=numIter,ncol=length(coefs))
  sigma2_save=rep(NA,numIter)
  smoothParam_save=rep(NA,numIter)
  acceptCoef=rep(NA,numIter)
  
  
  nobs=length(xobs)

  smoothParam=2
  
  begin=Sys.time()


  for(iter in 1:numIter){

    if(iter %% 500==0){
      cat('iter: ',iter,'\n')
      cat(Sys.time()-begin,'\n')
    }
    
    ### Update Coefs
    parms=updateCoefs(smoothParam,coefs,iter,bases,y_mu,knotseq,xobs,yobs,sigma2,coefMat)
    coefs=parms$coefs; y_mu=as.numeric(parms$y_mu); acceptCoef[iter]=parms$accept


    ### Update smoothness parameter
    parms=updateSmoothParm(smoothParam,coefs)
    smoothParam=parms$smooth
    
    ### Update error variance
    sigma2=updateVariance(xobs,yobs,y_mu,sigma2,smoothParam,coefs)

    # smoothParam_sav[iter]=smoothParam
    coefMat[iter,]=log(coefs)
    sigma2_save[iter]=sigma2
    smoothParam_save[iter]=smoothParam
    if(iter>burnin){
      fitMat[iter-burnin,]=coefs%*%t(designMatrixGrid)
    }

  }

  summary(smoothParam_save)
  summary(sigma2_save)
  
  tmp=data.frame(exp(coefMat),iter=1:nrow(coefMat))
  a1=gather(tmp,key,val,-iter)
  ggplot(a1,aes(x=iter,y=val))+geom_line()+facet_wrap(~key)
  
  fit=apply(fitMat,2,median)
  plot(xobs,yobs,pch=16)
  lines(xgrid,fit,col='green')


  return(list(MICDens=MICDens,fitMat=fitMat))

}






