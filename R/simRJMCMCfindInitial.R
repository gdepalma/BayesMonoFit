

simRJMCMC=function(xobs,yobs,xcens,ycens){
  
  nsim=8000
  burnin=4000
  
  xsig=sqrt(.5^2+.5^2)
  ysig=sqrt(1.5^2+1.5^2)
  
  a1=data.frame(xobs,yobs,xcens,ycens)
  a1=a1[order(a1$xobs,-a1$yobs),]
  xobs=a1$xobs
  yobs=a1$yobs
  xcens=a1$xcens
  ycens=a1$ycens
  nobs=length(xobs)
  xtrue1 = rnorm(length(xobs),xobs-.5,xsig)
  lowept=min(xobs)
  upperept=max(xobs)
  
  
  ### initalize Ispline
  lowept=min(xobs)-.5
  upperept=max(xobs)+.5
  if(sum(xcens==1)>0 | sum(ycens==1)>0)
    upperept=max(xobs)+4
  if(sum(xcens==-1)>0 | sum(ycens==-1)>0)
    lowept=min(xobs)-4
  
  ###find initial coef values and knot sequence
  parms=Ispline(c(-3,3),lowept,upperept)
  bases=parms$bases;
  knotseq=parms$knotseq
  designMatrix=getIsplineC(xtrue1,knotseq,bases)
  icoefs1=as.numeric(coef(lm(yobs~designMatrix-1)))
  ytrue1=as.numeric(icoefs1%*%t(designMatrix))
  
  intKnots=knotseq[c(-1,-2,-3,-4)]
  intKnots=intKnots[seq(-length(intKnots),-(length(intKnots)-3),by=1)]
  numInterior=length(intKnots)
  lowE=knotseq[1]; upE=knotseq[length(knotseq)]
  
  repeat{
    
    designMatrix=getIsplineC(xtrue1,knotseq,bases)
    tempY=as.numeric(icoefs1%*%t(designMatrix))  
    dataTemp=cbind(xtrue1,tempY)
    dataTemp=dataTemp[order(xtrue1),]
    ytemp=dataTemp[,2]
  
    #test mono
    isMono=1
    storage.mode(ytemp) <- "double"
    storage.mode(isMono) <- "integer"
    temp=.C("isMonotone",ytemp[xtrue1>-3.5 & xtrue1<2.5],length(xtrue1[xtrue1>-3.5 & xtrue1<2.5]),isMono)
    if(temp[[3]]==1) break
    
    #select new interior knot locations and new coefs
    repeat{
      newIntKnots=numeric(2)
      newIntKnots[1]=runif(1,lowE,upE-.5)
      newIntKnots[2]=runif(1,newIntKnots[1],upE)
      parms=Ispline(newIntKnots,lowE,upE)
      bases=parms$bases;
      knotseq=parms$knotseq
      designMatrix=getIsplineC(xtrue1,knotseq,bases)
      icoefs1=as.numeric(coef(lm(yobs~designMatrix-1)))
      if(sum(is.na(icoefs1))==0) break
    }
    print(knotseq)
  }
#   print('Initial Values Found')
#   print(knotseq)
#   print(icoefs1)
  
  #### DP
  alphaDP=1
  popgrp1=rep(1,nobs)
  parms=updateDP(popmn1=NULL,popstd1=NULL,popgrp1,xtrue1,nobs,alphaDP,ncomp=NULL)
  popgrp1=parms[[1]]; popmn1=parms[[2]]; popstd1=sqrt(parms[[3]]); ncomp=parms[[4]]; popprob1=parms[[5]]
  
  
  
  
  coefMat=matrix(nrow=nsim,ncol=length(icoefs1))
  fitMat=matrix(nrow=nsim-burnin,ncol=length(xgrid))
  MICDens=matrix(nrow=nsim-burnin,ncol=length(xgrid))
  DIABrkpts=matrix(nrow=nsim-burnin,ncol=2)
  accept=rep(0,nsim)
  accept2=rep(0,nsim)
  choice=rep(0,nsim)
  fitMat=matrix(nrow=nsim-burnin,ncol=length(xgrid))
  MICDens=matrix(nrow=nsim-burnin,ncol=length(xgrid))
  DIABrkpts=matrix(nrow=nsim-burnin,ncol=2)
  
  begin=Sys.time()
  
  
  for(isim in 1:nsim){
     
#     if(isim %% 1000==0) print(isim)
    
    ###DP
    if(isim %% 500==0){
      parms=updateDP(popmn1,popstd1,popgrp1,xtrue1,nobs,alphaDP,ncomp)
      popgrp1=parms[[1]]; popmn1=parms[[2]]; popstd1=sqrt(parms[[3]]); ncomp=parms[[4]]; popprob1=parms[[5]]
    }
    
    #get covariance
    designMatrix=getIsplineC(xtrue1,knotseq,bases)
    fit=lm(yobs~designMatrix-1)
    Sigma=vcov(fit)
  
    #joint adaptive update
    propCoef=as.numeric(rmvnorm(n=1,as.numeric(icoefs1),sigma=Sigma))
    designMatrix=getIsplineC(xtrue1,knotseq,bases)
    newY=as.numeric(propCoef%*%t(designMatrix))
    
    #test mono
    designMatrix=getIsplineC(xgrid,knotseq,bases)
    tempY=as.numeric(propCoef%*%t(designMatrix))  
    dataTemp=cbind(xgrid,tempY)
    dataTemp=dataTemp[order(xgrid),]
    ytemp=dataTemp[,2]
    isMono=1
    storage.mode(ytemp) <- "double"
    storage.mode(isMono) <- "integer"
    temp=.C("isMonotone",ytemp[xgrid>-3.5 & xgrid<2.5],length(xgrid[xgrid>-3.5 & xgrid<2.5]),isMono)
    
    if(temp[[3]]==1){
      #log likelihood and priors
      oldLLK=sum(getLLK(xobs,yobs,xtrue1,ytrue1,xcens,ycens,xsig,ysig))
      newLLK=sum(getLLK(xobs,yobs,xtrue1,newY,xcens,ycens,xsig,ysig))
      #accept/reject
      if(log(runif(1)) < newLLK-oldLLK){
        ytrue1=newY
        icoefs1=propCoef
        accept[isim]=1
      }
    }else{accept[isim]=-1}
    
    
    ### Update knot location                               
    parms=knotSelection(knotseq,bases,xobs,yobs,popmn1,popstd1,popprob1,xtrue1,xcens,ycens,
             xsig,ysig,ytrue1,icoefs1,lowept,upperept)
    knotseq=parms$knotseq
    bases=parms$bases 
    icoefs1=parms$coef
    ytrue1=parms$ytrue1
    accept2[isim]=parms$accept
    choice[isim]=parms$choice
    
    
    ###save
    if(isim>burnin){
      designMatrixGrid=getIsplineC(xgrid,knotseq,bases)
      fit=icoefs1%*%t(designMatrixGrid)
      fitMat[isim-burnin,]=fit
      MICDens[isim-burnin,]=findDens(xgrid,popmn1,popstd1,popprob1)
    }
    
  }
  
  #Thin
  MICDens=MICDens[seq(1,nrow(MICDens),by=2),]
  fitMat=fitMat[seq(1,nrow(fitMat),by=2),]
  
#   plot(accept)
#   plot(accept2)
#   plot(accept2[choice==3])
#   
  MICDens=apply(MICDens,2,median)
  fit=apply(fitMat,2,median)
  
#   plot(xgrid,fit)
#   points(xobs,yobs)
#   
  M1Test=-1; M2Test=1
  weights=MICDens/sum(MICDens)
  parms=findDIAC(yobs,xgrid,weights,fit,M1Test,M2Test,xsig,ysig,3,20)
  trueD1=parms[[1]]
  trueD2=parms[[2]]
  
  #   plot(xobs,yobs)
  #   lines(xgrid,apply(fitMat,2,median))
  
  return(list(D1=trueD1,D2=trueD2))
}
