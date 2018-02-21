priorNumberKnots=function(numIntKnots){
  
  lam=1
  
  x=seq(1:6)
  #   y=dpois(x,lam,log=TRUE)/(1-dpois(0,lam,log=TRUE))
  y=dpois(x,lam)/sum(dpois(1:6,lam))
  
  #   plot(x,y)
  
  
  #   return(dpois(numIntKnots,lam,log=TRUE)/(1-dpois(0,lam,log=TRUE)))
  return(log(y[numIntKnots]))
  
}

knotSelection=function(knotseq,bases,xobs,yobs,popmn1,popstd1,popprob1,xtrue1,xcens,ycens,
                       xsig,ysig,ytrue1,icoefs1,lowept,upperept){
  
  numInterior=length(knotseq)-8
  
  #only birth or relocate
  if(numInterior==1){
    choice=sample(c(1,3),1)
    if(choice==1){ parms=birthFunction(knotseq,bases,xobs,yobs,popmn1,popstd1,popprob1,xtrue1,xcens,ycens,
                                       xsig,ysig,ytrue1,icoefs1,lowept,upperept)
    }else{parms=relocateFunction(knotseq,bases,xobs,yobs,popmn1,popstd1,popprob1,xtrue1,xcens,ycens,
                                 xsig,ysig,ytrue1,icoefs1,lowept,upperept)}
    #birth, death, or relocate
  }else if(numInterior<=5){
    choice=sample(1:3,1)
    if(choice==1){ parms=birthFunction(knotseq,bases,xobs,yobs,popmn1,popstd1,popprob1,xtrue1,xcens,ycens,
                                       xsig,ysig,ytrue1,icoefs1,lowept,upperept)
    }else if(choice==2){ parms=deathFunction(knotseq,bases,xobs,yobs,popmn1,popstd1,popprob1,xtrue1,xcens,ycens,
                                             xsig,ysig,ytrue1,icoefs1,lowept,upperept)
    }else{ parms=relocateFunction(knotseq,bases,xobs,yobs,popmn1,popstd1,popprob1,xtrue1,xcens,ycens,
                                  xsig,ysig,ytrue1,icoefs1,lowept,upperept)}
    #death, relocate
  }else{
    choice=sample(2:3,1)
    if(choice==2){parms=deathFunction(knotseq,bases,xobs,yobs,popmn1,popstd1,popprob1,xtrue1,xcens,ycens,
                                      xsig,ysig,ytrue1,icoefs1,lowept,upperept)
    }else{ parms=relocateFunction(knotseq,bases,xobs,yobs,popmn1,popstd1,popprob1,xtrue1,xcens,ycens,
                                  xsig,ysig,ytrue1,icoefs1,lowept,upperept)}
  }
  #   
  #   choice=1
  #   parms=birthFunction(knotseq,bases,xobs,yobs,popmn1,popstd1,popprob1,xtrue1,xcens,ycens,
  #       xsig,ysig,ytrue1,icoefs1)
  #   
  knotseq=parms$knotseq
  bases=parms$bases
  accept=parms$accept
  ytrue1=parms$ytrue
  coef=parms$coef
  
  return(list(choice=choice,knotseq=knotseq,bases=bases,accept=accept,ytrue1=ytrue1,coef=coef))
  
}

deathFunction=function(knotseq,bases,xobs,yobs,popmn1,popstd1,popprob1,xtrue1,xcens,ycens,
                       xsig,ysig,ytrue1,icoefs1,lowept,upperept){
  
  #####delete knot
  lowE=knotseq[1]; upE=knotseq[length(knotseq)]
  knots=knotseq[c(-1,-2,-3)]
  knots=knots[seq(-length(knots),-(length(knots)-2),by=1)]
  numInterKnots=length(knots)-2
  t=sample(1:numInterKnots,1)
  deleteKnot=knots[t+1]
  oldRange=abs(max(lowept+.5,knots[t])-min(knots[t+2],upperept-.5))
  knots2=knots[-(t+1)]
  
  parms=Ispline(knots2[2:(length(knots2)-1)],lowE,upE)
  propBases=parms$bases;
  propKnotseq=parms$knotseq
  
  #   print(knotseq)
  #   print(knots)
  #   print(t)
  #   print(deleteKnot)
  #   print(propKnotseq)
  #   print(oldRange)
  
  
  
  #old coefficients
  designMatrix=getIsplineC(xtrue1,knotseq,bases)
  fit=lm(yobs~designMatrix-1)
  coefOldLS=as.numeric(coef(fit))
  covCoefOld=vcov(fit)
  oldY=as.numeric(icoefs1%*%t(designMatrix))
  
  
  #add noise to proposed coefficents and get y values
  designMatrix=getIsplineC(xtrue1,propKnotseq,propBases)
  fit=lm(yobs~designMatrix-1)
  coefNewLS=as.numeric(coef(fit))
  covCoefNew=vcov(fit)
  
  if(sum(is.na(coefNewLS))>0){
    return(list(knotseq=knotseq,bases=bases,accept=-1,ytrue=ytrue1,coef=icoefs1))    
  }  
  
  coefNew=tryCatch(as.numeric(rmvnorm(n=1,as.numeric(coefNewLS),sigma=covCoefNew)),
                   error=function(e){return(list(knotseq=knotseq,bases=bases,accept=-1,ytrue=ytrue1,coef=icoefs1))})
  newY=as.numeric(coefNew%*%t(designMatrix))
  
  ###test if new curve is monotonic
  designMatrix=getIsplineC(xgrid,propKnotseq,propBases)
  tempY=as.numeric(coefNew%*%t(designMatrix))  
  dataTemp=cbind(xgrid,tempY)
  dataTemp=dataTemp[order(xgrid),]
  ytemp=dataTemp[,2]
  isMono=1
  storage.mode(ytemp) <- "double"
  storage.mode(isMono) <- "integer"
  temp=.C("isMonotone",ytemp[xgrid>-4 & xgrid<4],length(xgrid[xgrid>-4 & xgrid<4]),isMono)
  if(temp[[3]]==0){
    return(list(knotseq=knotseq,bases=bases,accept=-1,ytrue=ytrue1,coef=icoefs1))    
  }
  
  
  #likelihoods and priors
  newLLK=sum(getLLK(xobs,yobs,xtrue1,newY,xcens,ycens,xsig,ysig))
  oldLLK=sum(getLLK(xobs,yobs,xtrue1,oldY,xcens,ycens,xsig,ysig))
  
  #dmnorm
  oldMV=dmnorm(as.numeric(icoefs1),as.numeric(coefOldLS),covCoefOld,log=TRUE)
  newMV=dmnorm(as.numeric(coefNew),as.numeric(coefNewLS),covCoefNew,log=TRUE)
  
  
  A=newLLK-oldLLK+oldMV-newMV+
    priorNumberKnots(numInterKnots-1)-priorNumberKnots(numInterKnots)+log(1/oldRange)
  
  if(is.na(A)==TRUE) A=-99
  
  #accept/reject
  if(log(runif(1)) < A){
    return(list(knotseq=propKnotseq,bases=propBases,accept=1,ytrue=newY,coef=coefNew))
  }else{
    return(list(knotseq=knotseq,bases=bases,accept=0,ytrue=ytrue1,coef=icoefs1))
  }
  
  
}

birthFunction=function(knotseq,bases,xobs,yobs,popmn1,popstd1,popprob1,xtrue1,xcens,ycens,
                       xsig,ysig,ytrue1,icoefs1,lowept,upperept){
  
  #####add knot
  lowE=knotseq[1]; upE=knotseq[length(knotseq)]
  knots=knotseq[c(-1,-2,-3)]
  knots=knots[seq(-length(knots),-(length(knots)-2),by=1)]
  numInterKnots=length(knots)-2
  interval=sample(1:(length(knots)-1),1)
  
  options(warn=-1)
  addKnot=runif(1,max(lowept+1.5,knots[interval]),min(knots[interval+1],upperept-1.5))
  newRange=abs(max(lowept+1.5,knots[interval])-min(knots[interval+1],upperept-1.5))
  rangeCoef=min(abs(knots[interval]-addKnot),abs(addKnot-knots[interval+1]))
  #   print(interval)
  #   print(addKnot)
  #   print(rangeCoef)
  #   print(addKnot)
  if(is.na(rangeCoef)==TRUE | rangeCoef<1) return(list(knotseq=knotseq,bases=bases,accept=-1,ytrue=ytrue1,coef=icoefs1))
  knots2=c(knots[1:interval],addKnot,knots[(interval+1):length(knots)])
  options(warn=1)
  
  parms=Ispline(knots2[2:(length(knots2)-1)],lowE,upE)
  propBases=parms$bases;
  propKnotseq=parms$knotseq
  
  #   print(knots)
  #   print(interval)
  #   print(addKnot)
  #   print(newRange)
  #   print(knots2)
  #   print(propKnotseq)
  #   print(numInterKnots)
  
  #old coefficients
  designMatrix=getIsplineC(xtrue1,knotseq,bases)
  fit=lm(yobs~designMatrix-1)
  coefOldLS=as.numeric(coef(fit))
  covCoefOld=vcov(fit)
  oldY=as.numeric(icoefs1%*%t(designMatrix))
  
  
  #add noise to proposed coefficents and get y values
  designMatrix=getIsplineC(xtrue1,propKnotseq,propBases)
  fit=lm(yobs~designMatrix-1)
  coefNewLS=as.numeric(coef(fit))
  covCoefNew=vcov(fit)
  
  if(sum(is.na(coefNewLS))>0){
    return(list(knotseq=knotseq,bases=bases,accept=-1,ytrue=ytrue1,coef=icoefs1))    
  }  
  
  coefNew=tryCatch(as.numeric(rmvnorm(n=1,as.numeric(coefNewLS),sigma=covCoefNew)),
                   error=function(e){return(list(knotseq=knotseq,bases=bases,accept=-1,ytrue=ytrue1,coef=icoefs1))})
  newY=as.numeric(coefNew%*%t(designMatrix))
  
  ###test if new curve is monotonic
  designMatrix=getIsplineC(xgrid,propKnotseq,propBases)
  tempY=as.numeric(coefNew%*%t(designMatrix))  
  dataTemp=cbind(xgrid,tempY)
  dataTemp=dataTemp[order(xgrid),]
  ytemp=dataTemp[,2]
  isMono=1
  storage.mode(ytemp) <- "double"
  storage.mode(isMono) <- "integer"
  temp=.C("isMonotone",ytemp[xgrid>-4 & xgrid<4],length(xgrid[xgrid>-4 & xgrid<4]),isMono)
  if(temp[[3]]==0){
    return(list(knotseq=knotseq,bases=bases,accept=-1,ytrue=ytrue1,coef=icoefs1))    
  }
  #   
  #   plot(xtrue1,newY)
  #   points(xtrue1,ytrue1,col=2)
  #   
  
  #likelihoods and priors
  newLLK=sum(getLLK(xobs,yobs,xtrue1,newY,xcens,ycens,xsig,ysig))
  oldLLK=sum(getLLK(xobs,yobs,xtrue1,oldY,xcens,ycens,xsig,ysig))
  
  #dmnorm
  oldMV=dmnorm(as.numeric(icoefs1),as.numeric(coefOldLS),covCoefOld,log=TRUE)
  newMV=dmnorm(as.numeric(coefNew),as.numeric(coefNewLS),covCoefNew,log=TRUE)
  
  
  A=newLLK-oldLLK+oldMV-newMV+
    priorNumberKnots(numInterKnots+1)-priorNumberKnots(numInterKnots)-log(1/newRange)
  
  if(is.na(A)==TRUE) A=-99
  
  #accept/reject
  if(log(runif(1)) < A){
    return(list(knotseq=propKnotseq,bases=propBases,accept=1,ytrue=newY,coef=coefNew))
  }else{
    return(list(knotseq=knotseq,bases=bases,accept=0,ytrue=ytrue1,coef=icoefs1))
  }
  
}

relocateFunction=function(knotseq,bases,xobs,yobs,popmn1,popstd1,popprob1,xtrue1,xcens,ycens,
                          xsig,ysig,ytrue1,icoefs1,lowept,upperept){
  
  
  #####sample knot
  lowE=knotseq[1]; upE=knotseq[length(knotseq)]
  intKnots=knotseq[c(-1,-2,-3,-4)]
  intKnots=intKnots[seq(-length(intKnots),-(length(intKnots)-3),by=1)]
  t=sample(1:(length(intKnots)),1)
  
  tempT=knotseq[c(-1,-2,-3,-length(knotseq),-(length(knotseq)-1),-(length(knotseq)-2))]
  options(warn=-1)
  propKnot=runif(1,max(lowept+1.5,tempT[t]),min(tempT[t+2],upperept-1.5))
  rangeCoef=min(abs(tempT[t]-propKnot),abs(propKnot-tempT[t+2]))
  if(is.na(propKnot)==TRUE | rangeCoef<1) return(list(knotseq=knotseq,bases=bases,accept=-1,ytrue=ytrue1,coef=icoefs1))
  options(warn=1)
  
  
  #get new knot sequence bases
  propKnotseq=knotseq
  propKnotseq[t+4]=propKnot
  parms=Ispline(c(propKnotseq[5:(length(intKnots)+4)]),lowE,upE)
  propBases=parms$bases;
  propKnotseq=parms$knotseq
  
  
  #old coefficients
  designMatrix=getIsplineC(xtrue1,knotseq,bases)
  fit=lm(yobs~designMatrix-1)
  coefOldLS=as.numeric(coef(fit))
  covCoefOld=vcov(fit)
  oldY=as.numeric(icoefs1%*%t(designMatrix))
  
  
  #add noise to proposed coefficents and get y values
  designMatrix=getIsplineC(xtrue1,propKnotseq,propBases)
  fit=lm(yobs~designMatrix-1)
  coefNewLS=as.numeric(coef(fit))
  covCoefNew=vcov(fit)
  
  if(sum(is.na(coefNewLS))>0){
    return(list(knotseq=knotseq,bases=bases,accept=-1,ytrue=ytrue1,coef=icoefs1))    
  }  
  
  coefNew=tryCatch(as.numeric(rmvnorm(n=1,as.numeric(coefNewLS),sigma=covCoefNew)),
                   error=function(e){return(list(knotseq=knotseq,bases=bases,accept=-1,ytrue=ytrue1,coef=icoefs1))})
  newY=as.numeric(coefNew%*%t(designMatrix))
  
  ###test if new curve is monotonic
  designMatrix=getIsplineC(xgrid,propKnotseq,propBases)
  tempY=as.numeric(coefNew%*%t(designMatrix))  
  dataTemp=cbind(xgrid,tempY)
  dataTemp=dataTemp[order(xgrid),]
  ytemp=dataTemp[,2]
  isMono=1
  storage.mode(ytemp) <- "double"
  storage.mode(isMono) <- "integer"
  temp=.C("isMonotone",ytemp[xgrid>-4 & xgrid<4],length(xgrid[xgrid>-4 & xgrid<4]),isMono)
  if(temp[[3]]==0){
    return(list(knotseq=knotseq,bases=bases,accept=-1,ytrue=ytrue1,coef=icoefs1))    
  }
  
  
  #likelihoods and priors
  newLLK=sum(getLLK(xobs,yobs,xtrue1,newY,xcens,ycens,xsig,ysig))
  oldLLK=sum(getLLK(xobs,yobs,xtrue1,oldY,xcens,ycens,xsig,ysig))
  
  #dmnorm
  oldMV=dmnorm(as.numeric(icoefs1),as.numeric(coefOldLS),covCoefOld,log=TRUE)
  newMV=dmnorm(as.numeric(coefNew),as.numeric(coefNewLS),covCoefNew,log=TRUE)
  
  
  A=newLLK-oldLLK+oldMV-newMV
  if(is.na(A)==TRUE) A=-99
  
  
  #accept/reject
  if(log(runif(1)) < A){
    return(list(knotseq=propKnotseq,bases=propBases,accept=1,ytrue=newY,coef=coefNew))
  }else{
    return(list(knotseq=knotseq,bases=bases,accept=0,ytrue=ytrue1,coef=icoefs1))
  }
}
