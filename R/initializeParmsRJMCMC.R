
initializeRJMCMC=function(xobs,yobs,xgrid){

  lowerept=min(xobs)
  upperept=max(xobs)
  nobs=length(xobs)
  
  ###find initial coef values and knot sequence
  tmp=(upperept-lowerept)/4
  parms=Ispline(c(tmp,tmp*3),lowerept,upperept)
  bases=parms$bases;
  knotseq=parms$knotseq
  designMatrix=getIsplineC(xobs,knotseq,bases)
  coefs=as.numeric(coef(lm(yobs~designMatrix-1)))
  y_mu=as.numeric(coefs%*%t(designMatrix))
  numInterior=length(knotseq)-8

  ### Loop over to find initial knotsequence and coefs that lead to monotonicity
  repeat{
    
    # test monotonicity
    tempY=y_mu 
    DT=data.table(xobs,tempY)
    DT=DT[order(xobs)]
    test_mono_y=DT[,2]
    if(all(test_mono_y == cummin(test_mono_y))==TRUE) break
    
    #select new interior knot locations and new coefs
    repeat{
      newIntKnots=numeric(2)
      newIntKnots[1]=runif(1,lowerept,upperept/2)
      newIntKnots[2]=runif(1,newIntKnots[1],upperept)
      parms=Ispline(newIntKnots,lowerept,upperept)
      bases=parms$bases;
      knotseq=parms$knotseq
      designMatrix=getIsplineC(xobs,knotseq,bases)
      coefs=as.numeric(coef(lm(yobs~designMatrix-1)))
      if(sum(is.na(coefs))==0) break
    }
  }


  designMatrixGrid=getIsplineC(xgrid,knotseq,bases)

  initializeValues=list(coefs=coefs,bases=bases,knotseq=knotseq,y_mu=y_mu,xgrid=xgrid,
                        designMatrix=designMatrix,designMatrixGrid=designMatrixGrid,
                        xobs=xobs,yobs=yobs)
  class(initializeValues)='initialRJMCMC'

  return(initializeValues)
}
