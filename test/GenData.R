

genData=function(type,nobs=200,sigma2=2){
  
  nobs=500
  
  popmn=c(-4.6,-2,1); popstd=c(1.1,1.5,1.5); popprob=c(.6,.2,.2)
#   popmn=c(-4,0,4); popstd=c(.7,.7,.7); popprob=c(1/3,1/3,1/3)
#   popmn=c(-6,-3,3); popstd=c(1.5,1.5,.7); popprob=c(.6,.3,.1)
#   popmn=c(-3,0,3); popstd=c(1,1,1); popprob=c(.5,.3,.2)
  
  # components <- sample(1:3,prob=popprob,size=nobs,replace=TRUE)
  # xobs <- rnorm(n=nobs,mean=popmn[components],sd=popstd[components])
  xobs=runif(nobs,0,20)
  xtrue=xobs
  
  ### X
  xgrid=seq(min(xobs),max(xobs),length=1000)
  
  ## Y
  
  #1
  if(type==1){
    ytrue=30-2*xtrue
    true=30-2*xgrid 
  }
  
  #2
  if(type==2){
    coef=c(29,1.17,.48)
    ytrue=coef[1]*(exp(coef[2]-coef[3]*xtrue)/(1+exp(coef[2]-coef[3]*xtrue)))
    true=coef[1]*(exp(coef[2]-coef[3]*xgrid)/(1+exp(coef[2]-coef[3]*xgrid)))
  }
  
  #3
  if(type==3){
    icoefsT=c(1,1,20,1,20,1)
    parms=Ispline(c(-3,0,10),-6,15)
    basesT=parms$bases;
    knotseqT=parms$knotseq
    designMatrix=getIsplineC(xtrue,knotseqT,basesT)
    ytrue = icoefsT%*%t(designMatrix)+5
    designMatrix=getIsplineC(xgrid,knotseqT,basesT)
    true = icoefsT%*%t(designMatrix)+5  
  }
  
  #4
  if(type==4){
    icoefsT=c(1,10,1,25,1,1,10,1)
    parms=Ispline(c(-4,-2,0,2,4),-6,6)
    basesT=parms$bases;
    knotseqT=parms$knotseq
    designMatrix=getIsplineC(xtrue,knotseqT,basesT)
    ytrue = icoefsT%*%t(designMatrix)
    designMatrix=getIsplineC(xgrid,knotseqT,basesT)
    true = icoefsT%*%t(designMatrix)
  }
  if(type==5){
    coef=c(35,1.17,.1,1.2)
    mb = (2*coef[3]*coef[4])/(coef[3]+coef[4])
    fx = 1/(1+exp(-mb*(coef[2]-xtrue)))
    ytrue=coef[1]*(fx*exp(coef[3]*(coef[2]-xtrue))+(1-fx)*exp(coef[4]*
        (coef[2]-xtrue)))/(1+fx*exp(coef[3]*(coef[2]-xtrue))+(1-fx)*
        exp(coef[4]*(coef[2]-xtrue)))
    fx = 1/(1+exp(-mb*(coef[2]-xgrid)))
    true=coef[1]*(fx*exp(coef[3]*(coef[2]-xgrid))+(1-fx)*exp(coef[4]*
          (coef[2]-xgrid)))/(1+fx*exp(coef[3]*(coef[2]-xgrid))+(1-fx)*
           exp(coef[4]*(coef[2]-xgrid)))
  }
  
  yobs=as.numeric(ytrue+rnorm(nobs,0,sigma2))
  
  a1=data.frame(xobs,yobs)

  plt=ggplot(a1,aes(x=xobs,y=yobs))+geom_point()+
    theme_fivethirtyeight(base_size=15)
  plot(plt)
  
  return(list(xobs=xobs,yobs=yobs,true=true,xtrue=xtrue,xgrid=xgrid))
}

