
plot.manyKnotsResults=function(resultData){

  fitMat=resultData$fitMat; xobs=resultData$xobs; yobs=resultData$yobs
  xgrid=resultData$xgrid

  gxDat=data_frame(xgrid,y=apply(fitMat,2,median))
  gxDat$lower=apply(fitMat,2,function(x) quantile(x,probs=c(.05)))
  gxDat$upper=apply(fitMat,2,function(x) quantile(x,probs=c(.95)))
  a1=data.frame(xobs,yobs)
  
  pltRel=ggplot(data=a1,aes(x=xobs,y=yobs))+geom_point()+
    geom_line(data=gxDat,aes(x=xgrid,y=y),color='darkblue',inherit.aes = FALSE)+
    geom_ribbon(data=gxDat,aes(x=xgrid,ymin=lower, ymax=upper),fill='cyan4',alpha=.5,inherit.aes = FALSE)+
    theme_fivethirtyeight(base_size=15)+
    theme(axis.line = element_line(colour = "black"))
  plot(pltRel)
  
}

trace <- function(x) UseMethod("trace")

trace.manyKnotsResults=function(resultData){
  
  coefMat=resultData$coefMat; sigma2=resultData$sigma2_save; smooth=resultData$smoothParam_save
  burnin=resultData$burnin
  
  a1=data.frame(coefMat,idx=1:nrow(coefMat))
  names(a1)=c(paste0("coef ",1:ncol(coefMat)),'idx')
  a1=a1 %>% gather(key,value,-idx)
  a1$key=factor(a1$key,levels=paste0("coef ",1:ncol(coefMat)))
  
  plotCoef=ggplot(data=a1,aes(x=idx,y=exp(value)))+geom_line()+facet_wrap(~key)+
    geom_vline(xintercept=burnin,color='red',linetype=2)+
    theme_fivethirtyeight(base_size=15)+
    annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf)+
    annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf)
  plot(plotCoef)
  
  
  a1=data.frame(VarianceParm=sigma2,SmoothParm=smooth,idx=1:length(sigma2))
  a1=a1 %>% gather(key,value,-idx)
  plotParms=ggplot(data=a1,aes(x=idx,y=value))+geom_line()+facet_wrap(~key,scale="free_y")+
    geom_vline(xintercept=burnin,color='red',linetype=2)+
    theme_fivethirtyeight(base_size=15)+
    annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf)+
    annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf)
  plot(plotParms)
  
}