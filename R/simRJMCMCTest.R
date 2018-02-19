# 
# nsim=20000
# burnin=4000
# nobs=600
# 
# xsig=.707; ysig=2.121
# xgrid=seq(-6,6,length=500)
# xcens=rep(0,nobs)
# ycens=rep(0,nobs)
# 
# # density
# densTrue=rep(0,length(xgrid))
# popmn=c(-3,0,3); popstd=c(.7,.7,.7); popprob=c(1/3,1/3,1/3)
# for(i in 1:length(xgrid))
#   densTrue[i]=sum(popprob*dnorm(xgrid[i],popmn,popstd))
# densTrue=densTrue/sum(densTrue)
# 
# popmn=c(-3,0,3); popstd=c(.7,.7,.7); popprob=c(1/3,1/3,1/3)
# xtrue=rnormmix(n=nobs,lambda=popprob,mu=popmn,sigma=popstd)
# 
# ### X
# xobs=ceiling(xtrue+rnorm(nobs,0,xsig))
# 
# # Y
# coef=c(29,1.17,.48)
# ytrue=coef[1]*(exp(coef[2]-coef[3]*xtrue)/(1+exp(coef[2]-coef[3]*xtrue)))
# true=coef[1]*(exp(coef[2]-coef[3]*xgrid)/(1+exp(coef[2]-coef[3]*xgrid)))
# yobs=round(ytrue+rnorm(nobs,0,ysig))
# 
# # a1=data.frame(xobs,yobs,xcens,ycens)
# # a1=a1[order(a1$xobs,-a1$yobs),]
# # xobs=a1$xobs
# # yobs=a1$yobs
# # xcens=a1$xcens
# # ycens=a1$ycens
# xtrue1 = rnorm(length(xobs),xobs-.5,xsig)
# # xtrue1 = xtrue
# lowept=min(xobs)
# upperept=max(xobs)
# 
# 
# ### initalize Ispline
# lowept=min(xobs)-.5
# upperept=max(xobs)+.5
# if(sum(xcens==1)>0 | sum(ycens==1)>0)
#   upperept=max(xobs)+4
# if(sum(xcens==-1)>0 | sum(ycens==-1)>0)
#   lowept=min(xobs)-4
# 
# ###initial coef values and knot sequence
# parms=Ispline(c(-4,-2,0,2,4),lowept,upperept)
# bases=parms$bases;
# knotseq=parms$knotseq
# designMatrix=getIsplineC(xtrue1,knotseq,bases)
# icoefs1=as.numeric(coef(lm(yobs~designMatrix-1)))
# ytrue1=as.numeric(icoefs1%*%t(designMatrix))
# 
# 
# #### DP
# alphaDP=1
# popgrp1=rep(1,nobs)
# parms=updateDP(popmn1=NULL,popstd1=NULL,popgrp1,xtrue1,nobs,alphaDP,ncomp=NULL)
# popgrp1=parms[[1]]; popmn1=parms[[2]]; popstd1=sqrt(parms[[3]]); ncomp=parms[[4]]; popprob1=parms[[5]]
# 
# coefMat=matrix(nrow=nsim,ncol=length(icoefs1))
# fitMat=matrix(nrow=nsim-burnin,ncol=length(xgrid))
# MICDens=matrix(nrow=nsim-burnin,ncol=length(xgridDens))
# DIABrkpts=matrix(nrow=nsim-burnin,ncol=2)
# accept=rep(0,nsim)
# accept2=rep(0,nsim)
# choice=rep(0,nsim)
# fitMat=matrix(nrow=nsim-burnin,ncol=length(xgrid))
# MICDens=matrix(nrow=nsim-burnin,ncol=length(xgrid))
# DIABrkpts=matrix(nrow=nsim-burnin,ncol=2)
# 
# begin=Sys.time()
# 
# 
# for(isim in 1:nsim){
#    
# #     if(isim %% 1000==0) print(isim)
#   
#   ###DP
#   if(isim %% 500==0){
#     parms=updateDP(popmn1,popstd1,popgrp1,xtrue1,nobs,alphaDP,ncomp)
#     popgrp1=parms[[1]]; popmn1=parms[[2]]; popstd1=sqrt(parms[[3]]); ncomp=parms[[4]]; popprob1=parms[[5]]
#   }
#   
#   ### Updating the xtrue vector
#   if(isim %% 25==0 | isim==1){
#     parms=updateSplineXtrue(xtrue1,ytrue1,knotseq,bases,lowept,upperept,icoefs1,popgrp1,popmn1,popstd1,popprob1,
#             xcens,ycens,xsig,ysig,xobs,yobs,nobs)
#     xtrue1=parms$xtrue1; ytrue1=as.numeric(parms$ytrue1)
#   }
# #   
#   #get covariance
#   designMatrix=getIsplineC(xtrue1,knotseq,bases)
#   fit=lm(yobs~designMatrix-1)
#   Sigma=vcov(fit)
# 
#   #### joint adaptive update
#   propCoef=as.numeric(rmvnorm(n=1,as.numeric(icoefs1),sigma=Sigma))
#   designMatrix=getIsplineC(xtrue1,knotseq,bases)
#   newY=as.numeric(propCoef%*%t(designMatrix))
#   
#   #test mono
#   designMatrix=getIsplineC(xgrid,knotseq,bases)
#   tempY=as.numeric(propCoef%*%t(designMatrix))  
#   dataTemp=cbind(xgrid,tempY)
#   dataTemp=dataTemp[order(xgrid),]
#   ytemp=dataTemp[,2]
#   isMono=1
#   storage.mode(ytemp) <- "double"
#   storage.mode(isMono) <- "integer"
#   temp=.C("isMonotone",ytemp[xgrid>-4 & xgrid<4],length(xgrid[xgrid>-4 & xgrid<4]),isMono)
#   
#   if(temp[[3]]==1){
#     #log likelihood and priors
#     oldLLK=sum(getLLK(xobs,yobs,xtrue1,ytrue1,xcens,ycens,xsig,ysig))
#     newLLK=sum(getLLK(xobs,yobs,xtrue1,newY,xcens,ycens,xsig,ysig))
#     #accept/reject
#     if(log(runif(1)) < newLLK-oldLLK){
#       ytrue1=newY
#       icoefs1=propCoef
#       accept[isim]=1
#     }
#   }else{accept[isim]=-1}
#   
#   
#   ### Update knot location                               
#   parms=knotSelection(knotseq,bases,xobs,yobs,popmn1,popstd1,popprob1,xtrue1,xcens,ycens,
#            xsig,ysig,ytrue1,icoefs1,lowept,upperept)
#   knotseq=parms$knotseq
#   bases=parms$bases 
#   icoefs1=parms$coef
#   ytrue1=parms$ytrue1
#   accept2[isim]=parms$accept
#   choice[isim]=parms$choice
#   
#   
#   ###save
#   if(isim>burnin){
#     designMatrixGrid=getIsplineC(xgrid,knotseq,bases)
#     fit=icoefs1%*%t(designMatrixGrid)
#     fitMat[isim-burnin,]=fit
#     MICDens[isim-burnin,]=findDens(xgrid,popmn1,popstd1,popprob1)
#   }
#   
# }
# 
# #Thin
# MICDens=MICDens[seq(1,nrow(MICDens),by=2),]
# fitMat=fitMat[seq(1,nrow(fitMat),by=2),]
# MICDens=apply(MICDens,2,median)
# fit=apply(fitMat,2,median)
# 
# plot(accept)
# plot(accept2)
# hist(xtrue1)
# 
# data1=data.frame(xgrid,MICDens,densTrue)
# a1=melt(data1,id=1)
# fit1=ggplot(a1,aes(xgrid,value,color=variable))+geom_line()
# plot(fit1)
# 
# trueFit=true
# data1=data.frame(xgrid,fit,trueFit)
# a1=melt(data1,id=1)
# fit2=ggplot(a1,aes(xgrid,value,color=variable))+geom_line()
# plot(fit2)
# 
# 
