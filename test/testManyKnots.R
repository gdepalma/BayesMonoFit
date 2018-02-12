
library(tidyverse)
library(mvtnorm)
parms=genData(3)
xobs=parms$xobs; yobs=parms$yobs; xgrid=parms$xgrid

parms=initializeSplineManyKnots(xobs,yobs,xgrid,dist_knots=1)
coefs=parms$coefs; bases=parms$bases; knotseq=parms$knotseq; y_mu=parms$y_mu
lowept=parms$lowept; upperept=parms$upperept; designMatrixGrid=parms$designMatrixGrid



summary(smoothParam_save)
summary(sigma2_save)

tmp=data.frame(exp(coefMat),iter=1:nrow(coefMat))
a1=gather(tmp,key,val,-iter)
ggplot(a1,aes(x=iter,y=val))+geom_line()+facet_wrap(~key)

fit=apply(fitMat,2,median)
plot(xobs,yobs,pch=16)
lines(xgrid,fit,col='green')