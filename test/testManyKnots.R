
library(tidyverse)
library(mvtnorm)
library(ggthemes)
parms=genData(3)
xobs=parms$xobs; yobs=parms$yobs; xgrid=parms$xgrid

initializeValues=initializeManyKnots(xobs,yobs,xgrid,dist_knots=2)
fit=runManyKnots(initializeValues)
plot(fit)
trace(fit)



summary(smoothParam_save)
summary(sigma2_save)

tmp=data.frame(exp(coefMat),iter=1:nrow(coefMat))
a1=gather(tmp,key,val,-iter)
ggplot(a1,aes(x=iter,y=val))+geom_line()+facet_wrap(~key)

fit=apply(fitMat,2,median)
plot(xobs,yobs,pch=16)
lines(xgrid,fit,col='green')

