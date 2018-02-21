
parms=genData(3)
xobs=parms$xobs; yobs=parms$yobs; xgrid=parms$xgrid

initializeValues=initializeManyKnots(xobs,yobs,xgrid,dist_knots=2)
fit=runManyKnots(initializeValues)
plot(fit)
trace(fit)




parms=genData(3)
xobs=parms$xobs; yobs=parms$yobs; xgrid=parms$xgrid

initializeValues=initializeRJMCMC(xobs,yobs,xgrid)

