



###Functions

GoFC<-function(M, y, na.rm=TRUE)
{
 # resS<-(y-M$summary.fitted.values$mean[1:length(y)])*(M$summary.hyperpar[1,]$mean)**0.5/M$summary.fitted.values$mean[1:length(y)]
  rmse<-(mean(((y-M$summary.fitted.values$mean[1:length(y)])/10)**2, na.rm=na.rm))**(0.5)
  cmsr<-c(M$dic$dic, M$waic$waic, rmse) 
  return(cmsr)
}



#########################################
#### Models with no spatial effects
#########################################


library(INLA)

#### Reference model with constant parameters

formula<-(IncS~1)
MM0<-inla(formula, data=datosvC2, 
          family="gamma", 
          control.compute = list(dic = TRUE,waic=TRUE),
          control.predictor=list(compute=T,link=1),
          control.family = list(control.link=list(model="log")))
summary(MM0)
round(GoFC(M=MM0, y=datosvC2$IncS),2)


formula<-(IncS~1)
MM0exp<-inla(formula, data=datosvC2, 
         family="exponential", 
         control.compute = list(dic = TRUE,waic=TRUE),
         control.predictor=list(compute=T,link=1),
         control.family = list(control.link=list(model="log")))
round(GoFC(M=MM0exp, y=datosvC2$IncS),2)

#### Model with covariates

formula<-(IncS~1+Cov)
MM1<-inla(formula, data=datosvC2, 
        family="gamma", 
        control.compute = list(dic = TRUE,waic=TRUE),
        control.predictor=list(compute=T,link=1),
        control.family = list(control.link=list(model="log")))
summary(MM1)
round(GoFC(M=MM1, y=datosvC2$IncS),2)




formula<-(IncS~1+Cov)
MM1exp<-inla(formula, data=datosvC2, 
         family="exponential", 
         control.compute = list(dic = TRUE,waic=TRUE),
         control.predictor=list(compute=T,link=1),
         control.family = list(control.link=list(model="log")))
summary(MM1exp)
round(GoFC(M=MM1exp, y=datosvC2$IncS),2)



#### Model with covariates and random seasonal effects

formula<-(IncS~1+Cov+f(dayInd, model="seasonal", season.length=92))

MM3exp<-inla(formula, data=datosvC2, 
          family="exponential", 
          control.compute = list(dic = TRUE,waic=TRUE),
          control.predictor=list(compute=T,link=1),
          control.family = list(control.link=list(model="log")))

MM3<-inla(formula, data=datosvC2, 
         family="gamma", 
         control.compute = list(dic = TRUE,waic=TRUE),
         control.predictor=list(compute=T,link=1),
         control.family = list(control.link=list(model="log")),
         control.mode = list(theta =c(MM1$mode$theta,MM3exp$mode$theta),restart=FALSE))
summary(MM3)
round(GoFC(M=MM3, y=datosvC2$IncS),2)


#########################################
#### Models with spatial effects
#########################################

  
bnd<-inla.nonconvex.hull(Ucoords)
mesh<-inla.mesh.2d(boundary=bnd,max.edge=c(1, 4))
plot(mesh)
points(coords, pch=16, cex=0.7)
lines(borderS, col="red", lw=3)

range0<-min(diff(range(coords[,1])),diff(range(coords[,2])))/2
spdepc<-inla.spde2.pcmatern(mesh=mesh, prior.range=c(range0,0.5),prior.sigma=c(1,0.1))
G<-spdepc$n.spde


#### Model with only spatial intercept

A<-inla.spde.make.A(mesh, loc=coords)
#dim(A) is  T*L*S  x G   61*92*40=224480 x 363

formula<-(IncS~-1+Intercept+f(spatial.field, model=spdepc))
MMS0<-inla(formula, data=list(IncS=datosvC2$IncS,Intercept=rep(1,G), spatial.field=c(1:G)), 
          family="gamma", 
          control.compute = list(dic = TRUE,cpo=TRUE,waic=TRUE),
          control.predictor=list(A=A,compute=T,link=1),
          control.family = list(control.link=list(model="log"))
)
summary(MMS0)
round(GoFC(M=MMS0, y=datosvC2$IncS),2)


#### Model with spatial intercept and covariates

s.index1<-inla.spde.make.index (name="spatialInd", n.spde=spdepc$n.spde)
stack1<-inla.stack(data=list(IncS=datosvC2$IncS), tag="est",
                      A=list(A,1),
                      effects=list(c(s.index1, list(Intercept=1)), 
                                   list(Cov=datosvC2$Cov, yearInd=datosvC2$yearInd,
                                    dayInd=datosvC2$dayInd,timeInd=datosvC2$timeInd,
                                    siteInd=datosvC2$siteInd))
                   )

# dim(stack1$A): [1]  224480 224566  ( nT*nL*nS=224480); 


formula<-(IncS~-1+Intercept+Cov+f(spatialInd, model=spdepc))
MMS1<-inla(formula, data=inla.stack.data(stack1,spde=spdepc), 
          family="gamma", 
          control.compute = list(dic = TRUE,cpo=TRUE,waic=TRUE),
          control.predictor=list(A=inla.stack.A(stack1),compute=T,link=1),
          control.family = list(control.link=list(model="log"))
)
summary(MMS1)
round(GoFC(M=MMS1, y=datosvC2$IncS),2)



#### Model with spatial intercept, covariates  and  iid year random effects

formula<-(IncS~-1+Intercept+Cov+f(spatialInd, model=spdepc)+f(yearInd, model="iid"))
MMS3<-inla(formula, data=inla.stack.data(stack1,spde=spdepc), 
          family="gamma", 
          control.compute = list(dic = TRUE,cpo=TRUE,waic=TRUE),
          control.predictor=list(A=inla.stack.A(stack1),compute=T,link=1),
          control.family = list(control.link=list(model="log"))
)
summary(MMS3)
round(GoFC(M=MMS3, y=datosvC2$IncS),2)


#### Model with spatial intercept, covariates  and  daily random effects

formula<-(IncS~-1+Intercept+Cov+f(spatialInd, model=spdepc)+
          f(timeInd, model="iid")) #very similar with seasonal instead of iid
MMS7<-inla(formula, data=inla.stack.data(stack1,spde=spdepc), 
          family="gamma", 
          control.compute = list(dic = TRUE,cpo=TRUE,waic=TRUE),
          control.predictor=list(A=inla.stack.A(stack1),compute=T,link=1),
          control.family = list(control.link=list(model="log"))
)
summary(MMS7)
round(GoFC(M=MMS7, y=datosvC2$IncS),2)


#### Model with spatial intercept, covariates  and  AR(1) daily random effects

formula<-(IncS~-1+Intercept+Cov+f(spatialInd, model=spdepc)+
            f(dayInd, model="ar1", replicate=yearInd))
MMS7arexp<-inla(formula, data=inla.stack.data(stack1,spde=spdepc), 
              family="exponential", 
              control.compute = list(dic = TRUE,waic=TRUE,cpo=TRUE),
              #,return.marginals = FALSE, return.marginals.predictor = FALSE),
              control.predictor=list(A=inla.stack.A(stack1),compute=T,link=1),
              control.family = list(control.link=list(model="log"))
)
summary(MMS7arexp)
round(GoFC(M=MMS7arexp, y=datosvC2$IncS),2)


formula<-(IncS~-1+Intercept+Cov+f(spatialInd, model=spdepc)+
            f(dayInd, model="ar1", replicate=yearInd))
thetaMS7ar<-c(MMS7$mode$theta[1], MMS7arexp$mode$theta)
MMS7ar<-inla(formula, data=inla.stack.data(stack1,spde=spdepc), 
           family="gamma", 
           control.compute = list(dic = TRUE,waic=TRUE,cpo=TRUE),
           control.predictor=list(A=inla.stack.A(stack1),compute=T,link=1),
           control.family = list(control.link=list(model="log")),
          control.mode = list(theta =thetaMS7ar,restart=FALSE)
)
summary(MMS7ar)
round(GoFC(M=MMS7ar, y=datosvC2$IncS),2)


# same models with arguments in control.compute to   obtain a model valid for generating samples

MSel<-MMS7ar

formula<-(IncS~-1+Intercept+Cov+f(spatialInd, model=spdepc)+
            f(dayInd, model="ar1", replicate=yearInd))
MSelT<-inla(formula, data=inla.stack.data(stack1,spde=spdepc), 
            family="gamma", 
            control.compute = list(config = TRUE),
            control.predictor=list(A=inla.stack.A(stack1),compute=T,link=1),
            control.family = list(control.link=list(model="log")),
            control.mode = list(theta =thetaMS7ar,restart=FALSE)
)

save(MSel, MSelT,datosvC2,mesh, A, stack1,nT, nL, nS,
    file="G:\\Mi unidad\\Espacial\\E9RecordsExpl\\Datos_Z\\Selected_Model.RData")


##############
### GOF Function for plotting residuals
#########################################


GOFPlots<-function(M, y=datosvC2$IncS,t=datosvC2$yearInd,l=datosvC2$dayInd,
	s=datosvC2$siteInd,vlim=NULL,yy=NULL, prow=c(3,3))
{

#Residuals and fitted aand observed values are re-scale to degrees C

	resS<-(y-M$summary.fitted.values$mean[1:length(y)])*(M$summary.hyperpar[1,]$mean)**0.5/M$summary.fitted.values$mean[1:length(y)]

	par(mfrow=prow)
	if (is.null(vlim)) vlim<-c(min(c(resS) , na.rm=T),max(c(resS) , na.rm=T))
	plot(t,resS, cex=0.5,pch=16, xlab="Year", ylab="Residuals", ylim=vlim,xaxt="n")
	axis(1,at=c(2,12,22,32,42,52,62), labels=seq(1961,2021, by=10))
	plot(l,resS, cex=0.5,pch=16, xlab="Day", ylab="Residuals", ylim=vlim)
	abline(v=c(30,61), col="gray")
	plot(s,resS, cex=0.5,pch=16, xlab="Site", ylab="Residuals", ylim=vlim)

	fitm<-M$summary.fitted.values$mean[1:length(y)]/10
	if (is.null(yy)) yy<-max(fitm, na.rm=T)
	print(summary(fitm))
	plot(t,fitm, cex=0.5,pch=16, xlab="Year", ylab="Posterior mean",ylim=c(0,yy),xaxt="n")
	axis(1,at=c(2,12,22,32,42,52,62), labels=seq(1961,2021, by=10))
	plot(l,fitm, cex=0.5,pch=16, xlab="Day",
	     ylab="Posterior mean", ylim=c(0,yy))
	abline(v=c(30,61), col="gray")
	plot(s,fitm, cex=0.5,pch=16, xlab="Site", 
	     ylab="Posterior mean",ylim=c(0,yy))

	plot(fitm,y/10, pch=16, cex=0.7,xlab="Posterior mean",ylab="Obs. Inc.",xlim=c(0,yy))
 	lines(y,y, col="grey")
}




#plotting the residuals

#pdf(file="G:/Mi unidad/Espacial/E9RecordsExpl/Graphs/ResidualsMS7ar.pdf")
GOFPlots(M=MMS7ar,y=datosvC2$IncS,t=datosvC2$yearInd,l=datosvC2$dayInd,s=datosvC2$siteInd, prow=c(3,3),vlim=c(0,10), yy=11))
#dev.off()
#pdf(file="G:/Mi unidad/Espacial/E9RecordsExpl/Graphs/ResidualsM1.pdf")
GOFPlots(M=MM1,y=datosvC2$IncS,t=datosvC2$yearInd,l=datosvC2$dayInd,s=datosvC2$siteInd, prow=c(3,3),vlim=c(0,10), yy=11))
#dev.off()



