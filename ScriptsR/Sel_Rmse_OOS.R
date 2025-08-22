
library(INLA)


mii<-cbind(c(8,28,29,34),c(6,18,22,40),c(7,15,26,38),c(13,16,35,39),c(2,9,10,23),
           c(3,11,19,27),c(5,12,30,37),c(1,4,17,24),c(25,31,33,36),c(14,20,21,32)  )

### Selecting non-spatial models

Sel_Rmse_OOS<-function(formula, ii,  ddatosvC2=datosvC2)
{
  ind<-((ddatosvC2$siteInd==ii[1])|(ddatosvC2$siteInd==ii[2])|(ddatosvC2$siteInd==ii[3])|(ddatosvC2$siteInd==ii[4]) )
  #ind is indicator of observations from out of sample locations. length: 61*92*40
  indS<-(ind==FALSE)
  #indS is indicator of observations from  sample locations. length: 61*92*40
  IncSC<-ddatosvC2$IncS  # Inc data complete (all the stations)

  ddatosvC2$IncS[ind]<-NA

  M<-inla(formula, data=ddatosvC2, 
          family="gamma", 
          control.compute = list(dic = FALSE,waic=FALSE),
          control.predictor=list(compute=T,link=1),
          control.family = list(control.link=list(model="log")))
  
  resS<-(IncSC[indS]-M$summary.fitted.values$mean[indS])/10 # to express the residuals in degrees
  yearS<-ddatosvC2$yearInd[indS] #year of the residuals S
  resOS<-(IncSC[ind]-M$summary.fitted.values$mean[ind])/10
  yearOS<-ddatosvC2$yearInd[ind] #year of the residuals OS
  mseOS<-sum(resOS**2, na.rm=TRUE)
  nOS<-sum(!is.na(resOS))
  #mse in 61 years in 4 OS locations
  mseS<-sum(resS**2, na.rm=TRUE)
  nS<-sum(!is.na(resS))
  #mse in 61 years in 36 S locations
 
  return(list(msee=c(mseOS, mseS),nn=c(nOS, nS), M=M))

}

###Application to non-spatial  models

#Model MM0 in Fit_Models.R
formula<-(IncS~1)
msumM0<-matrix(NA, nrow=2,ncol=10)  #10
nnM0<-matrix(NA, nrow=2,ncol=10) 
for(j in (1:10))  # 10
{ 
  cat("Step: ", j, fill=T)
  aux<-Sel_Rmse_OOS(formula=formula,  ii=mii[,j])
  msumM0[,j]<-aux$msee
  nnM0[,j]<-aux$nn
}

mmseM0<-msumM0/nnM0
sumM0<-apply(msumM0, MARGIN=1, FUN=sum)
nM0<-apply(nnM0, MARGIN=1, FUN=sum)
rmseM0<-(sumM0/nM0)**(0.5)
names(rmseM0)<-c("rmseOS", "rmseS")
rmseM0

#Model MM1 in Fit_Models.R
formula<-(IncS~1+Cov)
msumM1<-matrix(NA, nrow=2,ncol=10)  
nnM1<-matrix(NA, nrow=2,ncol=10) 
for(j in (1:10))  
{ 
  cat("Step: ", j, fill=T)
  aux<-Sel_Rmse_OOS(formula=formula,  ii=mii[,j])
  msumM1[,j]<-aux$msee
  nnM1[,j]<-aux$nn
}

mmseM1<-msumM1/nnM1
sumM1<-apply(msumM1, MARGIN=1, FUN=sum)
nM1<-apply(nnM1, MARGIN=1, FUN=sum)
rmseM1<-(sumM1/nM1)**(0.5)
names(rmseM1)<-c("rmseOS", "rmseS")
rmseM1



### Selecting spatial models

Sel_Rmse_OOS_S<-function(formula, ii,  ddatosvC2=datosvC2,  ttheta=NULL,
                 ss.index=s.index1, sspdepc=spdepc, mmesh=mesh,ccoords=coords)
{  
  IncSC<-ddatosvC2$IncS
  ind<-((ddatosvC2$siteInd==ii[1])|(ddatosvC2$siteInd==ii[2])|(ddatosvC2$siteInd==ii[3])|(ddatosvC2$siteInd==ii[4]) )
	#ind is indicator out of sample observations. length: 61*92*40
  indS<-(ind==FALSE)
	#indS is indicator of sample observations

  est.coords<-ccoords[indS,]
  val.coords<-ccoords[ind,]
  est.data<-ddatosvC2$IncS[indS]
  val.data<-ddatosvC2$IncS[ind]

  A.est<-inla.spde.make.A(mmesh, loc=est.coords)
  A.val<-inla.spde.make.A(mmesh, loc=val.coords)

  stack.est<-inla.stack(data=list(IncS=est.data), tag="est",
                      A=list(A.est,1),
                      effects=list(c(ss.index, list(Intercept=1)), list(Cov=ddatosvC2$Cov[indS,], 
                                                                        yearInd=ddatosvC2$yearInd[indS],dayInd=ddatosvC2$dayInd[indS],
                                                                        timeInd=ddatosvC2$timeInd[indS],siteInd=ddatosvC2$siteInd[indS]))
                   )
  
  stack.val<-inla.stack(data=list(IncS=NA), tag="val",
                      A=list(A.val,1),
                      effects=list(c(ss.index, list(Intercept=1)), list(Cov=ddatosvC2$Cov[ind,], #CovS=ddatosvC2$CovS[ind,],
                                                                        yearInd=ddatosvC2$yearInd[ind],dayInd=ddatosvC2$dayInd[ind],
                                                                        timeInd=ddatosvC2$timeInd[ind],siteInd=ddatosvC2$siteInd[ind]))
                   )
  join.stack<-inla.stack(stack.est, stack.val)
  index.val<-inla.stack.index(join.stack,"val")$data
  index.est<-inla.stack.index(join.stack,"est")$data
  
  M<-inla(formula, data=inla.stack.data(join.stack,spde=sspdepc), 
            family="gamma", 
            control.compute = list(dic = FALSE,cpo=FALSE,waic=FALSE),
            control.predictor=list(A=inla.stack.A(join.stack),compute=T,link=1),
            control.family = list(control.link=list(model="log")),
            control.mode = list(theta =ttheta,restart=FALSE)
  )
  
	resS<-(IncSC[indS]-M$summary.fitted.values$mean[index.est])/10
	yearS<-ddatosvC2$yearInd[indS] #year of the residuals S
	resOS<-(IncSC[ind]-M$summary.fitted.values$mean[index.val])/10
	yearOS<-ddatosvC2$yearInd[ind] #year of the residuals OS
	mseOS<-sum(resOS**2, na.rm=TRUE)
	nOS<-sum(!is.na(resOS))
#mse in 61 years in 4 OS locations
	mseS<-sum(resS**2, na.rm=TRUE)
	nS<-sum(!is.na(resS))
#mse in 61 years in 36 S locations
	return(list(msee=c(mseOS, mseS),nn=c(nOS, nS), M=M))#, indes.est=index.est))
}


###Aplication to spatial models

bnd<-inla.nonconvex.hull(Ucoords)
mesh<-inla.mesh.2d(boundary=bnd,max.edge=c(1, 4))
range0<-min(diff(range(coords[,1])),diff(range(coords[,2])))/2
spdepc<-inla.spde2.pcmatern(mesh=mesh, prior.range=c(range0,0.5),prior.sigma=c(1,0.1))
s.index1<-inla.spde.make.index (name="spatialInd", n.spde=spdepc$n.spde)

#Model MMS0 in Fit_Models
formula<-(IncS~-1+Intercept+f(spatialInd, model=spdepc))
msumMS0<-matrix(NA, nrow=2,ncol=10)  
nnMS0<-matrix(NA, nrow=2,ncol=10) 
for(j in (1:10))  
{ 
  cat("Step: ", j, fill=T)
  aux<-Sel_Rmse_OOS_S(formula=formula,  ii=mii[,j])
  msumMS0[,j]<-aux$msee
  nnMS0[,j]<-aux$nn
}  
mmseMS0<-msumMS0/nnMS0
sumMS0<-apply(msumMS0, MARGIN=1, FUN=sum)
nMS0<-apply(nnMS0, MARGIN=1, FUN=sum)
rmseMS0<-(sumMS0/nMS0)**(0.5)
names(rmseMS0)<-c("rmseOS", "rmseS")
rmseMS0
  
#Model MMS1 in Fit_Models
formula<-(IncS~-1+Intercept+Cov+f(spatialInd, model=spdepc))
msumMS1<-matrix(NA, nrow=2,ncol=10)  #10
nnMS1<-matrix(NA, nrow=2,ncol=10) 
for(j in (1:10))  # 10
{ 
   cat("Step: ", j, fill=T)
   aux<-Sel_Rmse_OOS_S(formula=formula,  ii=mii[,j])
   msumMS1[,j]<-aux$msee
   nnMS1[,j]<-aux$nn
}

mmseMS1<-msumMS1/nnMS1
sumMS1<-apply(msumMS1, MARGIN=1, FUN=sum)
nMS1<-apply(nnMS1, MARGIN=1, FUN=sum)
rmseMS1<-(sumMS1/nMS1)**(0.5)
names(rmseMS1)<-c("mseOS",  "mseS")
rmseMS1

#Model MMS3 in Fit_Models
formula<-(IncS~-1+Intercept+Cov+f(spatialInd, model=spdepc)+f(yearInd, model="iid"))
msumMS2<-matrix(NA, nrow=2,ncol=10)  #10
nnMS2<-matrix(NA, nrow=2,ncol=10) 
for(j in (1:10))  # 10
{ 
   cat("Step: ", j, fill=T)
   aux<-Sel_Rmse_OOS_S(formula=formula,  ii=mii[,j])
   msumMS2[,j]<-aux$msee
   nnMS2[,j]<-aux$nn
}

mmseMS2<-msumMS2/nnMS2
mmseMS2
sumMS2<-apply(msumMS2, MARGIN=1, FUN=sum)
nMS2<-apply(nnMS2, MARGIN=1, FUN=sum)
rmseMS2<-(sumMS2/nMS2)**(0.5)
names(rmseMS2)<-c("rmseOS", "rmseS")
rmseMS2

#Model MMS7ar in Fit_Models
formula<-(IncS~-1+Intercept+Cov+f(spatialInd, model=spdepc)+
            f(dayInd, model="ar1", replicate=yearInd))
msumMS3<-matrix(NA, nrow=2,ncol=10)  
nnMS3<-matrix(NA, nrow=2,ncol=10) 
theta<-c( 0.6050486,0.5850171,-1.5598324,2.6399777,1.3907978 )
  #c(MMS7ar$mode$theta[1], MMS7arexp$mode$theta)
for(j in (1:10))  # 10
{ 
   cat("Step: ", j, fill=T)
   aux<-Sel_Rmse_OOS_S(formula=formula,  ii=mii[,j],ttheta=theta)
   msumMS3[,j]<-aux$msee
   nnMS3[,j]<-aux$nn
}

mmseMS3<-msumMS3/nnMS3
mmseMS3
sumMS3<-apply(msumMS3, MARGIN=1, FUN=sum)
nMS3<-apply(nnMS3, MARGIN=1, FUN=sum)
rmseMS3<-(sumMS3/nMS3)**(0.5)
names(rmseMS3)<-c("rmseOS",  "rmseS")
rmseMS3






