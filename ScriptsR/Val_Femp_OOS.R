

library(INLA)

###Functions

# a is initial year of the period to compute the empirical distribution
#  b is final year of the period to compute the empirical distribution

Val_Femp_OOS<-function(ii, ddatosvC2=datosvC2,ccoords=coords, mmesh=mesh, 
	sspdepc=spdepc, tthetaMSel=thetaMSel,mseed=546, 
	nsamples=1000,nnS=nS,nnT=nT,nnL=nL,a=33, b=62, ddistance=distance)
{
	ddatosvC2$IncS<-ddatosvC2$IncS/10

	indOS<-(ddatosvC2$siteInd==ii)
	#indOS is indicator out of sample observations. length: 63*92*40=231840, sum: 63*92= 5796  
	indS<-(indOS==FALSE)
	#indS is indicator of sample observations

	###Estimation  of the model without location ii

	est.coords<-ccoords[indS,]
	val.coords<-ccoords[indOS,]
	est.data<-ddatosvC2$IncS[indS]
	val.data<-ddatosvC2$IncS[indOS]
  
	A.est<-inla.spde.make.A(mmesh, loc=est.coords)
	A.val<-inla.spde.make.A(mmesh, loc=val.coords)
	s.index<-inla.spde.make.index (name="spatialInd", n.spde=sspdepc$n.spde)
  
	stack.est<-inla.stack(data=list(IncS=est.data), tag="est",
                        A=list(A.est,1),
                        effects=list(c(s.index, list(Intercept=1)), list(Cov=ddatosvC2$Cov[indS,],
                              yearInd=ddatosvC2$yearInd[indS],dayInd=ddatosvC2$dayInd[indS],
                              timeInd=ddatosvC2$timeInd[indS],siteInd=ddatosvC2$siteInd[indS]))
	)
	stack.val<-inla.stack(data=list(IncS=NA), tag="val",
                        A=list(A.val,1),
                        effects=list(c(s.index, list(Intercept=1)), list(Cov=ddatosvC2$Cov[indOS,], 
                             yearInd=ddatosvC2$yearInd[indOS],dayInd=ddatosvC2$dayInd[indOS],
                             timeInd=ddatosvC2$timeInd[indOS],siteInd=ddatosvC2$siteInd[indOS]))
	)
	join.stack<-inla.stack(stack.est, stack.val)
	index.val<-inla.stack.index(join.stack,"val")$data
	index.est<-inla.stack.index(join.stack,"est")$data

	formula<-(IncS~-1+Intercept+Cov+f(spatialInd, model=spdepc)+
		f(dayInd, model="ar1", replicate=yearInd))
	M<-inla(formula, data=inla.stack.data(join.stack,spde=sspdepc), 
          family="gamma", 
          control.compute = list(config = TRUE),
          control.predictor=list(A=inla.stack.A(join.stack),compute=T,link=1),
          control.family = list(control.link=list(model="log")),
          control.mode = list(theta =tthetaMSel,restart=FALSE)
	)
  
	###Simulation of samples in ii using the fitted model (out of sample)

	set.seed(mseed)
	sampleOS<-inla.posterior.sample(nsamples, M)
	indiaux<-c(1:(nnS*nnT*nnL))
	indiOS<-indiaux[indOS] 
	aux<-lapply(sampleOS, function(x) x$latent[indiOS])
	sampPredOS<-matrix(unlist(aux), byrow=T, nrow=length(aux))
	colnames(sampPredOS)<-rownames(sampleOS[[1]]$latent)[indiOS]
#	print(dim(sampPredOS))#  dimension  must be nsamples x 92*61=5612
	sampPhiOS<-unlist(lapply(sampleOS, function(x) x$hyperpar[1]))
	sampScaleOS<-exp(sampPredOS)/sampPhiOS
	sampYOS<-t(apply(cbind(sampScaleOS, sampPhiOS), MARGIN=1, FUN=genaux))
	rownames(sampYOS)<-NULL
#	print(dim(sampYOS))# dimension must be nsamples x 92*61

	###Computation of the empirical distribution using the simulated  and the observed samples  and plot them

	#Empirical function with observed data
	Femp<-EmpF(IIncS=ddatosvC2$IncS[indOS], indNA=(!is.na(ddatosvC2$IncS[indOS])), 
        indgroup=ddatosvC2$yearInd[indOS], nR1=a, nR2=b, ddistance=ddistance)	
	#Empirical function with simulated  data
	lims<-length(Femp)-1
	pc=seq(0,lims, by=ddistance)

	mFsim<-apply(sampYOS, MARGIN=1, FUN=EmpF, pc=pc, nR1=a, nR2=b, 
        ddistance=ddistance, indNA=(!is.na(ddatosvC2$IncS[indOS])),indgroup=ddatosvC2$yearInd[indOS])
	
	Fsim2mean<-apply(mFsim, MARGIN=1, FUN=mean)
	Fsim2p025<-apply(mFsim, MARGIN=1, FUN=quantile, prob=0.025)
	Fsim2p975<-apply(mFsim, MARGIN=1, FUN=quantile, prob=0.975)
	return(list(Femp=Femp, FsimSummary=cbind(Fsim2mean,Fsim2p025,Fsim2p975),pc=pc))
}


EmpF<-function(IIncS, indNA,indgroup, nR1=2, nR2=8, ddistance=1,pc=NULL)
{
  iIncO<-IIncS[indNA]
  iindgroup<-indgroup[indNA]
  Rpos<-iIncO[(iindgroup<=nR2)&(iindgroup>=nR1)]

  #sample of records number nR in location 
   if (is.null(pc)) pc<-seq(0, (max(Rpos, na.rm=T)+ddistance), by=ddistance)
  estFO<-estF(x=Rpos,pcorte=pc)
	if (!is.null(pc)) return(estFO[,2]) else  return(estFO)
}

estF<-function(x,pcorte=pc)
{	
	n<-length(x)
	Fpc<-c(0,cumsum(table(cut(x,breaks=pcorte))))/n
	return(cbind(pcorte, Fpc))
}

genaux<-function(param)
{
	m<-length(param)-1
	vec<-rgamma(m, shape=param[(m+1)], scale=param[1:m])
	return(vec)
}


### Use of the function to  obtain plots in Figure 4

#Input data: Datos_Input.RData

bnd<-inla.nonconvex.hull(Ucoords)
mesh<-inla.mesh.2d(boundary=bnd,max.edge=c(1, 4))
range0<-min(diff(range(coords[,1])),diff(range(coords[,2])))/2
spdepc<-inla.spde2.pcmatern(mesh=mesh, prior.range=c(range0,0.5),prior.sigma=c(1,0.1))

thetaMSel<-c(0.6050486,0.5850447,-1.5617280,2.2620602, 1.3165265 ) 
#MS7ar$mode$theta   

obj<-Val_Femp_OOS(ii=18, ddatosvC2=datosvC2,mseed=546, 
	nsamples=1000,a=33, b=62, ddistance=1) #Bilbao
obj<-Val_Femp_OOS(ii=40, ddatosvC2=datosvC2,mseed=546, 
	nsamples=1000,a=33, b=62, ddistance=1) #Daroca
obj<-Val_Femp_OOS(ii=27, ddatosvC2=datosvC2,mseed=546, 
	nsamples=1000,a=33, b=62, ddistance=1) #Vitoria
obj<-Val_Femp_OOS(ii=33, ddatosvC2=datosvC2,mseed=546, 
	nsamples=1000,a=33, b=62, ddistance=1) #Huelva


titu<-"Bilbao" #to be changed
ccex<-1 
#pdf(file = "G:/Mi unidad/Espacial/E9RecordsExpl/Graphs/CumPBilbaoFinal.pdf", width = 4,height = 4)
plot(obj$pc,obj$Femp, pch=17,xlab="x", ylab="Cum. probability",cex=ccex,
	cex.axis=ccex, cex.lab=ccex)
points(obj$pc,obj$FsimSummary[,1],  col="red", pch=16, cex=ccex)
segments(obj$pc,obj$FsimSummary[,2], obj$pc,obj$FsimSummary[,3], col="red")
title(titu, cex.main=ccex)
legend("bottomright", bty="n", pch=c(16, 17), col=c("red", "black"), 
		legend=c("Model", "Observed"),cex=ccex)
#dev.off()




