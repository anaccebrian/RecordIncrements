

Val_Coverage_OOS<-function(ii, fformula=formula, tthetaMSel=thetaMSel, nsamples=1000,
		ddatosvC2=datosvC2,ccoords=coords, mmesh=mesh, sspdepc=spdepc,
		mseed=9541)
{

  obs=datosvC2$IncS
  indOS<-(ddatosvC2$siteInd%in%ii)
  #indOS is indicator out of sample observations. length: 61*92*40=224480, sum: 63*92*4= 22448 
  indS<-(indOS==FALSE)
  #indS is indicator of sample observations

  #Estimation  of the model without location ii
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
  M<-inla(fformula, data=inla.stack.data(join.stack,spde=sspdepc), 
          family="gamma", 
          control.compute = list(config = TRUE),
          control.predictor=list(A=inla.stack.A(join.stack),compute=T,link=1),
          control.family = list(control.link=list(model="log")),
          control.mode = list(theta =thetaMSel,restart=FALSE)
  )
  
  set.seed(mseed)
  sampleOS<-inla.posterior.sample(nsamples, M)
  indiaux<-c(1:(nS*nT*nL))
  indiOS<-indiaux[indOS] 
  aux<-lapply(sampleOS, function(x) x$latent[indiOS])
  sampPredOS<-matrix(unlist(aux), byrow=T, nrow=length(aux))
  colnames(sampPredOS)<-rownames(sampleOS[[1]]$latent)[indiOS]
 # print(dim(sampPredOS))#  dimension  must be nsamples x 92*61= 5612

  sampPhiOS<-unlist(lapply(sampleOS, function(x) x$hyperpar[1]))
  sampScaleOS<-exp(sampPredOS)/sampPhiOS
  sampYOS<-t(apply(cbind(sampScaleOS, sampPhiOS), MARGIN=1, FUN=genaux))
  rownames(sampYOS)<-NULL
 # print(dim(sampYOS))#  dimension  must be nsamples x 92*61= 5612
	
  #comparison of observed and simulated IC
  obsOS<-obs[indOS]
  Ii<-apply(sampYOS, MARGIN=2, FUN=quantile, probs=0.025, na.rm=T)
  Is<-apply(sampYOS, MARGIN=2, FUN=quantile, probs=0.975, na.rm=T)
  pert<-((Ii<=obsOS)&(obsOS<=Is))
#	cat(" length Ii, obs Is: ", length(Ii), length(Is), length(obs),length(pert), fill=T) 
#	cat(" length Ii, obs Is: ", sum(!is.na(Ii)), sum(!is.na(Is)), sum(!is.na(obs)),sum(!is.na(pert)), fill=T) 
  return(list(pert=pert, indOS=indOS))
}

genaux<-function(param)
{
	m<-length(param)-1
	set.seed(sum(param))
	vec<-rgamma(m, shape=param[(m+1)], scale=param[1:m])
	return(vec)
}

# Application

thetaMSel<-c(0.6052908,0.5904362,-1.5606298,2.2622344,1.3158043) # MMS7ar$mode$theta
bnd<-inla.nonconvex.hull(Ucoords)
mesh<-inla.mesh.2d(boundary=bnd,max.edge=c(1, 4))
range0<-min(diff(range(coords[,1])),diff(range(coords[,2])))/2
spdepc<-inla.spde2.pcmatern(mesh=mesh, prior.range=c(range0,0.5),prior.sigma=c(1,0.1))

formula<-(IncS~-1+Intercept+Cov+f(spatialInd, model=spdepc)+
            f(dayInd, model="ar1", replicate=yearInd)) 

pertT<-rep(NA, nT*nL*nS)

mii<-cbind(c(8,28,29,34),c(6,18,22,40),c(7,15,26,38),c(13,16,35,39),c(2,9,10,23),
           c(3,11,19,27),c(5,12,30,37),c(1,4,17,24),c(25,31,33,36),c(14,20,21,32)  )

for (ii in c(1:10))
{
	print(ii)
	aux<-Val_Coverage_OOS(mii[,ii]) #,mseed=1981)
	pertT[aux$indOS]<-aux$pert
 }

mean(pertT, na.rm=T)


##Plots of the coverage by year, by day within the year, by site
Iyear<-datosvC2$yearInd
covyear<-tapply(pertT, FUN=mean,INDEX=Iyear, na.rm=T)
sort(covyear)
contar<-function(x){	return(nn=sum(is.na(x)==FALSE)) }
tapply(pertT, FUN=contar,INDEX=Iyear)
#pdf(file="G:/Mi unidad/Espacial/E9RecordsExpl/articuloSERRA/Enviado2mayo/RecordApplication/Graphs/CovYear.pdf")
plot(c(1961:2021), covyear, pch=16, xlab="Year", ylab="Average coverage", ylim=c(0,1),
	cex=1.8, cex.axis=1.8,cex.lab=1.8)
lines(c(1961:2021), covyear)
abline(h=0.95, col="gray")
#dev.off()

covday<-tapply(pertT, FUN=mean,INDEX=datosvC2$dayInd, na.rm=T)
sort(covday)
tapply(pertT, FUN=contar,INDEX=datosvC2$dayInd)
#pdf(file="G:/Mi unidad/Espacial/E9RecordsExpl/articuloSERRA/Enviado2mayo/RecordApplication/Graphs/CovDay.pdf")
plot(c(1:92), covday, pch=16, xlab="Day", ylab="Average coverage", ylim=c(0,1),
	cex=1.8, cex.axis=1.8,cex.lab=1.8)
lines(c(1:92), covday)
abline(v=c(30,61), col="gray", lty=2)
abline(v=c(30,61), col="gray", lty=2)
abline(h=0.95, col="gray")
#dev.off()

sitios<-c(1:40)
covsite<-tapply(pertT, FUN=mean,INDEX=datosvC2$siteInd, na.rm=T)
sort(covsite)
tapply(pertT, FUN=contar,INDEX=datosvC2$siteInd)
#pdf(file="G:/Mi unidad/Espacial/E9RecordsExpl/articuloSERRA/Enviado2mayo/RecordApplication/Graphs/CovSites.pdf")
plot(sitios, covsite, pch=16, xlab="Sites", ylab="Average coverage", ylim=c(0,1),
cex=1.8, cex.axis=1.8,cex.lab=1.8)
lines(sitios, covsite)
abline(h=0.95, col="gray")
#dev.off()


save(pertT, covyear, covday, covsite, file="G:/Mi unidad/Espacial/E9RecordsExpl/Datos/Coverages.RData")


