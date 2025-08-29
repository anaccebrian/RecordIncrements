



###Defining covariates at  the prediction grid

#debe ser una matriz con filas nP*nL*(nT-2) y columnas 8:
#  "logt"  "logt2"  "Ixlag" "logAlt"  "logCoast" "logtlogCoast"  "IxlaglogCoast" "IxlagSp"  
#distP,elevP cts en el tiempo, varian en el espacio
#IxlagSp es cte espacialmente, varia en el tiempo cada 92*63 

Sim_Marg<-function(CCovGrid=CovGrid, IIsimlag=Isimlag, nnP=nP,nnT=nT, nnL=nL,
			ccoordsGrid=coordsGrid,mmesh=mesh,
			nSim=500,MSel=MSelT,mseed=NULL)
{
	#Defining projection matrix for defining the spatial random effects
	#Create Ap.s: matriz con filas: numero nP*nT*nL y columnas G
	Ap.s<-inla.spde.make.A(mmesh, loc=ccoordsGrid)

	#simulating the model parameters
	if (is.null(mseed)==F)
	{
		inla.seed = as.integer(runif(1)*.Machine$integer.max)
  		set.seed(mseed)
		samples<-inla.posterior.sample(n=nSim, result=MSel, add.names=FALSE,
			seed = inla.seed, num.threads="1:1")
	}else
	{
		samples<-inla.posterior.sample(n=nSim, result=MSel, add.names=FALSE)
	}
	# samples is a list with nSim elements, one for each simulation
	# each list has 3 elelments:names(samples[[1]]): "hyperpar" "latent"   "logdens" 
	# In each simulation samples[[1]]$latent contains the simulated values of
	# Apredictor (n=224480=61*92*40), Predictor (224564) spatialInd (spatial re,G=363)
	# dayInd ( 5704=92*62)  index t for records goes from 2 to 62 and generates from 1 to 62
	# b0(1) Covariates (8)

	xnames <- rownames(samples[[1]]$latent) ### collect the names
	idx <- lapply(c("Intercept","Cov","spatialInd","dayInd"), function(nam) 
	              which(substr(xnames, 1, nchar(nam))==nam))
	# idx is a list with 4 elements each of  different length: 1,8,363, 5796
	# that gives the indexes where those elements are in samples[[1]]$latent

      mat.samples1 <- sapply(samples, function(spl)
	                c(Intercept=spl$latent[idx[[1]]], Cov=spl$latent[idx[[2]]], 
	                spatialInd=spl$latent[idx[[3]]]))
	#dim(mat.samples1) #  372(=9+363)  x nSim
	mat.samples1a<-mat.samples1[-c(4,8),] # separate the samples of coefficients of Ilag and Ilag interaction
	Icoef<-mat.samples1[4,] #length(Icoef)=nSim
	Icoefint<-mat.samples1[8,]

	mat.samplesaux <- sapply(samples, function(spl)	dayInd=spl$latent[idx[[4]]])
	mat.samples2<-mat.samplesaux[(nnL+1):dim(mat.samplesaux)[1],] 
	# starting in (nL+1) to remove  obs first year since  samples are generated for
	#  1 to 62  but indYear goes from 2 a 62
	rm(mat.samplesaux,xnames,idx)

	mu.sampb1<-sapply(1:nSim, FUN=function(i, IIIsimlag,IIcoef)
			{return(as.vector(IIIsimlag[i,,,])*IIcoef[i])}, 
		    IIIsimlag=IIsimlag, IIcoef=Icoef)
	# dim(mu.sampb1)   nL*nT*nP x nSim    4736528=844*61*92 X 500
	mu.sampb2<-sapply(1:nSim, FUN=function(i, IIIsimlag,IIcoef, ccov)
			{return(as.vector(IIIsimlag[i,,,])*ccov*IIcoef[i])}, 
		    IIIsimlag=IIsimlag, IIcoef=Icoefint, ccov=CCovGrid[,4])
	# dim(mu.sampb2)   nL*nT*nP x nSim    4736528=844*61*92 X 500 
	mu.sampa<-cbind(Intercept=1, Cov=CCovGrid,s=Ap.s)%*%mat.samples1a
	# dim(mu.sampa)   nL*nT*nP x nSim    4736528=844*61*92 X 500
	mu.samp<-mu.sampa+mu.sampb1+mu.sampb2
	#warning: sparse-dense method
	rm(mu.sampa, mu.sampb1, mu.sampb2)

	lmu.samp<-lapply(seq_len(nSim), function(i,mm) mm[,i], mm=mu.samp)
	rm(mu.samp)
	lsamp2<-lapply(seq_len(nSim), function(i,mm) mm[,i], mm=mat.samples2)
	for (j in c(1:nSim))
	{
		print(j)
		for (i in seq(1,nnP*nnL*nnT, by=nnL*nnT))
		{
			lmu.samp[[j]][i:(i+nnL*nnT-1)]<-exp(lmu.samp[[j]][i:(i+nnL*nnT-1)]+lsamp2[[j]])
		}
	}
	rm(mat.samples2, mat.samples1, mat.samples1a)
	rm(lsamp2)

	sampPhiP<-unlist(lapply(samples, function(x) x$hyperpar[1]))
	set.seed(mseed)
	sampYP<-lapply(c(1:nSim), FUN=genaux2, llmu.samp=lmu.samp, ssampPhiP=sampPhiP,mm=nL*nT*nP)
	# sampYP is a list with nSim elements, each one a simulation
	rm(lmu.samp)
	return(sampYP)
}


genaux2<-function(i,llmu.samp, ssampPhiP, mm)
{
	ssampScaleP<-llmu.samp[[i]]/ssampPhiP[i]
	vec<-rgamma(mm, shape=ssampPhiP[i], scale=ssampScaleP)/10		
	# sample generated in ºC
	return(vec)
}

###Application of the functions

library(INLA)
load("G:\\Mi unidad\\Espacial\\E9RecordsExpl\\Datos_Z\\Grid_Spat_Cov.RData")
load("G:\\Mi unidad\\Espacial\\E9RecordsExpl\\Datos_Z\\Grid_Spat_Isim.RData")
load("G:\\Mi unidad\\Espacial\\E9RecordsExpl\\Datos_Z\\Selected_Model.RData") 

sampYP<-Sim_Marg(IIsimlag=Isimlag, mseed=2341)


#save(sampYP,file="C:/Users/PC/Desktop/AnaE9RecordsExpl/DatosSim/Sim_Data_Marg.RData")



