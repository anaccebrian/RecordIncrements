

###Functions

#Simulate samples of the increments using the conditional
#posterior distribution given that the previous day was a record or not



Sim_Cond<-function(Ilag, CCovGrid=CovGrid, nnP=nP,nnT=nT, nnL=nL,
			ccoordsGrid=coordsGrid,mmesh=mesh,
			nSim=500,MSel=MSelT,mseed=NULL)
{
	CovGridC<-cbind(CCovGrid[,1:2], Ilag, CCovGrid[,3:5],
		    CCovGrid[,4]*Ilag, CCovGrid[,6])

	#Defining projection matrix for defining the spatial random effects
	#Create Ap.s: matriz con filas: numero nP*nT*nL y columnas G
	Ap.s<-inla.spde.make.A(mmesh, loc=ccoordsGrid)

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
	#  samples is a list with nSim elements, one for each simulation
	# each list has 3 elelments:names(samples[[1]]): "hyperpar" "latent"   "logdens" 
	xnames <- rownames(samples[[1]]$latent) ### collect the names
	idx <- lapply(c("Intercept","Cov","spatialInd","dayInd"), function(nam) 
	              which(substr(xnames, 1, nchar(nam))==nam))
      #idx is a list with 4 elements each of  different length: 1,8,363, 5796
      # that give the indexes where those elements are in samples[[1]]$latent

      mat.samples1 <- sapply(samples, function(spl)
	                c(Intercept=spl$latent[idx[[1]]], Cov=spl$latent[idx[[2]]], 
	                spatialInd=spl$latent[idx[[3]]]))
	  #dim(mat.samples1) #  372(=9+363)  x nSim
      mat.samplesaux <- sapply(samples, function(spl)	dayInd=spl$latent[idx[[4]]])
      mat.samples2<-mat.samplesaux[(nnL+1):dim(mat.samplesaux)[1],] 
	  #nnL+1 because the first one must be removed since  indYear  goes from 2 to nT but generates from 1 to nT
      dim(mat.samples2) # 5612  x nSim  
		rm(mat.samplesaux,xnames,idx)

      mu.samp<-cbind(Intercept=1, Cov=CovGridC,s=Ap.s)%*%mat.samples1
      #dim(mu.samp) #  nnL*nnT*nnP x nSim  92*61*844  x 500
	  lmu.samp<-lapply(seq_len(nSim), function(i,mm) mm[,i], mm=mu.samp)
	  rm(mu.samp)
	  lsamp2<-lapply(seq_len(nSim), function(i,mm) mm[,i], mm=mat.samples2)

	  for (j in c(1:nSim))
	  {
		 print(j)
		 for (i in seq(1,nnP*nnL*nnT, by=(nnL*nnT)))
		 {
		 	lmu.samp[[j]][i:(i+nnL*nnT-1)]<-exp(lmu.samp[[j]][i:(i+nnL*nnT-1)]+lsamp2[[j]])
		 }
	  }
      rm(mat.samples2, mat.samples1, lsamp2)

	  sampPhiP<-unlist(lapply(samples, function(x) x$hyperpar[1]))
	  set.seed(mseed)
	  sampYP<-lapply(c(1:nSim), FUN=genaux2, llmu.samp=lmu.samp, ssampPhiP=sampPhiP,
					mm=nnL*nnT*nnP)
	  #sampYP is a list with nSim elements, each one a simulation
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
load("G:\\Mi unidad\\Espacial\\E9RecordsExpl\\Datos_Z\\Selected_Model.RData") 


sampYP1<-Sim_Cond(Ilag=1, mseed=2341)
sampYP0<-Sim_Cond(Ilag=0, mseed=7651)


save(sampYP0, sampYP1, file="C:/Users/PC/Desktop/AnaE9RecordsExpl/DatosSim/SimDataCond.RData")


