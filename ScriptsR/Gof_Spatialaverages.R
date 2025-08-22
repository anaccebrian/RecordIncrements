
### Functions

Gof_Spatialaverages<-function(index,fitted, Inc, ssampY,
				p=0.95, type=" ", vcex=1, vcex2=1.7, yylab=expression(bar(J)[t]),
				xaxis=NULL, marpar=c(5, 6, 4, 2))
{
	indNA=is.na(Inc)
	SMIncMaux<-apply(ssampY,MARGIN=1,FUN=IMY, iindex=index,iindNA=indNA)
	#SMIncMaux is a matrix of dimension the number of levels of index by the number of simulations
	SIncMaux<-apply(SMIncMaux,MARGIN=1,FUN=mean)
	SIncPIaux<-apply(SMIncMaux,MARGIN=1,FUN=quantile, probs=0.025, na.rm=TRUE)
	SIncPSaux<-apply(SMIncMaux,MARGIN=1,FUN=quantile, probs=0.975, na.rm=TRUE)

	#computes the means by grupos defined in  index of observed increments
	IncMaux<-tapply(Inc,INDEX=index,FUN=mean, na.rm=TRUE)

 	par(mar = marpar)
	if (is.null(xaxis)) xaxix<-unique(index)
	plot(xaxis, IncMaux,  pch=17, col="black", ylab=yylab,
		ylim=c(min(SIncPIaux, na.rm=T),max(SIncPSaux, na.rm=T)),xlab=type, 
		cex=vcex, cex.lab=vcex2,cex.axis=vcex2)
	points(xaxis, SIncMaux, cex=vcex, pch=16, col="red")
	segments(xaxis, SIncPIaux,xaxis, SIncPSaux, col="red")
	legend("topright", col=c("red","black"), pch=c(16,17), legend=c("Model", "Empirical"),cex=vcex2,bty="n")

	#coverage in sample by the levels defined by index
	coverage<-mean((IncMaux>=SIncPIaux)&(IncMaux<=SIncPSaux),na.rm=T)
	cat("Coverage: ", round(coverage,3), fill=T)

	return(coverage)
}

# IMY computes the means by grupos defined in  iindex of x, usually a simulated series of increments
IMY<-function(x,iindex, iindNA)
{
	x[iindNA]<-NA # only  simulated increments at observed positions are considered
	meangroup<-tapply(x,INDEX=iindex,FUN=mean, na.rm=TRUE)
	return(meangroup)
}

#Functions to generate samples

genParam<-function(M, seed=4123, nsim=1000)
{
	inla.setOption(inla.mode = "classic")
	Mm<-inla.hyperpar(M)
	set.seed(seed)
	sample<-inla.posterior.sample(nsim, Mm)
	return(sample)
}

genY<-function(sample, T=nT, iind=ind,seed=4563)
{
	#sample is a list where each element is a list corresponding to one sample, with elements latent and others
	#sample[[1]]$latent is a matrix with one column that is one simulatet observation for each latent field
	aux<-lapply(sample, function(x) x$latent[ind])
	sampPred<-matrix(unlist(aux), byrow=T, nrow=length(aux))
	colnames(sampPred)<-rownames(sample[[1]]$latent)[iind]
	#dim(sampPred)
	#sampPred is a matrix with number of rows equal to the number o simulations and each colum the fitted value at a point (time-space)
	#summary(sampPred[,1])

	#sample of Phi  and scale parameters of Gamma distribution for each response in the model
	sampPhi<-unlist(lapply(sample, function(x) x$hyperpar[1]))
	sampScale<-exp(sampPred)/sampPhi
	set.seed(seed)
	sampY<-t(apply(cbind(sampScale, sampPhi), MARGIN=1, FUN=genaux))
	rownames(sampY)<-NULL
	#sampY must be a matrix with n=1000 rows (sample size) and 224480 columns, one simulated increment
	# for each point at time and location
	return(sampY)
}

genaux<-function(param)
{
	#genaux generate a Gamma value for each parameter
	#param is a vctor containing  m scale parameters  and one common shape parameter
	m<-length(param)-1
	vec<-rgamma(m, shape=param[(m+1)], scale=param[1:m])
	return(vec)
}




###  Application  of Gof_Spatialaverages for Figure 5


library(INLA)

#generation of  samples of the incrments at observed positions
sampleP<-genParam(M=MSelT)
#sampleP is a list where each element is a list corresponding to one sample, with elements latent and others
ind<-c(1:(nS*nT*nL))  #40*92*61=224480
sampY<-genY(sample=sampleP,iind=ind)


fittedM<-MSelT$summary.fitted.values$mean[1:(nS*nT*nL)]

#pdf("G:/Mi unidad/Espacial/E9RecordsExpl/Graphs/PredYears.pdf")
nnyears<-Gof_Spatialaverages(index=datosvC2$yearInd, fitted=fittedM,
			Inc=datosvC2$IncS/10, ssampY=sampY/10, type="Years",xaxis=c(1961:2021))
#dev.off()

#pdf("G:/Mi unidad/Espacial/E9RecordsExpl/Graphs/PredDays.pdf")
nndays<-Gof_Spatialaverages(index=datosvC2$dayInd,fitted=fittedM,
 	Inc=datosvC2$IncS/10,ssampY=sampY/10, type="Days",xaxis=c(1:92))
#dev.off()

save(sampleP, sampY, file="C:/Users/PC/Desktop/AnaE9RecordsExpl/DatosSim/Sim_Samples_Obs.RData")

