

#############################################################
###Functions for computing  posterior summaries of the marginal 
###distribution of the increments and their use
#############################################################


### Loading data
load("G:\\Mi unidad\\Espacial\\E9RecordsExpl\\Datos_Z\\Grid_Spat_Cov.RData")
load("C:\\Users\\PC\\Desktop\\AnaE9RecordsExpl\\DatosSim\\Sim_Data_Marg.RData")
source("G:\\Mi unidad\\Espacial\\E9RecordsExpl\\Ficheros_Z\\Function_MapSpatial.R")





### Posterior distribution of the average increments over two periods
### at a spatial grid

Marg_Dist_Averages<-function(ssampYP=sampYP, ssGrid=sGrid, timeGrid)
{
	# it calculates posterior marginal means  and percentiles of the average increments over a time period defined
	# by timeGrid for each point of the spatial grid
	indexS=sGrid*timeGrid
	nSim<-length(ssampYP)
	meanDsamp<-sapply(c(1:nSim), FUN=function(i,lsamp, index) 
		{return(tapply(lsamp[[i]][index>0],FUN=mean,INDEX=index[index>0]))},
		 lsamp=ssampYP, index=indexS)
	#dim(meanDsamp) is nP=844 x nSim=500
	meanDmean<-apply(meanDsamp, MARGIN=1, FUN=mean)
	meanDp95<-apply(meanDsamp, MARGIN=1, FUN=quantile, p=0.95)
	meanDp05<-apply(meanDsamp, MARGIN=1, FUN=quantile, p=0.05)
	return(distAv=cbind(meanDmean, meanDp95, meanDp05))
}



decGrid<-1+floor((tGrid-3)/10)
decGrid[decGrid==0]<-1
decGrid<-factor(decGrid, levels=c(1:6), labels=c("D1", "D2", "D3","D4","D5","D6"))
# decGrid is an indicator of the decades, defined as 10 years starting by the  last year
# but  "D1" contais  the first 12 years

tGridD6<-(decGrid=="D6")
MargDistAvD6<-Marg_Dist_Averages(timeGrid=tGridD6)
tGridD5<-(decGrid=="D5")
MargDistAvD5<-Marg_Dist_Averages(timeGrid=tGridD5)
tGridD4<-(decGrid=="D4")
MargDistAvD4<-Marg_Dist_Averages(timeGrid=tGridD4)
tGridD3<-(decGrid=="D3")
MargDistAvD3<-Marg_Dist_Averages(timeGrid=tGridD3)
tGridD2<-(decGrid=="D2")
MargDistAvD2<-Marg_Dist_Averages(timeGrid=tGridD2)
tGridD1<-(decGrid=="D1")
MargDistAvD1<-Marg_Dist_Averages(timeGrid=tGridD1)

MargDistAvTot<-Marg_Dist_Averages(timeGrid=1)


require(tidyverse)
require("rnaturalearth")
require("sf")
require("sp")

# Illustrating how to map the summaries
mapSpatial(Z= MargDistAvD2[,3], 
		CoordS=Ucoords,
		ref =  1.5, 
          zlim = c(0,3), 
          picture.name="G:/Mi unidad/Espacial/E9RecordsExpl/Graphs/AIncP05D2.png",
          title = "D2   P05", #"D6   P95" "D6   P05"
          legend.name = "",
          save =TRUE,
          label = "",
          ggrid=gridd,
          bbackground=background)




### Posterior distribution of the difference of the average increments over two periods
### at a spatial grid


Marg_Dist_Dif_Averages<-function(ssampYP=sampYP, ssGrid=sGrid, timeGrid1,timeGrid2)
{
	# it calculates posterior marginal means  and percentiles of the difference of the average increments over a time period defined
	# by timeGrid for each point of the spatial grid
	
	nSim<-length(ssampYP)
	indexS1=sGrid*timeGrid1
	meanDsamp1<-sapply(c(1:nSim), FUN=function(i,lsamp, index) 
		{return(tapply(lsamp[[i]][index>0],FUN=mean,INDEX=index[index>0]))},
		 lsamp=ssampYP, index=indexS1)
	indexS2=sGrid*timeGrid2
	meanDsamp2<-sapply(c(1:nSim), FUN=function(i,lsamp, index) 
		{return(tapply(lsamp[[i]][index>0],FUN=mean,INDEX=index[index>0]))},
		 lsamp=ssampYP, index=indexS2)
	meandif<-apply(meanDsamp2-meanDsamp1, MARGIN=1, FUN=mean)
	p05dif<-apply(meanDsamp2-meanDsamp1, MARGIN=1, FUN=quantile, p=0.05)
	p95dif<-apply(meanDsamp2-meanDsamp1, MARGIN=1, FUN=quantile, p=0.95)
	return(distDifAv=cbind(meandif, p95dif, p05dif))
}

MargDistDifAvD6D2<-Marg_Dist_Dif_Averages(timeGrid1=tGridD2, timeGrid2=tGridD6)
MargDistDifAvD6D3<-Marg_Dist_Dif_Averages(timeGrid1=tGridD3, timeGrid2=tGridD6)
MargDistDifAvD6D4<-Marg_Dist_Dif_Averages(timeGrid1=tGridD4, timeGrid2=tGridD6)
MargDistDifAvD6D5<-Marg_Dist_Dif_Averages(timeGrid1=tGridD5, timeGrid2=tGridD6)

require(tidyverse)
require("rnaturalearth")
require("sf")
require("sp")

# Illustrating how to map the summaries
mapSpatial(Z= MargDistDifAvD6D5[,3], 
		CoordS=Ucoords,
		ref =  0, 
          zlim = c(-0.5,0.25), 
          picture.name="G:/Mi unidad/Espacial/E9RecordsExpl/Graphs/ADifIncP05D6D5.png",
          title = "P05  D6-D5", #"D6   P95" "D6   P05"
          legend.name = "",
          save =TRUE,
          label = "",
          ggrid=gridd,
          bbackground=background)




### Posterior probabilities of the average increments over one month being higher than in another
### at a spatial grid


Prob_Comp_Month_Averages<-function(ssampYP=sampYP, ssGrid=sGrid,timeGrid, monthGrid=MonthTSGrid,
	nnP=nP, nY=10)
{
	nSim<-length(ssampYP)
	meanMtssamp<-sapply(c(1:nSim), FUN=function(i,lsamp, index, indexSel) 
	{return(tapply(lsamp[[i]][indexSel],FUN=mean,INDEX=droplevels(index[indexSel])))},
	 lsamp=ssampYP, index=monthGrid, indexSel=timeGrid)

	meaJntssamp<-meanMtssamp[seq(1,nnP*nY*3,by=3),]
	meaJltssamp<-meanMtssamp[seq(2,nnP*nY*3,by=3),]
	meaAgtssamp<-meanMtssamp[seq(3,nnP*nY*3,by=3),]
	pCompAgJl<-sapply(1:nnP, FUN=function(i, mat, index){return(mean(mat[index==i,]))}, 
		mat=(meaAgtssamp>meaJltssamp), index=rep(1:nnP, each=nY))
	pCompAgJn<-sapply(1:nnP, FUN=function(i, mat, index){return(mean(mat[index==i,]))}, 
		mat=(meaAgtssamp>meaJntssamp), index=rep(1:nnP, each=nY))
	pCompJlJn<-sapply(1:nnP, FUN=function(i, mat, index){return(mean(mat[index==i,]))}, 
		mat=(meaJltssamp>meaJntssamp), index=rep(1:nnP, each=nY))
	return(probCompMonth=cbind(pCompAgJl,pCompAgJn,pCompJlJn))
}





monthGrid<-factor(1+(lGrid>30)+(lGrid>61),levels=c(1:3), labels=c("Jn", "Jl", "Ag"))
MonthTSGrid<-interaction(monthGrid,tGrid,sGrid)

probCompMonthD6<-Prob_Comp_Month_Averages(timeGrid=tGridD6, monthGrid=MonthTSGrid)
probCompMonthD5<-Prob_Comp_Month_Averages(timeGrid=tGridD5, monthGrid=MonthTSGrid)
probCompMonthD4<-Prob_Comp_Month_Averages(timeGrid=tGridD4, monthGrid=MonthTSGrid)
probCompMonthD3<-Prob_Comp_Month_Averages(timeGrid=tGridD3, monthGrid=MonthTSGrid)
probCompMonthD2<-Prob_Comp_Month_Averages(timeGrid=tGridD2, monthGrid=MonthTSGrid)
probCompMonthD1<-Prob_Comp_Month_Averages(timeGrid=tGridD1, monthGrid=MonthTSGrid)
probCompMonthTot<-Prob_Comp_Month_Averages(timeGrid=c(1:length(tGrid)), nY=nT, monthGrid=MonthTSGrid)

require(tidyverse)
require("rnaturalearth")
require("sf")
require("sp")

# Illustrating how to map the summaries
mapSpatial(Z=probCompMonthTot[,3], 
		CoordS=Ucoords,
		ref =  0.5,  
		zlim = c(0.38,0.6),  #c(0.45,0.55), #
		picture.name="G:/Mi unidad/Espacial/E9RecordsExpl/Graphs/ProbCompAllJlJn.png",
		title = "Prob. Inc. July-June  All period", 
		legend.name = "",
		save = TRUE,
		label = "",
		ggrid=gridd,
		bbackground=background)



### Averages over space


Marg_Samp_Averages<-function(ssampYP=sampYP, indGrid)
{
	# it calculates samples of the average increments levels defined in indGrid
	
	nSim<-length(ssampYP)
	meanYsamp<-sapply(c(1:nSim), FUN=function(i,lsamp, index) 
			{return(tapply(lsamp[[i]],FUN=mean,INDEX=index))},
	 		 lsamp=ssampYP, index=indGrid)  
	return(meanYsamp)
}

#Average for each year
meanTYsamp<-Marg_Samp_Averages(ssampYP=sampYP, indGrid=tGrid) #  dim(meanTYsamp) 61 500
#Average for each year and site
tsGrid<-interaction(tGrid,sGrid)
meanTSYsamp<-Marg_Samp_Averages(ssampYP=sampYP, indGrid=tsGrid) # dim(meanTSYsamp) 61*844=51484   500



#Summaries  and graphs in the text based on the samples of averages

Marg_Summaries_Averages<-function(mmeanYsamp, q1=0.05, q2=0.95)
{
	# it calculates mean and percentiles of the average increments levels defined in indGrid

	meanYmean<-apply(mmeanYsamp, MARGIN=1, FUN=mean)
	meanYP95<-apply(mmeanYsamp, MARGIN=1, FUN=quantile, p=q2)
	meanYP05<-apply(mmeanYsamp, MARGIN=1, FUN=quantile, p=q1)
	return(list(meanYmean=meanYmean, meanYP95=meanYP95,meanYP05=meanYP05))
}

summTAve<-Marg_Summaries_Averages(mmeanYsamp=meanTYsamp)

#Graph Fig 10a
years<-c(1961:2021)
png(file="G:/Mi unidad/Espacial/E9RecordsExpl/Graphs/Annualmean.png")
par(mar=c(5.1, 5, 4.1, 2.1))
plot(years, summTAve$meanYmean, type="l", xlab="Year", ylab=expression(bar(J)[t ~ "," ~ JJA](S)), lwd=3, 
	cex.lab=1.5, cex.axis=1.5)
polygon(c(years,rev(years)), c(summTAve$meanYP05,rev(summTAve$meanYP95)), col="gray",border=NA)
lines(years, summTAve$meanYmean,lw=2)
dev.off()



Marg_F_Averages<-function(mmeanYsamp, ref, ind)
{
	# it calculates P(J>ref) of the average increments levels defined in indGrid
	Pref<-apply(mmeanYsamp, MARGIN=2, FUN=function(x, iind, rref){
			return(tapply((x>rref), INDEX=iind, FUN=mean))}, iind=ind, rref=ref)
	return(Pref)
}



indD<-rep(c(rep(1,11),rep(2:6, each=10)), nP)
#inD defines 6 levels defined by the decades
PD20<-Marg_F_Averages(meanTSYsamp, ref=2.0, ind=indD) # dim(P20)  6 500
PD17<-Marg_F_Averages(meanTSYsamp, ref=1.75, ind=indD)
PD15<-Marg_F_Averages(meanTSYsamp, ref=1.5, ind=indD)
PD12<-Marg_F_Averages(meanTSYsamp, ref=1.25, ind=indD)
PD10<-Marg_F_Averages(meanTSYsamp, ref=1.0, ind=indD)
PD7<-Marg_F_Averages(meanTSYsamp, ref=0.75, ind=indD)
PD5<-Marg_F_Averages(meanTSYsamp, ref=0.5, ind=indD)

mPD<-matrix(NA, nrow=7,ncol=6)
#matrix containing the mean of the probabilities for the 7 considered range of values  by decade
mPD[1,]<-apply(PD5, MARGIN=1, FUN=mean)
mPD[2,]<-apply(PD7, MARGIN=1, FUN=mean)
mPD[3,]<-apply(PD10, MARGIN=1, FUN=mean)
mPD[4,]<-apply(PD12, MARGIN=1, FUN=mean)
mPD[5,]<-apply(PD15, MARGIN=1, FUN=mean)
mPD[6,]<-apply(PD17, MARGIN=1, FUN=mean)
mPD[7,]<-apply(PD20, MARGIN=1, FUN=mean)

#Graph Fig 10b
#png(file="G:/Mi unidad/Espacial/E9RecordsExpl/Graphs/mProbRef.png")
plot(seq(0.5, 2, by=0.25), mPD[,1], type="l", ylim=c(0,1),col="yellow", 
	ylab="P(J>x)", xlab="x",cex.lab=1.5, cex.axis=1.5)
lines(seq(0.5, 2, by=0.25), mPD[,2], col="orange")
lines(seq(0.5, 2, by=0.25), mPD[,3], col="green")
lines(seq(0.5, 2, by=0.25), mPD[,4], col="red")
lines(seq(0.5, 2, by=0.25), mPD[,5], col="blue")
lines(seq(0.5, 2, by=0.25), mPD[,6], col="purple")
legend("topright",legend=c("D1", "D2", "D3", "D4", "D5", "D6"),lwd=2,
	lty=1, col=c("yellow", "orange", "green","red", "blue", "purple"),bty="n")
#dev.off()


#Graph Fig 11a and 11b
#png(file="G:/Mi unidad/Espacial/E9RecordsExpl/Graphs/DenProbRef10.png")
plot(density(PD10[1,]),type="l", xlim=c(0,1), xlab="p", ylim=c(0,35),
	main="Av. Inc. >1",lwd=2, col="yellow")
lines(density(PD10[2,]),col="orange", lwd=2)
lines(density(PD10[3,]),col="green", lwd=2)
lines(density(PD10[4,]),col="red", lwd=2)
lines(density(PD10[5,]),col="blue", lwd=2)
lines(density(PD10[6,]),col="purple", lwd=2)
legend("topleft",legend=c("D1", "D2", "D3", "D4", "D5", "D6"),lwd=2,
	lty=1, col=c("yellow", "orange", "green","red", "blue", "purple"),bty="n")
#dev.off()
#png(file="G:/Mi unidad/Espacial/E9RecordsExpl/Graphs/DenProbRef15.png")
plot(density(PD15[1,]),type="l", xlim=c(0,1), xlab="p", ylim=c(0,35),
	main="Av. Inc. >1.5",lwd=2, col="yellow")
lines(density(PD15[2,]),col="orange", lwd=2)
lines(density(PD15[3,]),col="green", lwd=2)
lines(density(PD15[4,]),col="red", lwd=2)
lines(density(PD15[5,]),col="blue", lwd=2)
lines(density(PD15[6,]),col="purple", lwd=2)
legend("topright",legend=c("D1", "D2", "D3", "D4", "D5", "D6"),lwd=2,
	lty=1, col=c("yellow", "orange", "green","red", "blue", "purple"),bty="n")
#dev.off()

# RData file containing all the summaries compued in this file
save(MargDistAvD1, MargDistAvD2,MargDistAvD3,MargDistAvD4,MargDistAvD5,MargDistAvD6,MargDistAvTot,
	MargDistDifAvD6D2, MargDistDifAvD6D3, MargDistDifAvD6D4, MargDistDifAvD6D5, 
	probCompMonthD1,probCompMonthD2,probCompMonthD3,probCompMonthD4,probCompMonthD5,probCompMonthD6,probCompMonthTot,
	meanTYsamp, meanTSYsamp, summTAve,mPD,
	file="G:\\Mi unidad\\Espacial\\E9RecordsExpl\\Datos_Z\\AInc_Marg_Summaries.RData")



#Means over climate zones 

Marg_mean_Averages_list<-function(ssamYP, indGrid=tzGrid,nnT=nT)
	# it computes  the average over the  simulations of each element in the vector (sample of the responses)
	mean_lts<-Reduce(`+`,ssampYP)/length(ssampYP)
	mean_tz<-tapply(mean_lts, FUN=mean, INDEX=indGrid)  
	Mmean_tz<-do.call(cbind,split(mean_tz, rep(1:nnT, 4)))
	return(Mmean_tz)
}

tzGrid<-interaction(tGrid,zGrid)
MmeanZones<-Marg_mean_Averages_list(ssamYP=samYP, indGrid=tzGrid)


#Figure  9.b Supplementary material
#png(file="G:/Mi unidad/Espacial/E9RecordsExpl/Graphs/AnnualmeanZone.png")
par(mar=c(5.1, 5, 4.1, 2.1))
plot(year, meanYmean, type="n", xlab="Year", ylab=expression(bar(J)[JJA ~ "," ~ t](Z)),
	ylim=range(MmeanZones), lwd=1)
sapply(c(1:4), FUN=function(i, mat, colZ=c("black", "blue", "red", "green"))
       {lines(year, mat[i,],col=colZ[i],lwd=2)}, mat=MmeanZones)
legend("topright", legend=c("Inland", "North C.", "East C", "South C."), col=c(colZ), lty=1, 
	lwd=2, bty="n")
#dev.off()

















