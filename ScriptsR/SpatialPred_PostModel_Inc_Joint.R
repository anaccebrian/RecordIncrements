#####################################################
###Generating functions
#####################################################

Sim_Joint<-function(ssampYP=sampYP,IIsim=Isim)
{
	#Generate samples from the joint distribution using samples from the value of the increments (from the marginal
	# distribution) and samples of the occurrence of records
	#sampYP[[1]] is a list with 500 elements ( corresponding to 1 simulation), each one vector of length  844*92*61
	#dim(Isim31) is an array  with dim: 500  92  61 844 formed by binary variables
	nSim<-length(ssampYP)
	sampYPJoint<-lapply(c(1:nSim), FUN=function(i,lsamp, indRecord) 
	   {return(lsamp[[i]]*as.vector(indRecord[i,,,])) },	
	lsamp=ssampYP, indRecord=IIsim)
	#sampYPJoint is a list  with nSim elements each one correpsonding to a simulation
	#length(sampYPJoint[[1]]) 4736528= nP*nT*nL
	return(sampYPJoint)
}


#Code to apply the function

load("C:\\Users\\PC\\Desktop\\AnaE9RecordsExpl\\DatosSim\\SimDataMarg.RData")
load("G:\\Mi unidad\\Espacial\\E9RecordsExpl\\Datos_Z\\Grid_Spat_Isim.RData")
sampYPJ<-Sim_Joint(ssampYP=sampYP,IIsim=Isim)
save(sampYPJ, file="C:/Users/PC/Desktop/AnaE9RecordsExpl/DatosSim/SimDataJoint.RData")



#######################################################
### Post model functions
#######################################################


Mean_Averages<-function(ssampYP=sampYPJ, iindexS=indexS,nnP=nP, nnL=nL)
{
	#compute the samples and the corresponding posterior mean of the average (over days within the year) of the 	# cumulative increments during the last 30 years at each point of the spatial grid
	nSim<-length(ssampYP)
	sumsampJaux<-sapply(c(1:nSim), FUN=function(i,lsamp, index) 
	   {return(tapply(lsamp[[i]],FUN=sum,INDEX=index))},
	   lsamp=ssampYP, index=iindexS)	
	#sumsampJaux contains nSim simulations of  the sum of the incrments in a given period at each dal l and site s
	sumsampJ<-sumsampJaux[2:(nnP*nnL+1),] 
	#removes the first position corresponding to the sum of observations with level0 in indexS 0
	AveJsamp<-apply(sumsampJ, MARGIN=2, FUN=function(x, ind)
		{return(tapply(x,INDEX=ind, FUN=mean))},
		ind=rep(1:nP, each=nL), simplify=TRUE)
	#For each simulation, computes the average,  over the days of the year, of the sum of the increments is calculated for each site s
	#The result is a  matrix nP x nSim
	AveJmean<-apply(AveJsamp, MARGIN=1, FUN=mean)
	#Computes the mean of the averages ( the mean over the simulations)

	return(list(AveJsamp=AveJsamp, AveJmean=AveJmean))
}


#Code to apply the function

load("C:\\Users\\PC\\Desktop\\AnaE9RecordsExpl\\DatosSim\\SimDataJoint.RData")
load("G:\\Mi unidad\\Espacial\\E9RecordsExpl\\Datos_Z\\Grid_Spat_Cov.RData")
source("G:\\Mi unidad\\Espacial\\E9RecordsExpl\\Ficheros_Z\\Function_MapSpatial.R")

lsGrid<-interaction(lGrid,sGrid)
indexS=factor(lsGrid, levels=c(0, levels(lsGrid)))
indexS[(tGrid<=31)]<-0 #keep last 30 years
#indexS  takes different values for each day and site s (but equal for all the observations in the year)

Avesum30J<-Mean_Averages(ssampYP=sampYPJ,  iindexS=indexS)

library(tidyverse)
library("sp")

mapSpatial(Z= Avesum30J$AveJmean,
		CoordS=Ucoords,
            ref =1.5, 
           	zlim =  c(0,3), # c(0.8, 4.1),
           	picture.name="G:/Mi unidad/Espacial/E9RecordsExpl/Graphs/AveJmean30.png",
           	title = "Period: 1992-2021 ", 
           	legend.name = "",
           	save = TRUE,
           	label = "",
           	ggrid=gridd,
           	bbackground=background)


# RData file containing all the summaries compued in this file
save(,file="G:\\Mi unidad\\Espacial\\E9RecordsExpl\\Datos_Z\\AInc_Joint_Summaries.RData")


