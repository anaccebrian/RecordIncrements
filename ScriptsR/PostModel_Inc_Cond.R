
### Functions

CondMeanInc<-function(ssampYP,period, nnSim=nSim,ttGrid=tGrid, ssGrid=sGrid)
{
	nSim<-length(ssampYP)
	indexS=ssGrid*((ttGrid>period[1])&(ttGrid<=period[2]))
	meanDSsamp<-sapply(c(1:nnSim), FUN=function(i,lsamp, index) 
	{return(tapply(lsamp[[i]][index>0],FUN=mean,INDEX=index[index>0]))},
	      lsamp=ssampYP, index=indexS)
	meanDmean<-apply(meanDSsamp, MARGIN=1, FUN=mean)
	return(meanDmean)
}


###Running the functions

require(tidyverse)
require("rnaturalearth")
require("sf")
require("sp")


load("G:\\Mi unidad\\Espacial\\E9RecordsExpl\\Datos_Z\\Grid_Spat_Cov.RData")
load("C:\\Users\\PC\\Desktop\\AnaE9RecordsExpl\\DatosSim\\SimDataCond.RData")
source("G:\\Mi unidad\\Espacial\\E9RecordsExpl\\Ficheros_Z\\Function_MapSpatial.R")

meanDmean0D6<-CondMeanInc(ssampYP=sampYP0,period=c(52,62))
meanDmean1D6<-CondMeanInc(ssampYP=sampYP1,period=c(52,62))
mapSpatial(Z=((meanDmean1D6-meanDmean0D6)), 
            ref = 0.5,   
            zlim=c(0.23, 0.73),    
            picture.name = "G:/Mi unidad/Espacial/E9RecordsExpl/Graphs/DifCondMeanD6.png",
            title = "Decade 2012-21", 
            legend.name = "",
            save = FALSE,
            label = "",
            ggrid=gridd,
            bbackground=background)

meanDmean0D5<-CondMeanInc(ssampYP=sampYP0,period=c(42,52))
meanDmean1D5<-CondMeanInc(ssampYP=sampYP1,period=c(42,52))
mapSpatial(Z=((meanDmean1D5-meanDmean0D5))
            ref = 0.5,  
            zlim=c(0.23, 0.73),  
            picture.name = "G:/Mi unidad/Espacial/E9RecordsExpl/Graphs/DifCondMeanD5.png",
            title = "Decade 2002-11", 
            legend.name = "",
            save = TRUE,
            label = "",
            ggrid=gridd,
            bbackground=background)


meanDmean0D4<-CondMeanInc(ssampYP=sampYP0,period=c(32,42))
meanDmean1D4<-CondMeanInc(ssampYP=sampYP1,period=c(32,42))
mapSpatial(Z=((meanDmean1D4-meanDmean0D4)), 
            ref = 0.5,   
            zlim=c(0.23, 0.73),
            picture.name = "G:/Mi unidad/Espacial/E9RecordsExpl/Graphs/DifCondMeanD4.png",
            title = "Decade 1992-2001", 
            legend.name = "",
            save = TRUE,
            label = "",
            ggrid=gridd,
            bbackground=background)

meanDmean0D3<-CondMeanInc(ssampYP=sampYP0,period=c(22,32))
meanDmean1D3<-CondMeanInc(ssampYP=sampYP1,period=c(22,32))
mapSpatial(Z=((meanDmean1D3-meanDmean0D3)), 
            ref = 0.5,  
            zlim=c(0.23, 0.73),
            picture.name = "G:/Mi unidad/Espacial/E9RecordsExpl/Graphs/DifCondMeanD3.png",
            title = "Decade 1982-1991", 
            legend.name = "",
            save = TRUE,
            label = "",
            ggrid=gridd,
            bbackground=background)

mmeanDmean<-cbind(meanDmean0D3,meanDmean1D3,meanDmean0D4,meanDmean1D4,
                  meanDmean0D5,meanDmean1D5,meanDmean0D6,meanDmean1D6)
save(mmeanDmean, file="G:\\Mi unidad\\Espacial\\E9RecordsExpl\\Datos_Z\\Con_Mean_Inc.RData")

