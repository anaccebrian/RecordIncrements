

###########################
#Functions for maping features
#############################

mapSpatial <- function(coords_limit = data.frame(X = c(-10, 4), Y = c(35.5, 44)), 
                     Z,
                     Zobs,
			   CoordS=Ucoords,
                     ref = 40,  
                     zlim = NULL, 
                     picture.name = "photo.png",
                     title = "Mean increments", 
                     legend.name = "",
                     save = TRUE,
                     label = "",
                     ggrid=grid,
                     bbackground=background)
{
    
  coords_limit <- as(
      SpatialPointsDataFrame(coords = coords_limit, data = coords_limit,
                   proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")),'sf')
  
  map <- ggplot(data = bbackground) + 
    geom_sf(fill="antiquewhite") + 
    xlab("Longitude") + ylab("Latitude") + ggtitle(title) +
    theme(panel.background = element_rect(fill = "aliceblue"),
          axis.text.x=element_text(size=6),
          axis.text.y=element_text(size=6,angle=90),
          axis.title=element_text(size=10,face="bold")) + 
    geom_tile(data = ggrid, ggplot2::aes(x = st_coordinates(ggrid)[,1], y = st_coordinates(ggrid)[,2], fill = Z)) + 
    scale_fill_gradient2(midpoint = ref, low = scales::muted("blue"), mid = "white", high = scales::muted("red"), 
	space = "Lab", limits = zlim, name = legend.name)
  
  
  map <- map + coord_sf(xlim = st_coordinates(coords_limit)[,1], 
		ylim = st_coordinates(coords_limit)[,2])# +
  #ggplot2::geom_point(ggplot2::aes(x = X, y = Y), data = data.frame(coords * 1000), color = "black")
  


 if (!missing(Zobs)) {
  stations=data.frame(CoordS, Zobs) 
    map <- map +
      geom_point(data=stations,aes(x = CoordS[,1], 
		 y = CoordS[,2], fill = Zobs),
        color = "black", pch=21, size=3)
  } 

  if (save) {
    ggplot2::ggsave(picture.name, map, width = 8.27 / 2, height = 11.69 / 4)
  }
  map
}

