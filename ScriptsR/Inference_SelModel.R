

#####################################################
###Map of the spatial random effects
######################################################

library("sf")
library("sp")
library("tidyverse")
library("rnaturalearth")
library(INLA)


mapSpatialE <- function(coords_limit = data.frame(X = c(-10, 4), Y = c(35.5, 44)), 
                     ref = .5,  
                     zlim = NULL, 
                     picture.name = "photo.png",
                     title = expression(p["tl"](s)), 
                     legend.name = "",
                     save = FALSE,
                     breaks = 8:12 / 10,
                     label = "",
                     dist,
			   csize=0.1,
			   mmesh=mesh,
			   Msel) {
    
  coords_limit <- as(
      SpatialPointsDataFrame(coords = coords_limit, 
                             data = coords_limit,
                             proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")),
      'sf')
  
  #define grid of Spain
  spain <- ne_countries(scale = "large", country = "Spain", returnclass = "sf")
  
  spain_coords <- Polygons(
    list(Polygon(st_coordinates(spain)[st_coordinates(spain)[,"L2"] == 3,1:2])),
    ID = "spain")
  spain_coords <- SpatialPolygons(list(spain_coords))
  spain_coords <- as(spain_coords, "sf")
  st_crs(spain_coords) <- st_crs(spain)
  
  gridaux <- st_make_grid(spain, cellsize = csize, what = "centers")
  grid <- st_intersection(gridaux, spain_coords)
  mgrid<-t(sapply(grid, FUN=as.vector, simplify = TRUE))

  gprojSp<-inla.mesh.projector(mmesh,loc=mgrid)
  Z<-inla.mesh.project(gprojSp, Msel$summary.random$spatialInd$mean)

  # background
  world_map <- ne_countries(scale = "large", returnclass = 'sf')
  european_union <- c("Andorra","France","Morocco","Portugal","Spain","Gibraltar","Algeria")
  background <- 
    world_map %>% 
    filter(name %in% european_union) 
  
  map <- ggplot(data = background) + 
    geom_sf(fill="antiquewhite") + 
    xlab("Longitude") + ylab("Latitude") + ggtitle(title) +
    theme(panel.background = element_rect(fill = "aliceblue"),
          axis.text.x=element_text(size=6),
          axis.text.y=element_text(size=6,angle=90),
          axis.title=element_text(size=10,face="bold")) + 
    geom_tile(data = grid, ggplot2::aes(x = st_coordinates(grid)[,1], y = st_coordinates(grid)[,2], fill = Z)) + 
    scale_fill_gradient2(midpoint = ref, low = scales::muted("blue"), mid = "white", high = scales::muted("red"), space = "Lab", limits = zlim, name = legend.name)
  
  
  map <- map + coord_sf(xlim = st_coordinates(coords_limit)[,1], ylim = st_coordinates(coords_limit)[,2])# +
  #ggplot2::geom_point(ggplot2::aes(x = X, y = Y), data = data.frame(coords * 1000), color = "black")
  
  if (save) {
    ggplot2::ggsave(picture.name, map, width = 8.27 / 2, height = 11.69 / 4)
  }
  
  map
}

### Applciation of the mapSpatialE

mapSpatialE(coords_limit = data.frame(X = c(-10, 4), Y = c(35.5, 44)), 
         ref = 0,#sum(1 / 1:62),  
         zlim = c(-0.4, 0.55), 
         save = FALSE,
#        picture.name = "G:/Mi unidad/Espacial/E9RecordsExpl/Graphs/SpatialEffect.png",
         title = "Spatial random effect", 
         legend.name = "",
         breaks = 1 + 1:4 / 20,
         label = paste0(10 * 1:4 / 2, "%"),
         csize=0.1,
	   mmesh=mesh,
	   Msel=MSelT)



###############################
### Plot of the coefficients of the model
###############################


plot_coef<- function (M=MSelT,params=NULL, 
			    coef.names=NULL,
                      ci_level = 0.95, 
                      intercept = 0,
                      colors = "black", 
                      point.size = 5, 
                      legend.title = "Model") {
   
if (is.null(coef.names))  coef.names<-rownames(M$summary.fixed)[-1]

  tidies <- data.frame(estimate =M$summary.fixed[-1,1])
  tidies$conf.low  <- M$summary.fixed[-1,3]
  tidies$conf.high <- M$summary.fixed[-1,5]
  tidies$model     <- as.factor("M1")
  tidies$term      <- factor(1:length(coef.names))# factor(coef.names, levels = rev(coef.names))
print(tidies)
  
  dh <- as.numeric(!TRUE) * 0.5
  
  p <- ggplot(data = tidies, 
              aes(y = term, x = estimate, xmin = conf.low, xmax = conf.high, 
				colour = model)) 
  
  p <- p + ggplot2::geom_pointrange(
    aes(y = term, x = estimate, xmin = conf.low, xmax = conf.high, colour = model, shape = model), 
    position = ggplot2::position_dodge(width = dh), 
    fill = "white", fatten = 1)
  
  p <- p + geom_vline(xintercept = intercept, linetype = 2) + 
    scale_colour_manual(values = colors, limits = rev(levels(tidies$model)), 
                        breaks = rev(levels(tidies$model)), 
                        labels = rev(levels(tidies$model))) + 
    drop_y_gridlines()
  
  p <- p + 
    scale_y_discrete(limits = levels(tidies$term), name = "Parameters", labels=coef.names) +
    scale_x_continuous(breaks=seq(-2,2,.2)) 
  
  yrange <- ggplot_build(p)$layout$panel_params[[1]]$y.range
  xrange <- ggplot_build(p)$layout$panel_params[[1]]$x.range
  
  if (yrange[2] <= (length(unique(tidies$term)) + 0.8)) {
    upper_y <- length(unique(tidies$term)) + 0.8
  } else {
    upper_y <- yrange[2]
  }
  
  lower_y <- 0.8
  p <- p + coord_cartesian(ylim = c(lower_y, upper_y), xlim = xrange, expand = FALSE) +
    theme_minimal() + theme(legend.position = "none",axis.title.y = element_blank(), axis.title.x = element_text(size = point.size),
		axis.text.y = element_text(size = point.size), axis.text.x = element_text(size = point.size)) + 
    xlab(ifelse(FALSE, no = "Estimate", yes = "exp(Estimate)"))  
  
  return(p)
}


### Aplication of the function plot_coef to the selected model

library("ggplot2")
library("jtools")
library(latex2exp)

coefnames<-c("log(t)",TeX("$(log(t))^2$"),TeX("$I_{t,l-1}(s)$"),"logElev","logDC","log(t)*logDC",
		TeX("$I_{t,l-1}(s)*logDC$"), TeX("$I^{Sp}_{t,l-1}$"))
aux.ggplot<-plot_coef(M=MSelT,params=NULL, 
		coef.names=coefnames,
            ci_level = 0.95, 
            intercept = 0,
            colors = "red", 
            point.size = 20, 
            legend.title = "Model") 
aux.ggplot
#gplot2::ggsave("G:/Mi unidad/Espacial/E9RecordsExpl/Graphs/coefmodel.pdf", aux.ggplot, 
#                width = 1.2 * 8.27, height = 1.2 * 11.69 / 2)


