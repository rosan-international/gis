# GPS Displacement Procedure
# This function is an adaptation to R from the original Python code from Brendan Collis (Blue Raster) and released by DHS: https://dhsprogram.com/pubs/pdf/SAR7/SAR7.pdf. 

displace <- function(gps.vars, admin, samp_num, other_num){
  start <- Sys.time()
  map1 <- ggplot() + geom_polygon(data=countrymap, aes(x=long, y=lat, group=group), fill="darkgrey", alpha=0.3) +
    geom_point(data=mydata, aes_string(x=gps.vars[1], y=gps.vars[2])) + 
    ggtitle(paste0("Original GPS points: ", colnames=(gps.vars[1]), ", ", colnames=(gps.vars[2])))
  plot(map1)
  print("Summary Long/Lat statistics before displacement")
  print(summary(mydata[gps.vars]))
  coords <- as.data.frame(mydata[complete.cases(mydata[gps.vars]),gps.vars])
  URBAN_RURA <- as.data.frame(rep("R", length(coords[,1])))
  colnames(URBAN_RURA) <- "URBAN_RURA"
  crs_proj    <- admin@proj4string
  pts <- SpatialPointsDataFrame(coords, data=cbind(URBAN_RURA), proj4string = crs_proj)
  
  n <- length(pts)
  offset.dist <- ifelse(pts$URBAN_RURA == "U", 0.025, 0.05) # Tweak for offset distance
  rural <- which(pts$URBAN_RURA == "R")
  rur.n <- floor(0.01*length(rural))
  offset.dist[sample(rural, rur.n, replace = FALSE)] <- 0.1 # Tweak for offset distance
  r.pts0 <- list(0)
  for(i in 1:nrow(pts)){
    r.pts0[[i]]<-matrix(0,nrow=samp_num,ncol=2)
    #-- Buffer around point --#
    pdsc <- disc(radius = offset.dist[i], centre = c(coordinates(pts)[i,1],
                                                     coordinates(pts)[i,2]))
    pdsc <- as(pdsc, "SpatialPolygons")
    proj4string(pdsc) <- crs_proj
    #-- Intersection with admin --#
    int <- gIntersection(pdsc, admin)
    #-- Generating random point
    if(!is.null(int)){
      rpt <- csr(int@polygons[[1]]@Polygons[[1]]@coords, other_num)
      probs<-1/rdist(coordinates(pts[i,]),rpt)
      rpt<-rpt[sample(c(1:other_num),size=samp_num,prob=(probs/sum(probs))),]
      r.pts0[[i]] <- rpt
    }
    if(is.null(int)){
      rpt <- csr(pdsc@polygons[[1]]@Polygons[[1]]@coords, other_num)
      probs<-1/rdist(coordinates(pts[i,]),rpt)
      rpt<-rpt[sample(c(1:other_num),size=samp_num,prob=(probs/sum(probs))),]
      r.pts0[[i]] <- rpt
    }
  }
  #Arranging the Output
  if(samp_num==1){
    r.pts<-list(0)
    r.pts[[1]]<-matrix(0,nrow=n,ncol=2)
    for(k in 1:n){
      r.pts[[1]][k,]<-c(r.pts0[[k]])
    }
    r.pts[[1]]<- SpatialPoints(r.pts[[1]], CRS(as.character(crs_proj)))
  }
  if(samp_num>1){
    r.pts<-list(0)
    for(j in 1:samp_num){
      r.pts[[j]]<-matrix(0,nrow=n,ncol=2)
      for(k in 1:n){
        r.pts[[j]][k,]<-r.pts0[[k]][j,]
      }
      r.pts[[j]]<- SpatialPoints(r.pts[[j]], CRS(as.character(crs_proj)))
    }
  }
  
  mydata_displaced <- displace.merger(r.pts, gps.vars)
  plot(map1 + geom_point(alpha = 0.4, data=mydata_displaced, 
                         aes_string(x=gps.vars[1], 
                                    y=gps.vars[2]), 
                         colour = "red")
       + ggtitle(paste0("Displaced GPS points: ", colnames=(gps.vars[1]), ", ", colnames=(gps.vars[2]))))
  
  # Plot detail by focusing map on 1st and 3rd lon/lat quartiles 
  
  summary_long <- summary(mydata[[gps.vars[1]]])
  summary_lat <- summary(mydata[[gps.vars[2]]])
  plot(map1 + geom_point(alpha = 0.4, data=mydata_displaced, 
                         aes_string(x=gps.vars[1], 
                                    y=gps.vars[2]), 
                         colour = "red")       + 
         coord_map(xlim = c(summary_long[[2]], summary_long[[5]]),
                   ylim = c(summary_lat[[2]], summary_lat[[5]])) +
         ggtitle(paste0("Displaced GPS points (Interquartile Detail): ", 
                        colnames=(gps.vars[1]), ", ", colnames=(gps.vars[2]))))
  
  print("Summary Long/Lat statistics after displacement")
  print(summary(mydata_displaced[gps.vars]))
  end <- Sys.time()
  print(paste0("Processing time = ", round(end-start, 0)))
  return(mydata_displaced)
}
