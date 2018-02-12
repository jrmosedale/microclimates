vyds<-fortify(vineyards.epsg)
sel<-which(vyds$id=="1")
cv<-vyds[sel,]

leaflet() %>% setView(lng = -4.5, lat = 50.75, zoom = 8) %>% 
  addTiles() %>%
  #addProviderTiles("Esri.WorldImagery")  %>%  
  addRasterImage(cornwall.dem,opacity=0.4) %>%
  addPolygons(data=vineyards.ll,opacity=0.75,weight=1,fill=TRUE,fillColor="black",fillOpacity=0.75)

#Esri.WorldImagery Thunderforest.Lanscape

leaflet(vineyards.epsg)   %>% 
  addTiles()  %>%
  addPolygons()

m<-leaflet()  %>% setView(lng = -4.75, lat = 50.4, zoom = 10) 
m %>%  addTiles()  
m %>% addRasterImage(cornwall.dem,opacity=0.4) 
m %>% addPolygons(data=vineyards.epsg)
m

plot(dem.cornwall)
plot(vineyards.epsg,add=TRUE)

leaflet() %>% setView(lng = -4.5, lat = 50.75, zoom = 8) %>% 
  addTiles()  %>%  
  addPolygons(lng=c(-4.4,-4.45,-4.45),lat=c(50.8,50.6,50.8))
  addMarkers(lng=-4.4,lat=50.8)
  addPolygons(data=vineyards.epsg,stroke=TRUE,opacity=1,weight=5,fill=TRUE,color="black",fillColor="black") %>%
  #addRasterImage(cornwall.dem,opacity=0.4) %>%
    
  
  # PRINT COORDINATES OF A POLYGON
    x<-1 # polygon number
    vineyards.epsg@polygons[[x]]@Polygons[[1]]@coords
  
    coordinates(vineyards.epsg) # central coordinates for each polygon
    fortify(vineyards.epsg)  # create normal df
    
   # http://stackoverflow.com/questions/29803253/r-extracting-coordinates-from-spatialpolygonsdataframe  
    extractCoords <- function(sp.df)
    {
      results <- list()
      for(i in 1:length(sp.df@polygons[[1]]@Polygons))
      {
        results[[i]] <- sp.df@polygons[[1]]@Polygons[[i]]@coords
      }
      results <- Reduce(rbind, results)
      results
    }
