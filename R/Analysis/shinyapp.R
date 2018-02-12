library(raster)
library(leaflet)
library(magrittr)
library(shiny)

# Link map data to files and prepare for plotting
dem<-raster(paste(dir_dem,"dem.tif",sep=""))
slope<-crop(raster(paste(dir_terrain,"slope.tif",sep="")),dem)
aspect<-crop(raster(paste(dir_terrain,"aspect.tif",sep="")),dem)

e<-extent(120000,240000,10000,140000)
dem<-crop(dem,e)
slope<-crop(slope,e)
aspect<-crop(aspect,e)

epsg<-"+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext +no_defs"
#epsg<-"+init=epsg:3857" # mercator
dem.epsg<-projectRaster(dem,crs=epsg,res=c(100,100),filename=paste(dir_dem,"demepsg.tif"),overwrite=TRUE)
slope.epsg<-projectRaster(slope,crs=epsg,res=c(100,100),filename=paste(dir_dem,"slopeepsg.tif"),overwrite=TRUE)
aspect.epsg<-projectRaster(aspect,crs=epsg,res=c(100,100),filename=paste(dir_dem,"aspectepsg.tif"),overwrite=TRUE)

map.s<-stack(dem.epsg,slope.epsg,aspect.epsg)


# Prepare user interface 
ui<-fluidPage(
  titlePanel("Demo Maps") ,
  selectInput(inputId = "county",label="Select Geographical area", 
              choices=  c("Cornwall" = "cornwall","Devon" = "devon","Somerset" = "other") ,
              selectize=TRUE  ) ,
  selectInput("variable","Select Map", c("Elevation" = "elev","Slope" = "slope", "Aspect" = "aspect") ,
              selectize=TRUE  ) ,
  leafletOutput("plotmap") ,
  actionButton("Re-display map","remap")
)

# Server functions using ui inputs 
server<-function(input,output) {
  mapdata<-reactive({
    
  }) ,
  output$map<-renderLeaflet({
    leaflet() %>% addTiles() %>%
      addRasterImage(dem.epsg,  opacity = 0.7)
  } ,
  #  output$table<-renderTable({summary(mapdata())})
}
  
 
shinyApp(ui=ui, server=server) 

 