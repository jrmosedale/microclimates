library(raster)
library(leaflet)
library(magrittr)
library(shiny)

# Loads raster stack and displays - allows modification of visibility
dir_shinydata<-paste(root,"shinydata/",sep="")
pal <-colorQuantile(c("#B2FF66","#66CC00","#4C9900","#336600","#193300"),
                    NULL, n = 5, na.color="#FFFFFF")
res<-c(100,100)

lonlat_to_rasterxy<-function(x,y,zone=0){ #http://stackoverflow.com/questions/18639967/converting-latitude-and-longitude-points-to-utm
  xy <- data.frame(ID = 1:length(x), X = x, Y = y)
  coordinates(xy) <- c("X", "Y")
  proj4string(xy) <- CRS("+init=epsg:4326")  ## lat lon
  res <- spTransform(xy, CRS(epsg)) # convert to raster xy
  return(as.data.frame(res))
}

# Function that will return file nam
getmap<-function(var){ 
  if (var=="elevation") raster(paste(dir_shinydata,"demepsg.tif",sep="")) else 
    if (var=="slope") raster(paste(dir_shinydata,"slopeepsg.tif",sep="")) else
      if (var=="aspect") raster(paste(dir_shinydata,"aspectepsg.tif",sep=""))
}

getcoord<-function(county){ # returns list of xmin,xmax,ymin,ymax coordinates
  if (var=="cornwall") c(1,2,3,4) else 
    if (var=="devon") c(1,2,3,4) else
      if (var=="other") c(1,2,3,4) 
}

getfile<-function(var){ 
  if (var=="elevation") paste(dir_shinydata,"demepsg.tif",sep="") else 
    if (var=="slope") paste(dir_shinydata,"slopeepsg.tif",sep="") else
      if (var=="aspect") paste(dir_shinydata,"aspectepsg.tif",sep="")
}



# Prepare user interface 
ui<-fluidPage(
  fluidRow(titlePanel("Demo Maps")) ,
  fluidRow(
    column(3,
           selectInput(inputId = "county",label="Select Geographical area", 
                       choices=  c("Cornwall" = "cornwall","Devon" = "devon","Somerset" = "other") ,
                       selected="cornwall",selectize=TRUE  ) ,
           selectInput("variable","Select data", c("Elevation" = "elevation","Slope" = "slope", "Aspect" = "aspect") ,
                       selected="elevation",selectize=TRUE  ) ,
           checkboxInput("legend", "Show legend", FALSE)
    ),
    column(9,
           leafletOutput("map")  
    ) ),
  fluidRow(
    column(3,offset=3,
           plotOutput("hist1")
    ),
    column(3,
           plotOutput("hist2")
    ),
    column(3,
           plotOutput("hist3")
    )  )
  #actionButton("Re-display map","remap")
)


# Server functions using ui inputs 
server<-function(input,output) {
  # display base tiles - CHECK starting raster image displayed is same as selected in UI
  output$map<-renderLeaflet({
    leaflet() %>% addTiles()  %>% 
      addRasterImage(map.s,  color=pal, opacity = 0.7,project=FALSE) %>% 
      addRasterImage(getmap("elevation"),  color=pal, opacity = 0.7,project=FALSE) %>% 
      addRasterImage(getmap("elevation"),  color=pal, opacity = 0.7,project=FALSE) 
      
  })
  
  # set layer to be mapped to that chosen in UI
  chosenlayer <- reactive({getmap(input$variable)}) 
  
  # Observe to change data if UI changes
  observe({
    leafletProxy("map") %>% 
      clearImages() %>% 
      addRasterImage(chosenlayer(),color=pal, opacity = 0.7,project=FALSE)  %>% addScaleBar()
  })
  
  # Observe to change legend if requested
  observe({
    proxy <- leafletProxy("map")
    # Remove any existing legend, and only if the legend is
    # enabled, create a new one.
    proxy %>% clearControls()
    if (input$legend){
      pal<-pal
      proxy %>% addLegend(position = "bottomright",pal=pal,values=values(chosenlayer()) )
    }
  })  
  
  # Calculate value for mapclick - observeEvent(input$mouseover) observeEvent(input$MAPID_click)
  # Display lat lon of mouse click
  observe({  #Observer to show Popups on click
    click <- input$map_click
    if (!is.null(click)) {
      showpos(lon=click$lng, lat=click$lat,chosenlayer())
      
    }
  })
  
  showpos <- function(lon=NULL, lat=NULL,r) {#Show popup on clicks
    #Get value of the given cell
    xy<-lonlat_to_rasterxy(lon,lat)
    #print(xy$X)
    #print(xy$Y)
    cell <- cellFromXY(r, c(xy$X,xy$Y))
    #print(cell)
    value = as.numeric(r[cell])
    #print(value)
    content <- paste("Lon=", round(lon, 2),"; Lat=", round(lat, 2),"; ", input$variable,"=",round(value,2),sep="")
    proxy <- leafletProxy("map")
    #add Popup
    proxy %>% clearPopups() %>% addPopups(lon, lat, popup = content)
    # PLot hist
    # Calculate dem, slope and aspect values
    dem.value<-dem.epsg[cell]
    slope.value<-slope.epsg[cell]
    aspect.value<-aspect.epsg[cell]
    
    output$hist1<-renderPlot(dem.hist + geom_vline(xintercept = dem.value,color = "red", size=1))
    output$hist2<-renderPlot(slope.hist + geom_vline(xintercept = slope.value,color = "red", size=1))
    output$hist3<-renderPlot(aspect.hist + geom_vline(xintercept = aspect.value,color = "red", size=1))
  }
  
  output$hist1<-renderPlot(dem.hist)
  
}


shinyApp(ui=ui, server=server) 