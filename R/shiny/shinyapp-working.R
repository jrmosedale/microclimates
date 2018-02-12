# Install dev version of leaflet for R if required
#if (!require('devtools')) install.packages('devtools')
#devtools::install_github('rstudio/leaflet')

library(shinydashboard)
library(raster)
library(leaflet)
library(magrittr)
library(shiny)
library(colorspace)
library(rgeos)
library(ggplot2)
#library(grDevices)

# To use with uploaded app set work directory to directory of app and modify data locations accordingly
# setwd( "/Users/yimizhao/Desktop/Study/Data_science/Shiny/census_app")

root<-"~/Documents/Exeter/Data2015/"; in.root<-"~/Documents/Exeter/Data2015/"

dir_shinydata<-paste(root,"shinydata/",sep="")
#dir_swfiles<-paste(root,"shinydata/swfiles/",sep="")
dir_cornwall<-paste(dir_shinydata,"cornwall/",sep="")
dir_devon<-paste(dir_shinydata,"devon/",sep="")
dir_dorset<-paste(dir_shinydata,"dorset/",sep="")
dir_somerset<-paste(dir_shinydata,"somerset/",sep="")
#dir_input<-paste(dir_shinydata,"shinydata/input/",sep="")  

latlong = "+init=epsg:4326"
#epsg<-"+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext +no_defs"

#vineyards.shp<- readOGR(dsn = "/Users/jonathanmosedale/Documents/Exeter/Data2015/vineyards", layer = "vineyardplots")
#vineyards.shp<-spTransform(vineyards.shp,latlong)

load("/Users/jonathanmosedale/Documents/Exeter/Data2015/shinydata/vineyards.RData") # loads vineyards.shp

res<-c(100,100)

t1<-theme(                              
  plot.background = element_blank(), 
  panel.grid.major = element_blank(), 
  panel.grid.minor = element_blank(), 
  panel.border = element_blank(), 
  panel.background = element_blank(),
  axis.line = element_line(size=1)
)
      

############################################
# Define UI
############################################

ui <- dashboardPage(
  dashboardHeader(title="Viticulture: climate risk maps", titleWidth=300),
  
  dashboardSidebar(
    sidebarMenuOutput("menu")
  ),
  
  dashboardBody(
    tabItems(
      tabItem(tabName="about",
              h2("Something about this app")
      ),
      
      tabItem(tabName="explore",
              fluidRow(height=100,
                box(width=3 ,
                    selectInput(inputId = "county",label="Geographical area", 
                                   choices=  c("Cornwall" = "cornwall","Devon" = "devon","Dorset" = "dorset", "Somerset" = "somerset") ,
                                   selected="cornwall",selectize=TRUE  ) ,
                    selectInput(inputId = "year",label="Time period", 
                                   choices=  c("2009" = "2009","2010" = "2010","2011" = "2011") ,
                                   selected="2009",selectize=TRUE  ) 
                ), # box
                tabBox(id="tabVariables", width=8 ,
                       
                	tabPanel("Growing season temperature",
                    selectInput("variable1","Select data", c("Choose one"="","Degree days (base 10C)" = "gdd10","Degree days (base 5C)" = "gdd5", "Mean temperature" = "meant", "Max temperature"="maxt","Min temperature"="mint","No. days > 20C"="num20","No. days > 25C"="num25","No. days > 30C"="num30") ,selected="",selectize=TRUE  ),
                    tags$p("Growing season: 1 April to 31 October."),
                    value="temp" 
                  ) , # tabPanel
                    
                  tabPanel("Frost Risks",  
                  	 selectInput("variable2","Select data", c("Choose one"="","Last spring frost" = "spfr","First autumn frost" = "autfr", "Frost free period" = "frfree") ,selected="",selectize=TRUE  ),
                  	  tags$p("Frost event defined as < 1C. Dates given as days of year (1 to 366)"),
                  	 value="frost"
             		  ), # tabPanel
             		
             		  tabPanel("Flowering Risks",  
                  	 selectInput("variable3","Select data", c("Choose one"="","Flowering conditions" = "flcond") ,selected="",selectize=TRUE  ),
                  	  tags$p("Bad flowering defined as ..."),
                  	 value="flower"
             		  ), # tabPanel

             		  tabPanel("Terrain",  
                     selectInput("variable4","Select data", c("Choose one"="","Elevation" = "elevation","Slope" = "slope", "Aspect" = "aspect"), selected="",selectize=TRUE  ), 
                     value="terrain"
                    )  # tabPanel
                	
                  )  #  # tabBox
              ) ,  #fluidRow
              
              fluidRow( height=200,
                  box( leafletOutput("map"), width=10
                  ) , # box
                  box(title="Map controls", width=2,
                      sliderInput("visibility","Set transparency:",min=0,max=1,value=0.7),
                      checkboxInput("legend", "Show legend", FALSE) ,
                      checkboxInput("vineyards", "Show vineyards", FALSE),
                      textOutput("maptext")
                  )  # box"
                ) , # fluidRow
                
              fluidRow(height=75,
                box(title="Selected cell data", width=3,
                    textOutput("celltext") ),
                box(title="Distribution of values", width=3,
                    plotOutput("hist",height="200px") ),
                box(title="Seasonal variation", width=3,
                    plotOutput("timeseriesplot",height="200px") )
              ) # fluidRow
      ), 
      
      tabItem(tabName="search",
              fluidRow(height=30, h4("  Define area & time period"), 
                            box(width=4,
                                selectInput(inputId = "searchcounty",label="Geographical area", 
                                        choices=  c("Cornwall" = "cornwall","Devon" = "devon","Dorset" = "dorset", "Somerset" = "somerset") ,
                                        selected="cornwall",selectize=TRUE  ) 
                                ),
                            box(width=4, 
                            selectInput(inputId = "searchyear",label="Time period", 
                                        choices=  c("2009" = "2009","2010" = "2010","2011" = "2011") ,
                                        selected="2009",selectize=TRUE ) 
                            )
                        ), # fluidRow
              fluidRow(height=60, h4("  Choose search variables & range:"),
                        box(width=4,
                            selectInput("search1",NULL, c("None chosen"="","Degree days (base 10C)" = "gdd10","Degree days (base 5C)" = "gdd5", 
                                                                   "Mean temperature" = "meant", "Max temperature"="maxt","Min temperature"="mint",
                                                                   "No. days > 20C"="num20","No. days > 25C"="num25","No. days > 30C"="num30",
                                                                   "Last spring frost" = "spfr","First autumn frost" = "autfr", "Frost free period" = "frfree",
                                                                   "Flowering conditions" = "flcond",
                                                                   "Elevation" = "elevation","Slope" = "slope", "Aspect" = "aspect"), 
                                                                    selected="",selectize=FALSE  ) ,
                            sliderInput("slider1", label=NULL, min = 0, 
                                        max = 100, value = c(20,50))
                        ) ,  # box
                       box( width=4,
                            selectInput("search2",NULL, c("None"="","Degree days (base 10C)" = "gdd10","Degree days (base 5C)" = "gdd5", 
                                                                   "Mean temperature" = "meant", "Max temperature"="maxt","Min temperature"="mint",
                                                                   "No. days > 20C"="num20","No. days > 25C"="num25","No. days > 30C"="num30",
                                                                   "Last spring frost" = "spfr","First autumn frost" = "autfr", "Frost free period" = "frfree",
                                                                   "Flowering conditions" = "flcond",
                                                                   "Elevation" = "elevation","Slope" = "slope", "Aspect" = "aspect"), 
                                        selected="",selectize=FALSE  ) ,
                           sliderInput("slider2", label=NULL, min = 0, 
                                       max = 100, value = c(20,50))
                       ) , # box
                       box( width=4,
                            selectInput("search3",NULL,c("None"="","Degree days (base 10C)" = "gdd10","Degree days (base 5C)" = "gdd5", 
                                                                   "Mean temperature" = "meant", "Max temperature"="maxt","Min temperature"="mint",
                                                                   "No. days > 20C"="num20","No. days > 25C"="num25","No. days > 30C"="num30",
                                                                   "Last spring frost" = "spfr","First autumn frost" = "autfr", "Frost free period" = "frfree",
                                                                   "Flowering conditions" = "flcond",
                                                                   "Elevation" = "elevation","Slope" = "slope", "Aspect" = "aspect"), 
                                        selected="",selectize=FALSE  ) ,
                            sliderInput("slider3", label=NULL, min = 0, 
                                        max = 100, value = c(20,50))
                      ) # box
              ) ,  #fluidRow
              fluidRow( height=200,
                        box( leafletOutput("searchmap"), width=10
                        ) , # box
                        box(title="Map controls", width=2,
                            actionButton(inputId="startsearch",icon("refresh"),label=" REFRESH MAP "),
                            p(),p(),
                            textOutput("searchmaptext"),
                            p(),p(),
                            sliderInput("visibility","Set transparency:",min=0,max=1,value=0.7),
                            # checkboxInput("legend", "Show legend", FALSE) ,
                            checkboxInput("vineyards2", "Show vineyards", FALSE) #,
                            #textOutput(paste(input$search1," between ", min(input$slider1)," and ",max(input$slider1),sep="")),
                            #textOutput(paste(input$search1," between ", min(input$slider2)," and ",max(input$slider2),sep="")),
                            #textOutput(paste(input$search1," between ", min(input$slider3)," and ",max(input$slider3),sep=""))
                            )  # box
              )  # fluidRow
      ) # tabItem
    ) # tabItems
  ) # dashboardBody
) # dashboardPage


####################################################################################################################################
# Server function
####################################################################################################################################
server<-function(input,output,session) {
  # side menu
  output$menu <- renderMenu({
    sidebarMenu(
      menuItem("About", tabName="about", icon=icon("info")),
      menuItem("Explore Maps", tabName="explore", icon=icon("map")),
      menuItem("Search Maps", tabName="search", icon=icon("search"))    )
  })
  
  #################################################
  # EXPLORE MAPS PAGE
  #################################################  
  # Outline of Map
  output$map<-renderLeaflet({
    leaflet() %>% setView(lng = -4.5, lat = 50.75, zoom = 8) %>% 
    addTiles()  %>%  
    clearImages() %>% 
    # addRasterImage(chosenlayer(),color=pal, opacity = 0.7,project=FALSE) %>%
    addScaleBar()
  })
  # Display what map shows
  output$maptext<-renderText({
    req(input$county)
    req(input$year)
    req(chosenvar())
    paste("Map showing ",chosenvar(),"for",input$county, "in",input$year,sep=" ")
  })
  
  # Observe to change map data if UI changes
  observe({
    leafletProxy("map") %>% 
      clearPopups() %>%
      clearImages() %>% 
      addRasterImage(chosenlayer(),color=mapcolour(chosenvar()), opacity = visibility(),project=FALSE)  
  })
  
  # Observe to change legend if requested
  observe({
    proxy <- leafletProxy("map")
    # Remove any existing legend, and only if the legend is
    # enabled, create a new one.
    proxy %>% clearControls()
    if (input$legend){
      #pal<-pal
      proxy %>% addLegend(position = "bottomright",pal=mapcolour(chosenvar()),values=values(chosenlayer()) )
    }
  })  
  
  # Observe to add Vineyards
  observe({
    proxy <- leafletProxy("map")
    # Remove any existing legend, and only if the legend is
    # enabled, create a new one.
     proxy %>% clearShapes()
    if (input$vineyards){
      #pal<-pal
      proxy %>% addPolygons(data=vineyards.shp, fillOpacity=0.9,color=("black"),,label = ~as.character(vineyard))
    } 
  })  
  
  # TO display in leaflet with labels try:
  #leaflet(vineyards) %>% addTiles() %>%
  #  addPolygons(lng = ~Long, lat = ~Lat,label = ~as.character(Name))
  
  # Calculate value for mapclick - Display lat lon of mouse click
  observe({  #Observer to show Popups on click
    click <- input$map_click
    if (!is.null(click)) {
      lon=click$lng
      lat=click$lat
      cellnumber<-lonlat_to_cellnumber(chosenlayer(),lon,lat)
      value = as.numeric(chosenlayer()[cellnumber])
      show_popup(cellnumber,chosenvar(),value,lon,lat)
      # Plot hist
      show_hist(getValues(chosenlayer()),value)
      show_timeseries(chosendata(),cellnumber)
    }
  })
  
  # REACTIVE FUNCTIONS
  # Set variable to be mapped - dependent on tab panel selected - CHECK if works without!
  chosenvar <- reactive({
    if (input$tabVariables == "temp") input$variable1 else 
      if (input$tabVariables == "frost") input$variable2 else 
        if (input$tabVariables == "flower") input$variable3 else 
          if (input$tabVariables == "terrain") input$variable4 
  })
  
  chosenlayer <- reactive({
    req(input$county)
    req(input$year)
    req(chosenvar())
    getmap(input$county,chosenvar(),input$year)
    }) 
  
  chosendata<-reactive({
    req(input$county)
    req(chosenvar())
    get_timeseries(input$county,chosenvar())
    }) 

  visibility<-reactive({
    req(input$visibility)
    input$visibility
  }) 
  
#################################################
# SEARCH MAPS PAGE
#################################################
# Modify slider input range according to variables chosen
  observe({ 
    minmax1 <- switch(input$search1,
                      "gdd10" = list(0,2000),
                      "gdd5" = list(0,3000),
                      "meant" = list(5,25),
                      "maxt" = list(20,50),
                      "mint" =list(-10,10),
                      "num20"=list(0,200),
                      "num25"=list(0,200),
                      "num30"=list(0,100),
                      "spfr"=list(0,150),
                      "autfr"=list(260,366),
                      "frfree"=list(0,366),
                      "elevation"=list(0,500),
                      "slope"=list(0,30),
                      "aspect"=list(0,360)
                   )
    minmax2 <- switch(input$search2,
                      "gdd10" = list(0,2000),
                      "gdd5" = list(0,3000),
                      "meant" = list(5,25),
                      "maxt" = list(20,50),
                      "mint" =list(-10,10),
                      "num20"=list(0,200),
                      "num25"=list(0,200),
                      "num30"=list(0,100),
                      "spfr"=list(0,150),
                      "autfr"=list(260,366),
                      "frfree"=list(0,366),
                      "elevation"=list(0,500),
                      "slope"=list(0,30),
                      "aspect"=list(0,360)
    )
    minmax3 <- switch(input$search3,
                      "gdd10" = list(0,2000),
                      "gdd5" = list(0,3000),
                      "meant" = list(5,25),
                      "maxt" = list(20,50),
                      "mint" =list(-10,10),
                      "num20"=list(0,200),
                      "num25"=list(0,200),
                      "num30"=list(0,100),
                      "spfr"=list(0,150),
                      "autfr"=list(260,366),
                      "frfree"=list(0,366),
                      "elevation"=list(0,500),
                      "slope"=list(0,30),
                      "aspect"=list(0,360)
    )
    updateSliderInput(session, "slider1", min = minmax1[1], max=minmax1[2])
    updateSliderInput(session, "slider2", min = minmax2[1], max=minmax2[2])
    updateSliderInput(session, "slider3", min = minmax3[1], max=minmax3[2])
    }) # observe 
 
  # Display basic map 
  output$searchmap<-renderLeaflet({
    leaflet() %>% setView(lng = -4.5, lat = 50.75, zoom = 8) %>% 
      addTiles()  %>%  
      clearImages() %>% 
      addScaleBar()
  })
  
  # Observe to change data if UI changes
  observeEvent(input$startsearch,{
    vars<-c(input$search1,input$search2,input$search3)
    mins<-c(min(input$slider1),min(input$slider2),min(input$slider3))
    maxs<-c(max(input$slider1),max(input$slider2),max(input$slider3))
    
    # Observe to add Vineyards
    observe({
      proxy <- leafletProxy("searchmap")
      # Remove any existing legend, and only if the legend is
      # enabled, create a new one.
      proxy %>% clearShapes()
      if (input$vineyards2){
        #pal<-pal
        proxy %>% addPolygons(data=vineyards.shp, fillOpacity=0.75,color=("black"),,label = ~as.character(vineyard))
      } 
    })  
    
    #searchmap<-create_search_map(vars,mins,maxs,input$searchcounty,input$searchyear)
    leafletProxy("searchmap") %>% 
      clearPopups() %>%
      clearImages() %>% 
      addRasterImage(create_search_map(vars,mins,maxs,input$searchcounty,input$searchyear), 
                     color=colorBin(c("grey","green"),domain=c(0,1),bins=2),opacity = visibility(),project=FALSE)   
    
    create_search_map_text(vars,mins,maxs,input$searchcounty,input$searchyear)
  })

#################################################
# FUNCTIONS CALLED BY SERVER
#################################################
  
# Function that will return file name of raster according to input choices
getmap<-function(county,var,year){ 
  if (county=="cornwall") dir_data<-dir_cornwall else
    if (county=="devon") dir_data<-dir_devon else
      if (county=="dorset") dir_data<-dir_dorset else
        if (county=="somerset") dir_data<-dir_somerset
        
        if (year=="2009")  dir_data2<-paste(dir_data,"2009/",sep="")   
        if (year=="2010")  dir_data2<-paste(dir_data,"2010/",sep="")   
        if (year=="2011")  dir_data2<-paste(dir_data,"2011/",sep="")   
        
        if (var=="elevation") raster(paste(dir_data,"dem.tif",sep="")) else 
          if (var=="slope") raster(paste(dir_data,"slope.tif",sep="")) else
            if (var=="aspect") raster(paste(dir_data,"aspect.tif",sep="")) else
              
              if (var=="gdd10") raster(paste(dir_data2,"gdd10-",year,".tif",sep="")) else 
                if (var=="gdd5") raster(paste(dir_data2,"gdd5-",year,".tif",sep="")) else
                  if (var=="meant") raster(paste(dir_data2,"meant-",year,".tif",sep="")) else
                    if (var=="mint") raster(paste(dir_data2,"mint-",year,".tif",sep="")) else 
                      if (var=="maxt") raster(paste(dir_data2,"maxt-",year,".tif",sep="")) else
                        if (var=="num20") raster(paste(dir_data2,"days20-",year,".tif",sep="")) else
                          if (var=="num25") raster(paste(dir_data2,"days25-",year,".tif",sep="")) else 
                            if (var=="num30") raster(paste(dir_data2,"days30-",year,".tif",sep="")) else
                              
                              if (var=="spfr") raster(paste(dir_data2,"spfrost-",year,".tif",sep="")) else 
                                if (var=="autfr") raster(paste(dir_data2,"autfrost-",year,".tif",sep="")) else
                                  if (var=="frfree") raster(paste(dir_data2,"frostfree-",year,".tif",sep="")) else
                                    
                                    if (var=="flcond") raster(paste(dir_data2,"flcond.tif",sep=""))
} # getmap

# Set colour scheme according to variable mapped
mapcolour <- function(var) {
  if (var == "gdd10") colorBin("YlOrRd",domain=c(0,1800),bins=9,pretty=TRUE) else    # 300 to 1800
    if (var == "gdd5") colorBin("YlOrRd",domain=c(1000,3000),bins=9,pretty=TRUE) else      #1000 to 3000
      if (var == "maxt") colorBin("YlOrRd",domain=c(10,50),bins=8,pretty=TRUE) else   #15 to 40
        if (var == "mint") colorBin(c("blue","white"),domain=c(-10,12),bins=10,pretty=TRUE) else    # -10 to 10
          if (var == "meant") colorBin("YlOrRd",domain=c(8,24),bins=7,pretty=TRUE) else   # 10-20
            if (var == "num20") colorBin(c("white","yellow","red"),domain=c(0,200),bins=10,pretty=TRUE) else  # 0 to 200
              if (var == "num25") colorBin(c("yellow","red"),domain=c(0,100),bins=10,pretty=TRUE) else #0 to 100
                if (var == "num30") colorBin("OrRd",domain=c(0,50),bins=10,pretty=TRUE) else # 0 to 30
                  if (var == "spfr") colorBin("Blues",domain=c(1,150),bins=7,pretty=TRUE) else      # range 1 to 150 = 150 
                    if (var == "autfr") colorBin("Blues",domain=c(370,240),bins=7,pretty=TRUE) else    # range 240 to 366 = 126
                      if (var == "frfree") colorBin("BuGn",domain=c(90,370),bins=9,pretty=TRUE) else   # 90-366
                        if (var == "flcond") colorBin("PuRd",domain=c(300,1800),bins=7,pretty=TRUE) else
                          if (var == "elevation") colorBin(c("blue","green"),domain=c(0,450),bins=9,pretty=TRUE) else # 1 to 450
                            if (var == "slope") colorBin("Purples",domain=c(0,27),bins=9,pretty=TRUE) else  # 1-30
                              if (var == "aspect") colorBin("PRGn",domain=c(0,360),bins=7,pretty=TRUE)     # 1 to 360
} # mapcolour

# Returns filename of matrix var.m[cellnumber,year] holding values of variable for each cell and year
get_timeseries<-function(county,var){ 
  paste(dir_shinydata,county,"/",var,"-timeseries.R",sep="")
} 

#get_vineyards<-function(county){ 
#  load(paste(dir_shinydata,county,"/vineyards.RData",sep=""))
#  return(vineyards.epsg)
#} 

# Display popup on mouse click of lon,lat and cell value
show_popup <- function(cellnumber,var,value,lon=NULL, lat=NULL) { 
  #print(value)
  content <- paste("Lon=", round(lon, 2),"; Lat=", round(lat, 2),"; ", var,"=",round(value,0),sep="")
  proxy <- leafletProxy("map")
  #add Popup
  proxy %>% clearPopups() %>% addPopups(lon, lat, popup = content)
  output$celltext<-renderText({
    paste("Lat: ",round(lat, 2),"; Lon: ",round(lon,2), "; Value: ",round(value,0),sep="")
  })
}

# Display freq (as %of all cells) histogram of all cell values and indicate chosen cell value
show_hist <- function(data.v,value){
  output$hist<-renderPlot( plothist(data.v) + geom_vline(xintercept = value,color = "red", size=1) )
}

show_timeseries<-function(filename,cellnumber){
  # plot timeseries
  load(filename) # load var.m
  values<-var.m[cellnumber,]
  years<-c(2009,2010,2011)
  y.labels<-c("2009","2010","2011")
  df = data.frame(years,values) 
  output$timeseriesplot<-renderPlot(ggplot(df, aes(x = df$years, y = df$values))  + 
                                      geom_line(color="red") +
                                      scale_y_continuous("Value") +
                                      scale_x_continuous("Years")
  ) # renderPlot
}

# Function to create raster of cells meeting search criteria

create_search_map<-function(vars,mins,maxs,searchcounty,searchyear) {
  r<-getmap(searchcounty,"elevation",searchyear)
  r<-calc(r,fun=function(x){ifelse(is.na(x),NA,1)}) # set default - all land cells selected
  if (vars[1]!="") {
    new.r<-calc(getmap(searchcounty,vars[1],searchyear),fun=function(x){ifelse(x>=mins[1] & x<=maxs[1],1,0)})      
    r<-mask(r,new.r,maskvalue=0,updatevalue=0)  
  }
  if (vars[2]!="") {
    new.r<-calc(getmap(searchcounty,vars[2],searchyear),fun=function(x){ifelse(x>=mins[2] & x<=maxs[2],1,0)})      
    r<-mask(r,new.r,maskvalue=0,updatevalue=0) 
  }
  if (vars[3]!="") {
    new.r<-calc(getmap(searchcounty,vars[3],searchyear),fun=function(x){ifelse(x>=mins[3] & x<=maxs[3],1,0)})      
    r<-mask(r,new.r,maskvalue=0,updatevalue=0)
  }
return(r)
  
}# create_search_map

create_search_map_text<-function(vars,mins,maxs,searchcounty,searchyear) {
  text1<-paste(vars[1],"between",mins[1],"and",maxs[1])
  if (vars[2]!="") text2<-paste("; ",vars[2],"between",mins[2],"and",maxs[2]) else text2<-""
  if (vars[3]!="") text3<-paste("; ",vars[3],"between",mins[3],"and",maxs[3]) else text3<-""
  
  output$searchmaptext<-renderText(
    paste("Highlighted locations where:",text1,text2,text3)
  )
} # end search_map_text

# Calculate cell number from lon lat values
lonlat_to_cellnumber<-function(r,x,y,zone=0){ #http://stackoverflow.com/questions/18639967/converting-latitude-and-longitude-points-to-utm
xy <- data.frame(ID = 1:length(x), X = x, Y = y)
coordinates(xy) <- c("X", "Y")
proj4string(xy) <- CRS("+init=epsg:4326")  ## lat lon
xy <- spTransform(xy, CRS(epsg)) # convert to raster xy
cell <- cellFromXY(r, c(xy$X,xy$Y))
return(cell)
}

roundUpNice <- function(x, nice=c(1,2,4,5,6,8,10)) {
  if(length(x) != 1) stop("'x' must be of length 1")
  10^floor(log10(x)) * nice[[which(x <= 10^floor(log10(x)) * nice)[[1]]]]
}

# Plot histogram
plothist<-function(values){
  sel<-which(!is.na(values))
  values<-values[sel]
  numcells<-length(values)
  mnx<-min(values)
  mxx<-max(values)
  #f <- hist(values, maxpixels=length(values),breaks=50) # calculate from all cells of raster
  f <- hist(values,breaks=50) 
  #abline(v=100, col="blue",lwd=3)
  dat <- data.frame(counts= ((f$counts/numcells)*100),breaks = f$mids)
  ggplot(dat, aes(x = breaks, y = counts, fill =counts)) + ## Note the new aes fill here
    geom_bar(stat = "identity",alpha = 0.8,fill="blue")+
    xlab("value")+ ylab("%")+
    scale_x_continuous(breaks = seq(0,mxx,(roundUpNice(mxx/5))),
                       labels = seq(0,mxx,(roundUpNice(mxx/5))) )
}
  
} # server

shinyApp(ui=ui, server=server) 


# Colour Functions
# Create a continuous palette function
#pal <- colorNumeric(palette = "Blues",domain = countries$gdp_md_est)
#qpal <- colorQuantile("Blues", countries$gdp_md_est, n = 7)
##binpal <- colorBin("Blues", countries$gdp_md_est, 6, pretty = FALSE)
#factpal <- colorFactor(topo.colors(5), countries$category)
#color = ~qpal(gdp_md_est)


# NOTES - FUTURE OPPS
# NEXT
  # Add timescale selector: Indivudal Year OR Risk timescale 
  # Relocate scale
# Use conditional panel to provide relevatn info for viewed plot (eg on GDD etc )
# Display values with mouse click - is it poss for multiple layers of info?? - use observeEvent(input$MAPID_click)
# Multiple layer display??
# Find locations UI screen - set gdd and other risk ranges to locate
# Add vineyard locations shape file???

# Side Panel
# Display summary table of ??


#### Add little shape county map to allow counties to be identified ?
#sidebarPanel()
#mainPanel(
 # tabsetPanel(
  #tabPanel("title",output,actionbutton etc)
#)
#)
#navbarPage - whole page layers 

#tabPanel("title",)
