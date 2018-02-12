# Install dev version of leaflet for R - required
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

# Define directories and files
# To use with uploaded app set work directory to directory of app and modify data locations accordingly
# setwd( "/Users/yimizhao/Desktop/Study/Data_science/Shiny/census_app")
# 
root<-"/Volumes/Big Sam/"; in.root<-"/Volumes/Big Sam/wineclim/"
dir_shinydata<-paste(root,"wineclim/",sep="")
dir_rasters<-paste(dir_shinydata,"rasters/",sep="")
dir_timeseries<-paste(dir_shinydata,"timeseries/",sep="")

# Link to vineyards file
load(paste(dir_shinydata,"vineyards.RData",sep="") ) # loads vineyards.shp

res<-c(100,100)

t1<-theme(                              
  plot.background = element_blank(), 
  panel.grid.major = element_blank(), 
  panel.grid.minor = element_blank(), 
  panel.border = element_blank(), 
  panel.background = element_blank(),
  axis.line = element_line(size=1)
)

#################################################
# FUNCTIONS CALLED BY SERVER
#################################################

# Function that will return file name of raster according to input choices
getmap<-function(county,var,year){ 
  if (var %in% c("elevation","slope","aspect")) paste(dir_rasters,var,"_",county,".tif",sep="") else 
    paste(dir_rasters,var,"_",year,"_",county,".tif",sep="")
} # getmap

# Set colour scheme according to variable mapped
# mnmx = 2 element array holding min and max
mapcolour <- function(var,mnmx) {
  mnmx<-unlist(mnmx)
  switch(var,
         "gdd10_gs" = colorQuantile("YlOrRd",domain=mnmx,n=10),
         "gdd5_gs" = colorBin("YlOrRd",domain=mnmx,bins=9,pretty=TRUE) ,
         "tmean_gs" = colorBin("YlOrRd",domain=mnmx,bins=7,pretty=TRUE) ,
         "tmax_year" = colorBin("YlOrRd",domain=mnmx,bins=8,pretty=TRUE),
         "tmin_year" = colorBin(c("Blues"),domain=mnmx,bins=10,pretty=TRUE),
         "t20_gsdays"= colorBin(c("white","red"),domain=mnmx,bins=10,pretty=TRUE),
         "t25_gsdays"= colorBin(c("white","red"),domain=mnmx,bins=10,pretty=TRUE),
         "t30_gsdays"= colorBin(c("white","red"),domain=mnmx,bins=10,pretty=TRUE),
         "lastspfr_doy"=colorBin("Blues",domain=mnmx,bins=7,pretty=TRUE) ,
         "firstautfr_doy"=colorBin("Blues",domain=mnmx,bins=7,pretty=TRUE),
         "frostfree_days"=colorBin("BuGn",domain=mnmx,bins=9,pretty=TRUE),
         "fl_tmean"=colorBin("YlOrRd",domain=mnmx,bins=7,pretty=TRUE),
         "fl_numday"=colorBin("Greens",domain=mnmx,bins=7,pretty=TRUE),
         "elevation"=colorBin(c("blue","green"),domain=mnmx,bins=9,pretty=TRUE),
         "slope"=colorBin("Purples",domain=mnmx,bins=7,pretty=TRUE),
         "aspect"=colorBin("PRGn",domain=mnmx,bins=7,pretty=TRUE) 
  )
} # mapcolour

# Returns filename of matrix var.m[cellnumber,year] holding values of variable for each cell and year
# or "" if not approoriate data (eg terrain)
get_timeseries<-function(county,var){ 
  if (var %in% c("elevation","slope","aspect")) "" else 
    paste(dir_timeseries,var,"_",county,"_timeseries.R",sep="")  
} 



# Return min max values for different variables
set_min_max<-function(var){
  switch(var,
         "gdd10_gs" = list(0,1500),
         "gdd5_gs" = list(0,2000),
         "tmean_gs" = list(10,20),
         "tmax_year" = list(20,45),
         "tmin_year" =list(-10,15),
         "t20_gsdays"=list(0,200),
         "t25_gsdays"=list(0,100),
         "t30_gsdays"=list(0,100),
         "lastspfr_doy"=list(0,150),
         "firstautfr_doy"=list(260,366),
         "frostfree_days"=list(0,366),
         "fl_tmean"=list(5,40),
         "fl_numday"=list(0,14),
         "elevation"=list(0,500),
         "slope"=list(0,35),
         "aspect"=list(0,360)
  )
} # set_min_max

# Return long text label for variable
var_longlabel<-function(var){
  switch(var,
         "gdd10_gs" = "Growing degree days (base 10C)",
         "gdd5_gs" = "Growing degree days (base 5C)",
         "tmean_gs" = "Mean Growing Season Temperature (C)",
         "tmax_year" ="Max Temperature (C)",
         "tmin_year" ="Min Temperature (C)",
         "t20_gsdays"="No. days where max temperature>20C",
         "t25_gsdays"="No. days where max temperature>25C",
         "t30_gsdays"="No. days where max temperature>30C",
         "lastspfr_doy"="Last spring frost (day of year)",
         "firstautfr_doy"="First autumn frost (day of year)",
         "frostfree_days"="Frost free days in year",
         "fl_tmean"="Mean flowering temperature",
         "fl_numday"="No. of good flowering days",
         "elevation"="Elevation (metres)",
         "slope"="Slope (degrees)",
         "aspect"="Aspect (degrees)"
  )
} # var_longlabel

# Return long text label for variable
var_unit<-function(var){
  switch(var,
         "gdd10_gs" = "degree days",
         "gdd5_gs" = "degree days",
         "tmean_gs" = "degrees C",
         "tmax_year" ="degrees C",
         "tmin_year" ="degrees C",
         "t20_gsdays"="days",
         "t25_gsdays"="days",
         "t30_gsdays"="days",
         "lastspfr_doy"="day of year",
         "firstautfr_doy"="day of year",
         "frostfree_days"="days",
         "fl_tmean"="degrees C",
         "fl_numday"="days",
         "elevation"="metres",
         "slope"="degrees",
         "aspect"="degrees"
  )
} # var_unit


# Function to create raster of cells meeting search criteria

create_search_map<-function(vars,mins,maxs,searchcounty,searchyear) {
  r<-raster(getmap(searchcounty,"elevation",searchyear))
  r<-calc(r,fun=function(x){ifelse(is.na(x),NA,1)}) # set default - all land cells selected
  if (vars[1]!="") {
    new.r<-raster(getmap(searchcounty,vars[1],searchyear))
    new.r<-calc(new.r,fun=function(x){ifelse(x>=mins[1] & x<=maxs[1],1,0)})      
    r<-mask(r,new.r,maskvalue=0,updatevalue=0)  
  }
  if (vars[2]!="") {
    new.r<-raster(getmap(searchcounty,vars[2],searchyear))
    new.r<-calc(new.r,fun=function(x){ifelse(x>=mins[2] & x<=maxs[2],1,0)})      
    r<-mask(r,new.r,maskvalue=0,updatevalue=0) 
  }
  if (vars[3]!="") {
    new.r<-raster(getmap(searchcounty,vars[3],searchyear))
    new.r<-calc(new.r,fun=function(x){ifelse(x>=mins[3] & x<=maxs[3],1,0)})      
    r<-mask(r,new.r,maskvalue=0,updatevalue=0)
  }
  return(r)
  
}# create_search_map



# Calculate cell number from lon lat values
lonlat_to_cellnumber<-function(r,x,y,zone=0){ #http://stackoverflow.com/questions/18639967/converting-latitude-and-longitude-points-to-utm
  xy <- data.frame(ID = 1:length(x), X = x, Y = y)
  coordinates(xy) <- c("X", "Y")
  proj4string(xy) <- CRS("+init=epsg:4326")  ## lat lon
  xy <- spTransform(xy, CRS("+init=epsg:3857")) # convert to raster xy
  cell <- cellFromXY(r, c(xy$X,xy$Y))
  print(cell)
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
  mnx<-min(values,na.rm=TRUE)
  mxx<-max(values,na.rm=TRUE)
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
                                   choices=  c("1996" = "1996","1997" = "1997","1998" = "1998","1999" = "1999","2000" = "2000",
                                               "2001" = "2001","2002" = "2002","2003" = "2003","2004" = "2004","2005" = "2005",
                                               "2006" = "2006","2007" = "2007","2008" = "2008","2009" = "2009","2010" = "2010",
                                               "2011" = "2011","2012" = "2012") ,
                                   selected="2012",selectize=TRUE  ) 
                ), # box
                tabBox(id="tabVariables", width=8 ,
                      
                	tabPanel("Growing season temperature",
                    selectInput("variable1","Select data", c("Choose one"="","Degree days (base 10C)" = "gdd10_gs","Degree days (base 5C)" = "gdd5_gs", 
                                                             "Mean temperature" = "tmean_gs", "Max temperature"="tmax_year","Min temperature"="tmin_year",
                                                             "No. days > 20C"="t20_gsdays","No. days > 25C"="t25_gsdays","No. days > 30C"="t30_gsdays") ,
                                selected="",selectize=TRUE  ),
                    tags$p("Growing season: 1 April to 31 October."),
                    value="temp" 
                  ) , # tabPanel
                    
                  tabPanel("Frost Risks",  
                  	 selectInput("variable2","Select data", c("Choose one"="","Last spring frost" = "lastspfr_doy","First autumn frost" = "firstautfr_doy", 
                  	                                          "Frost free period" = "frostfree_days") ,
                  	             selected="",selectize=TRUE  ),
                  	  tags$p("Frost event defined as < 1C. Dates given as days of year (1 to 366)"),
                  	 value="frost"
             		  ), # tabPanel
             		
             		  tabPanel("Flowering Risks",  
                  	 selectInput("variable3","Select data", c("Choose one"="","Mean flowering temperature" = "fl_tmean",
                  	                                          "Good flowering days" = "fl_numday") ,
                  	             selected="",selectize=TRUE  ),
                  	  tags$p("Flowering: 22 June-5 July. Good flowering day: mean temperature>15 C"),
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
                    htmlOutput("celltext") ),
                box(title="% Distribution of map values", width=3,
                    plotOutput("hist",height="200px") ),
                box(title="Historic variation", width=3,
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
                                        choices=  c("1996" = "1996","1997" = "1997","1998" = "1998","1999" = "1999","2000" = "2000",
                                                    "2001" = "2001","2002" = "2002","2003" = "2003","2004" = "2004","2005" = "2005",
                                                    "2006" = "2006","2007" = "2007","2008" = "2008","2009" = "2009","2010" = "2010",
                                                    "2011" = "2011","2012" = "2012") ,
                                        selected="2012",selectize=TRUE ) 
                            )
                        ), # fluidRow
              fluidRow(height=60, h4("  Choose search variables & range:"),
                        box(width=4,
                            selectInput("search1",NULL, c("None chosen"="","Degree days (base 10C)" = "gdd10_gs","Degree days (base 5C)" = "gdd5_gs", 
                                                                   "Mean temperature" = "tmean_gs", "Max temperature"="tmax_year","Min temperature"="tmin_year",
                                                                   "No. days > 20C"="t20_gsdays","No. days > 25C"="t25_gsdays","No. days > 30C"="t30_gsdays",
                                                                   "Last spring frost" = "lastspfr_doy","First autumn frost" = "firstautfr_doy", "Frost free period" = "frostfree_days",
                                                                   "Mean flowering temp" = "fl_tmean","Good flowering days"="fl_numday",
                                                                   "Elevation" = "elevation","Slope" = "slope", "Aspect" = "aspect"), 
                                                                    selected="",selectize=FALSE  ) ,
                            sliderInput("slider1", label=NULL, min = 0, 
                                        max = 100, value = c(20,50))
                        ) ,  # box
                       box( width=4,
                            selectInput("search2",NULL, c("None"="","Degree days (base 10C)" = "gdd10_gs","Degree days (base 5C)" = "gdd5_gs", 
                                                                   "Mean temperature" = "tmean_gs", "Max temperature"="tmax_year","Min temperature"="tmin_year",
                                                                   "No. days > 20C"="t20_gsdays","No. days > 25C"="t25_gsdays","No. days > 30C"="t30_gsdays",
                                                                   "Last spring frost" = "lastspfr_doy","First autumn frost" = "firstautfr_doy", "Frost free period" = "frostfree_days",
                                                                   "Mean flowering temp" = "fl_tmean","Good flowering days"="fl_numday",
                                                                   "Elevation" = "elevation","Slope" = "slope", "Aspect" = "aspect"), 
                                        selected="",selectize=FALSE  ) ,
                           sliderInput("slider2", label=NULL, min = 0, 
                                       max = 100, value = c(20,50))
                       ) , # box
                       box( width=4,
                            selectInput("search3",NULL,c("None"="","Degree days (base 10C)" = "gdd10_gs","Degree days (base 5C)" = "gdd5_gs", 
                                                                   "Mean temperature" = "tmean_gs", "Max temperature"="tmax_year","Min temperature"="tmin_year",
                                                                   "No. days > 20C"="t20_gsdays","No. days > 25C"="t25_gsdays","No. days > 30C"="t30_gsdays",
                                                                   "Last spring frost" = "lastspfr_doy","First autumn frost" = "firstautfr_doy", "Frost free period" = "frostfree_days",
                                                                   "Mean flowering temp" = "fl_tmean","Good flowering days"="fl_numday",
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
  
  # Display popup on mouse click of lon,lat and cell value
  show_popup <- function(cellnumber,var,value,lon=NULL, lat=NULL) { 
    #print(value)
    content <- paste("Lon=", round(lon, 2),"; Lat=", round(lat, 2),"; Value=",round(value,0)," ",var_unit(var),sep="")
    proxy <- leafletProxy("map")
    #add Popup
    proxy %>% clearPopups() %>% addPopups(lon, lat, popup = content)
    output$celltext<-renderUI({
      line1<-paste("Latitude:",round(lat, 2),"; Longitude:",round(lon,2))
      line2<-paste(var_longlabel(var)," = ",round(value,0)," ",var_unit(var),sep="")
      HTML(paste(line1, line2, sep = '<br/>'))
    })
  }
  
  # Display freq (as %of all cells) histogram of all cell values and indicate chosen cell value
  show_hist <- function(data.v,value){
    output$hist<-renderPlot( plothist(data.v) + geom_vline(xintercept = value,color = "red", size=1) )
  }
  
  show_timeseries<-function(filename,cellnumber){
    # plot timeseries if appropriate or blank figure if not
    if (filename==""){
      years<-c() 
      values<-c()
    }
    if (filename!=""){
      load(filename) # load var.m
      values<-var.m[cellnumber,]
      years<-c(1983:2013)
      #y.labels<-c("2012")
    }
    values.df = data.frame(years,values) 
    print(paste("Number of NA =",length(which(is.na(values.df)))))
    maxy<-max(values.df$values,na.rm=TRUE); miny<-min(values.df$values,na.rm=TRUE)
    maxx<-max(values.df$years,na.rm=TRUE); minx<-min(values.df$years,na.rm=TRUE)
    print(paste("X:",minx,maxx,"Y:",miny,maxy))
    output$timeseriesplot<-renderPlot(ggplot(values.df, aes(x = values.df$years, y = values.df$values))  + 
                                        geom_line(color="red") +
                                        scale_y_continuous(limits=c(miny,maxy),"Value") +
                                        scale_x_continuous(limits=c(minx,maxx),"Years")
    ) # renderPlot
  }
  
  create_search_map_text<-function(vars,mins,maxs,searchcounty,searchyear) {
    text1<-paste(vars[1],"between",mins[1],"and",maxs[1])
    if (vars[2]!="") text2<-paste("; ",vars[2],"between",mins[2],"and",maxs[2]) else text2<-""
    if (vars[3]!="") text3<-paste("; ",vars[3],"between",mins[3],"and",maxs[3]) else text3<-""
    
    output$searchmaptext<-renderText(
      paste("Highlighted locations where:",text1,text2,text3)
    )
  } # end search_map_text
  
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
    paste("Map showing ",var_longlabel(chosenvar()),"for",input$county, "in",input$year,sep=" ")
  })
  
  # Observe to change map data if UI changes
  observe({
    leafletProxy("map") %>% 
      clearPopups() %>%
      clearImages() %>% 
      addRasterImage(chosenlayer(),color=mapcolour(chosenvar(),set_min_max(chosenvar()) ), 
                     opacity = visibility(),project=FALSE)  
  })
  
  # Observe to change legend if requested
  observe({
    proxy <- leafletProxy("map")
    # Remove any existing legend, and only if the legend is
    # enabled, create a new one.
    proxy %>% clearControls()
    if (input$legend){
      #pal<-pal
      proxy %>% addLegend(position = "bottomright",pal=mapcolour(chosenvar(),set_min_max(chosenvar())),values=values(chosenlayer()) )
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
      proxy %>% addPolygons(data=vineyards.shp, fillOpacity=0.9,color=("black"),label = ~as.character(vineyard))
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
      # Plot timeseries if appropriate
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
    raster(getmap(input$county,chosenvar(),input$year))
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
    minmax1 <- set_min_max(input$search1)
    minmax2 <- set_min_max(input$search2)
    minmax3 <- set_min_max(input$search3)
    updateSliderInput(session, "slider1", min = minmax1[[1]], max=minmax1[[2]])
    updateSliderInput(session, "slider2", min = minmax2[[1]], max=minmax2[[2]])
    updateSliderInput(session, "slider3", min = minmax3[[1]], max=minmax3[[2]])
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
        proxy %>% addPolygons(data=vineyards.shp, fillOpacity=0.75,color=("black"),label = ~as.character(vineyard))
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
