## Notes and trial code for Shiny
# APP development
#Function
latlon_to_rasterxy<-function(x,y,zone=0){
  xy <- data.frame(ID = 1:length(x), X = x, Y = y)
  coordinates(xy) <- c("X", "Y")
  proj4string(xy) <- CRS("+init=epsg:4326")  ## lat lon
  res <- spTransform(xy, CRS(epsg))
  return(as.data.frame(res))
}



mapdata<-reactive({raster(file=getmap(input$variable))}) 
observe({
  newmap<-mapdata()
  leafletProxy("map") %>% 
    addRasterImage(getmap(input$variable),  opacity = 0.7,project=FALSE)
})    


observe(input$variable {
leafletProxy("map") %>% removeRasterImage()
leafletProxy("map") %>% addRasterImage(getmap(input$variable))
})



leafletProxy("map", data = filteredData()) %>%
  clearShapes() %>%
  addCircles(radius = ~10^mag/10, weight = 1, color = "#777777",
             fillColor = ~pal(mag), fillOpacity = 0.7, popup = ~paste(mag)

             
  # to add legend           
             observe({
               proxy <- leafletProxy("map")
               # Remove any existing legend, and only if the legend is
               # enabled, create a new one.
               proxy %>% clearControls()
               if (input$legend) proxy %>% addLegend(position = "bottomright")
             })             
             
             
             
#  output$table<-renderTable({summary(mapdata())})


# sliderInput(inputId = "",label = "Select budbreak time (day of year)", value= 180, min=  , max = ) ,

output$plot<-renderPlot({image})
plot.r<-reactive({})
isolate
observeEvent(input$replot{})
observe({})
data<-eventReactive(input$go,{do something/select data})
output$plot<-renderPlot({plot(data)})

rv$data<-reactiveValues(arguments = )
observeEvetnt
output$plot<-renderPlot({plot(rv$data)})

ui.select<-fluidPage(
  selectInput(inputId = "county",label="Select Geographical area", 
              choices=  c("Cornwall" = "cornwall",
                          "Devon" = "devon",
                          "Somerset" = "other")  )
  selectInput(inputId = "risk",label="Select Risk", 
              choices=  c("GDD" = "gdd",
                          "Frost" = "frost",
                          "Flowering" = "floer") , )
  
  
  
  
  
  input$MAPID_click is an event that is sent when the map background or basemap is clicked. The value is a list with lat and lng.
  
  input$MAPID_bounds provides the latitude/longitude bounds of the currently visible map area; the value is a list() that has named elements north, east, south, and west.
  
  input$MAPID_zoom is an integer that indicates the zoom level.
  
  
  # Show a popup at the given location
  showZipcodePopup <- function(zipcode, lat, lng) {
    selectedZip <- allzips[allzips$zipcode == zipcode,]
    content <- as.character(tagList(
      tags$h4("Score:", as.integer(selectedZip$centile)),
      tags$strong(HTML(sprintf("%s, %s %s",
                               selectedZip$city.x, selectedZip$state.x, selectedZip$zipcode
      ))), tags$br(),
      sprintf("Median household income: %s", dollar(selectedZip$income * 1000)), tags$br(),
      sprintf("Percent of adults with BA: %s%%", as.integer(selectedZip$college)), tags$br(),
      sprintf("Adult population: %s", selectedZip$adultpop)
    ))
    leafletProxy("map") %>% addPopups(lng, lat, content, layerId = zipcode)
  }
  
  # When map is clicked, show a popup with city info
  observe({
    leafletProxy("map") %>% clearPopups()
    event <- input$map_shape_click
    if (is.null(event))
      return()
    
    isolate({
      showZipcodePopup(event$id, event$lat, event$lng)
    })
  
    observe({
      click<-input$map_marker_click
      if(is.null(click))
        return()
      text<-paste("Lattitude ", click$lat, "Longtitude ", click$lng)
      text2<-paste("You've selected point ", click$id)
      map$clearPopups()
      map$showPopup( click$lat, click$lng, text)
      output$Click_text<-renderText({
        text2
      })
      
      
      
      
      output$location <- renderPrint({
        validate(need(input$map_click, FALSE))
        str(input$map_click)
      })
      
      input$map_mouseover
      
      
      value<-extract(chosenlayer(),cell,method="simple")
      
      #Translate Lat-Lon to cell number using the unprojected raster
      cell <- cellFromXY(chosenlayer(), c(x, y))
      if (!is.na(cell)) {#If the click is inside the raster...
        xy <- xyFromCell(chosenlayer(), cell) #Get the center of the cell
        x <- xy[1]
        y <- xy[2]
        #Get row and column, to print later
        rc <- rowColFromCell(lldepth, cell)
        #Get value of the given cell
        val = depth[cell]
        content <- paste0("X=",rc[2],
                          "; Y=",rc[1],
                          "; Lon=", round(x, 5),
                          "; Lat=", round(y, 5),
                          "; Value=", round(val, 1), " m")