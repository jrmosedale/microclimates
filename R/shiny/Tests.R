ui <- dashboardPage(
  dashboardHeader(title="Viticulture: climate risk maps", titleWidth=300),
  
  dashboardSidebar(
    sidebarMenu(
      menuItem("About", tabName="about", icon=icon("info")),
      menuItem("Explore Maps",  icon=icon("map"),
               collapsible =
                 menuSubItem("Summary stats (1983-2013)", tabName = "explore_summary"),
               menuSubItem("Individual years", tabName = "explore_year")
      ),
      menuItem("Search Maps", icon=icon("search"),    
               menuSubItem("Summary (1983-2013)", tabName = "search_summary"),
               menuSubItem("Individual years", tabName = "search_year")
      )
    )
  ),
  
  dashboardBody(
    tabItems(
      tabItem(tabName="about",
              h2("Something about this app")
      ),
      tabItem(tabName="explore_summary",
              fluidRow(height=30, h4("  Define area & time period"), 
                       h2("Explore_summary"),
                       box(width=4,
                           selectInput(inputId = "searchcounty",label="Geographical area", 
                                       choices=  c("Cornwall" = "cornwall","Devon" = "devon","Dorset" = "dorset", "Somerset" = "somerset") ,
                                       selected="cornwall",selectize=TRUE  ),
                           textOutput("text")
                       )
      ) ),
      tabItem(tabName="explore_year",
              h2("Explore years !")
      ),
      tabItem(tabName="search_summary",
              h2("search Summary")
      ),
      tabItem(tabName="search_year",
              h2("search year")
      )
    ) # tabItems
  ) # dashboardBody
) # dashboardPage

server<-function(input,output,session) {
output$text<-renderText(chosenvar())
  chosenvar <- reactive({
    req(input$county)
    input$county
  })
  }
shinyApp(ui, server)
