# Plot uk vineyard maps

# Vineyard hectares
plot_fig1a<-function(ukstats) {
  ggplot(ukstats, aes(x = ukstats$Year, y = ukstats$Total_Vines))  + 
    geom_line(color="red") +
    scale_y_continuous(limits=c(0,2000),"Hectares of Vines", breaks=c(0,500,1000,1500,2000) ) +
    scale_x_continuous(limits=c(1990,2015),"Year")
  
}

# Mean size of UK vineyards
plot_fig1b<-function(ukstats){
  ggplot(ukstats, aes(x = ukstats$Year, y = ukstats$Mean_Size_Vineyards))  + 
    geom_line(color="red") +
    scale_y_continuous(limits=c(0,4),"Mean size of Vineyards (Ha)") +
    scale_x_continuous(limits=c(1990,2015),"Year")
}

# Production bottles
plot_fig1c<-function(ukstats) {
  ggplot(ukstats, aes(x = ukstats$Year, y = ukstats$Production_Bottles))  + 
    geom_line(color="red") +
    scale_y_continuous(limits=c(0,8),"Production (million bottles)") +
    scale_x_continuous(limits=c(1990,2015),"Year")
  
}

# Timeseries plot of south-west mean growing season temperature
plot_fig2a<-function(swclimate){
  ggplot(swclimate, aes(x = swclimate$Year, y = swclimate$Temp))  + 
    geom_line(color="red") +
    scale_y_continuous(limits=c(11,15),"Mean Growing Season Temperature") +
    scale_x_continuous(limits=c(1920,2015),breaks=c(1920,1930,1940,1950,1960,1970,1980,1990,2000,2010,2020),"Year")
}

# Yield per Ha by Year
plot_fig2b<-function(ukstats) {
  ggplot(ukstats, aes(x = ukstats$Year, y = ukstats$Yield_per_Ha))  + 
    geom_line(color="red") +
    scale_y_continuous("Yield per Hectare of Vines", breaks=c(0,10,20,30,40)) +
    scale_x_continuous(limits=c(1990,2015),"Year")
}

# Plot yield per Ha vs SE mean GS temp
plot_fig2c<-function(seclimate,ukstats){
  plotdata<-merge(ukstats,seclimate,by=c("Year"))
  ggplot(plotdata, aes(y = plotdata$Yield_per_Ha, x = plotdata$Temp))  + 
    geom_point(color="red") +
    scale_y_continuous(limits=c(0,50),"Yield per Hectare") +
    scale_x_continuous(limits=c(10,15),"Mean Growing Season Temperature")
}

website_theme = theme(panel.grid.major = element_line(size = 0.5, color = "grey"), 
                      axis.line = element_line(size = 0.7, color = "black"), 
                      legend.position = c(0.85, 0.7), text = element_text(size = 14),
                      panel.background = element_blank())

# Load data
swfilein<-paste(dir_shinydata,"about/SW_GS_MeanT_1910-2014.csv",sep="")
swclimate<-read.csv(swfilein,header=TRUE)
sefilein<-paste(dir_shinydata,"about/SE_GS_MeanT_1910-2014.csv",sep="")
seclimate<-read.csv(swfilein,header=TRUE)
filein<-paste(dir_shinydata,"about/vineyards_stats.csv",sep="")
ukstats<-read.csv(filein,header=TRUE)

# Plot figs
plot_fig1a(ukstats) + website_theme
plot_fig1b(ukstats) + website_theme
plot_fig1c(ukstats) + website_theme

plot_fig2a(swclimate)+website_theme # Trend in GS mean T for SW
plot_fig2b(ukstats) + website_theme # Yield per Ha by Year
plot_fig2c(seclimate,ukstats) + website_theme # Yield per Ha vs Mean GS temp in SE England




