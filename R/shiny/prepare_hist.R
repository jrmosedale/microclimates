# Prepare files required for shinyapp
# Reproject  and save in appropriate directories 

library(rgdal)
library(rgeos)
library(raster)
library(ggplot2)

root<-"~/Documents/Exeter/Data2015/"; in.root<-"~/Documents/Exeter/Data2015/"

# Input dirs
dir_dem<-paste(in.root,"DEM/",sep="")
dir_terrain<-paste(in.root,"Terrain/",sep="")
dir_results<-paste(root,"Outputs/",sep="") # dir of results rasters for whole of SW

# Output dirs
dir_shinydata<-paste(root,"shinydata/",sep="")
dir_rasters<-paste(root,"shinydata/rasters/",sep="")
dir_timeseries<-paste(root,"shinydata/timeseries/",sep="")
dir_hist<-paste(root,"shinydata/hist/",sep="")

#dir_swfiles<-paste(root,"shinydata/swfiles/",sep="")
#dir_cornwall<-paste(root,"shinydata/cornwall/",sep="")
#dir_devon<-paste(root,"shinydata/devon/",sep="")
#dir_dorset<-paste(root,"shinydata/dorset/",sep="")
#dir_somerset<-paste(root,"shinydata/somerset/",sep="")

# Define projections
latlong = "+init=epsg:4326"
ukgrid <- "+init=epsg:27700"
new.crs<-"+init=epsg:3857"
#new.crs<-'+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext +no_defs'

# Histogram functions - 
roundUpNice <- function(x, nice=c(1,2,4,5,6,8,10)) {
  if(length(x) != 1) stop("'x' must be of length 1")
  10^floor(log10(x)) * nice[[which(x <= 10^floor(log10(x)) * nice)[[1]]]]
}

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
#####################################################################
# Read shape file for county boundary 
#####################################################################
hcounties.shp<- readOGR(dsn = "/Users/jonathanmosedale/Documents/Exeter/Data2015/OSdata/bdline_essh_gb/Data/Supplementary_Historical", layer = "Boundary-line-historic-counties_region")
hcounties.shp<-spTransform(hcounties.shp,new.crs) # reproject

# extract county polygons - use historic counties data =OSuk projection
sel <- which(hcounties.shp$Name == "Cornwall")
cornwall<-hcounties.shp[sel,]
sel <- which(hcounties.shp$Name == "Devon")
devon<-hcounties.shp[sel,]
sel <- which(hcounties.shp$Name == "Dorset")
dorset<-hcounties.shp[sel,]
sel <- which(hcounties.shp$Name == "Somerset")
somerset<-hcounties.shp[sel,]

# PLOT - overlay counties 
plot(devon)
plot(somerset,add=TRUE)
plot(cornwall,add=TRUE)
plot(dorset,add=TRUE)

# Link terrain map data to files
dem<-raster(paste(dir_dem,"dem.tif",sep=""))
slope<-crop(raster(paste(dir_terrain,"slope.tif",sep="")),dem)
aspect<-crop(raster(paste(dir_terrain,"aspect.tif",sep="")),dem)

######################################
# Extract county rastersfor each statistic
######################################

# Define result years to be processed
years<-c("2012")

# County list and label
county.txt<-c("cornwall","devon","dorset","somerset")
county<-list(cornwall,devon,dorset,somerset) 

# Results variables 
# Statnames used in results files from analyse_years_carson.R
statnames<-c("gdd10_gs","gdd5_gs","tmean_gs","tmin_year","tmax_year",
             "t20_gsdays","t25_gsdays","t30_gsdays", 
             "lastspfr_doy", "firstautfr_doy", "frostfree_days",
             "fl_tmean", "fl_numday")
# var names used in  shiny app
vars<-c("gdd10","gdd5","maxt","mint","meant","num20","num25","num30", "spfr","autfr","frfree")
labs<-c("Degree days","Degree days","degrees C","degrees C","degrees C","Number of days","Number of days","Number of days",
        "Day of Year (1:366)","Day of Year (1:366)","Number of days")

# First set raster template in new projection for whole results area based on dem
  e<-extent(dem)
  template.r<-projectExtent(dem,crs=new.crs)
  res(template.r)<-c(100,100)
  print(template.r)

# 1 Reproject terrain data then mask to extract county rasters
for (c in 1:length(county)){
  print (county.txt[c])
  
  dem.county<-raster(paste(dir_rasters,"elevation_",county.txt[c],".tif",sep=""))
  slope.county<-raster(paste(dir_rasters,"slope_",county.txt[c],".tif",sep=""))
  aspect.county<-raster(paste(dir_rasters,"aspect_",county.txt[c],".tif",sep=""))
  # Plot and write histograms
  dem.h<-plothist(getValues(dem.out))
  slope.h<-plothist(getValues(slope.out))
  aspect.h<-plothist(getValues(aspect.out))
  
  save(dem.h,file=paste(dir_hist,"hist_elevation_",county.txt[c],sep=""))
  save(slope.h,file=paste(dir_hist,"hist_slope",county.txt[c],sep=""))
  save(aspect.h,file=paste(dir_hist,"hist_aspect",county.txt[c],sep=""))
  
}

# 2 Crop and reproject results rasters

for (y in 1:length(years)){
  for (n in 1:length(statnames)){
    print(years[y])
    print(statnames[n])
    results.file<-paste(dir_results,statnames[n],"_",year,".tif",sep="")
    print(results.file)
    r<-raster(results.file) # raster for whole area
    crs(r)<-ukgrid # old crs (OS)
    new.r<-projectRaster(r,to=template.r) # reproject 
    #plot(r)
    for (c in 1:length(county)){
      print (county.txt[c])
      ### Reproject & write county rasters
      e<-extent(county[[c]])
      county.r<-mask(crop(new.r,e), county[[c]])
      plot(county.r)
      # Write Raster -DECIDE OUTPUT DIR
      r.filename<-paste(dir_rasters,statnames[n],"_",year,"_",county.txt[c],".tif",sep="")
      print(r.filename)
      writeRaster(county.r,file=r.filename,overwrite=TRUE)
    } # for c
  } # for n
} # for y

  
# 2b Write time series files for each var and county
for (c in 1:length(county)){
  for (n in 1:length(statnames)){
    numcells<-ncell(raster(paste(dir_rasters,"elevation_",county.txt[c],".tif",sep="")))
    var.m<-matrix(data=NA,ncol=length(years),nrow=numcells)
    for (y in 1:length(years)){
      r.filename<-paste(dir_rasters,statnames[n],"_",year,"_",county.txt[c],".tif",sep="")
      r<-raster(r.filename)
      values<-getValues(r)
      var.m[,y]<-values
    } # year
    var.filename<-paste(dir_timeseries,statnames[n],"_",county.txt[c],"_timeseries.R",sep="")
    print(var.filename)
    save(var.m,file=var.filename)
  } # stat
} # county
      
  ### Save histogram (if required?)
  #h<-plothist(values)
  
  
# 3. Load vineyards shape file for whole region convert to lat lon and save file
vineyards.shp<- readOGR(dsn = "/Users/jonathanmosedale/Documents/Exeter/Data2015/vineyards", layer = "vineyardplots")
vineyards.shp<-spTransform(vineyards.shp,latlong) # CONVERT to lat lon
save(vineyards.shp,file=paste(dir_shinydata,"vineyards.RData",sep=""))



