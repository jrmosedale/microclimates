# Prepare files required for shinyapp
# Reproject to esprg and save in appropriate directories 


library(rgdal)
library(rgeos)
library(raster)

root<-"~/Documents/Exeter/Data2015/"; in.root<-"~/Documents/Exeter/Data2015/"

# Input dirs
dir_cornwall_in<-paste(root,"proxyt100/riskmaps/cornwall/",sep="")
dir_devon_in<-paste(root,"proxyt100/riskmaps/devon/",sep="")
dir_dorset_in<-paste(root,"proxyt100/riskmaps//dorset/",sep="")
dir_somerset_in<-paste(root,"proxyt100/riskmaps/somerset/",sep="")

dir_dem<-paste(in.root,"DEM/",sep="")
dir_terrain<-paste(in.root,"Terrain/",sep="")

# Output dirs
dir_shinydata<-paste(root,"shinydata/",sep="")
dir_swfiles<-paste(root,"shinydata/swfiles/",sep="")
dir_cornwall<-paste(root,"shinydata/cornwall/",sep="")
dir_devon<-paste(root,"shinydata/devon/",sep="")
dir_dorset<-paste(root,"shinydata/dorset/",sep="")
dir_somerset<-paste(root,"shinydata/somerset/",sep="")


latlong = "+init=epsg:4326"
epsg<-"+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext +no_defs"
#epsg<-"+init=epsg:3857" # mercator

#####################################################################
# Read shape file for county boundary 
#####################################################################
#counties.shp<- readOGR(dsn = "/Users/jonathanmosedale/Documents/Exeter/Data2015/OSdata/bdline_essh_gb/Data/GB", layer = "county_region")
#divisions.shp<- readOGR(dsn = "/Users/jonathanmosedale/Documents/Exeter/Data2015/OSdata/bdline_essh_gb/Data/GB", layer = "district_borough_unitary_region")
#west.shp<- readOGR(dsn = "/Users/jonathanmosedale/Documents/Exeter/Data2015/OSdata/bdline_essh_gb/Data/GB", layer = "westminster_const_region")
hcounties.shp<- readOGR(dsn = "/Users/jonathanmosedale/Documents/Exeter/Data2015/OSdata/bdline_essh_gb/Data/Supplementary_Historical", layer = "Boundary-line-historic-counties_region")
hwater.shp<- readOGR(dsn = "/Users/jonathanmosedale/Documents/Exeter/Data2015/OSdata/bdline_essh_gb/Data/GB", layer = "high_water_polyline")
ogrInfo(dsn = "/Users/jonathanmosedale/Documents/Exeter/Data2015/OSdata/bdline_essh_gb/Data/GB", layer = "county_region")
#print(divisions.shp$POLYGON_ID)

#plot(divisions.shp,main="district_borough_unitary_region")
# extract SW coast outlines
swcoast.shp<-crop(hwater.shp,dem)

# extract county polygons - use historic counties data 
sel <- which(hcounties.shp$Name == "Cornwall")
cornwall<-hcounties.shp[sel,]
#plot(cornwall)
sel <- which(hcounties.shp$Name == "Devon")
devon<-hcounties.shp[sel,]
#plot(devon)
sel <- which(hcounties.shp$Name == "Dorset")
dorset<-hcounties.shp[sel,]
#plot(dorset)
sel <- which(hcounties.shp$Name == "Somerset")
somerset<-hcounties.shp[sel,]
#plot(somerset)

# PLOT - overlay counties 
plot(swcoast.shp)
plot(devon,add=TRUE)
plot(somerset,add=TRUE)
plot(cornwall,add=TRUE)
plot(dorset,add=TRUE)

# Load vineyards shape file and crop to cornwall
vineyards.shp<- readOGR(dsn = "/Users/jonathanmosedale/Documents/Exeter/Data2015/vineyards", layer = "vineyardplots")
#vineyards.os<-spTransform(vineyards.shp,crs(hcounties.shp)) # project to OS - same as county shapes
#vineyards<-crop(vineyards.os,cornwall)
#vineyards.epsg<-spTransform(vineyards,epsg)
vineyards.shp<-spTransform(vineyards.shp,latlong) # CONVERT to lat lon
save(vineyards.shp,file="/Users/jonathanmosedale/Documents/Exeter/Data2015/shinydata/vineyards.RData")



######################################

######################################
new.crs<-"+init=epsg:3857"

# Link terrain map data to files and prepare for plotting
dem<-raster(paste(dir_dem,"dem.tif",sep=""))
slope<-crop(raster(paste(dir_terrain,"slope.tif",sep="")),dem)
aspect<-crop(raster(paste(dir_terrain,"aspect.tif",sep="")),dem)

dem.ll<-projectRaster(dem,crs=new.crs,filename=paste(dir_swfiles,"demepsg.tif",sep=""),overwrite=TRUE)
slope.ll<-projectRaster(slope,crs=new.crs,filename=paste(dir_swfiles,"slopeepsg.tif",sep=""),overwrite=TRUE)
aspect.ll<-projectRaster(aspect,crs=new.crs,filename=paste(dir_swfiles,"aspectepsg.tif",sep=""),overwrite=TRUE)

#map.s<-stack(dem.epsg,slope.epsg,aspect.epsg)

######################################

######################################
# Extract county rasters - quite slow  and result same size as whole rasters
# for every raster in input directory assumes all files .tif
# Input files: crs=
######################################
#new.crs<-'+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext +no_defs'
new.crs<-"+init=epsg:3857"
years<-c("2009","2010","2011")
county.txt<-"cornwall"
county<-cornwall 
e<-extent(county)

# Crop and reproject terrain data for each county
#plot(mask(crop(dem,county), county))
dem.out<-projectRaster(mask(crop(dem,county), county),crs=new.crs)
slope.out<-projectRaster(mask(crop(slope,county), county),crs=new.crs)
aspect.out<-projectRaster(mask(crop(aspect,county), county),crs=new.crs)

writeRaster(dem.out,file=paste(dir_shinydata,county.txt,"/dem.tif",sep=""),overwrite=TRUE)
writeRaster(slope.out,file=paste(dir_shinydata,county.txt,"/slope.tif",sep=""),overwrite=TRUE)
writeRaster(aspect.out,file=paste(dir_shinydata,county.txt,"/aspect.tif",sep=""),overwrite=TRUE)


for (y in 1:length(years)){
  print(years[y])
  dir_input<-paste(dir_cornwall_in,years[y],"/",sep="")
  filelist<-list.files(dir_input)
  print(filelist)
  
  dir_output<-paste(dir_shinydata,county.txt,"/",years[y],"/",sep="")
  print(dir_output)
  
  for (n in 1:length(filelist)){
    print(filelist[n])
    r<-raster(paste(dir_input,filelist[n],sep=""))
    filename<-substr(filelist[n],1,nchar(filelist[n])-4)
    #plot(r,main=filename)
  
    # mask and write county rasters
    #new.r<-mask(crop(r,e), cornwall)
    new.r<-projectRaster(r,crs=new.crs)
    #projection(new.r)<-CRS(new.crs)
    #plot(new.r)
    
    # write raster
    writeRaster(new.r,file=paste(dir_output,filename,".tif",sep=""),overwrite=TRUE)
    
  } # for n
} # for y
  

######################################
# Save matrices for timeseries lookup
# Save histograms of cell values
library(ggplot2)
# Function that will return file name
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
roundUpNice <- function(x, nice=c(1,2,4,5,6,8,10)) {
  if(length(x) != 1) stop("'x' must be of length 1")
  10^floor(log10(x)) * nice[[which(x <= 10^floor(log10(x)) * nice)[[1]]]]
}

plothist<-function(lab,values){
  sel<-which(!is.na(values))
  values<-values[sel]
  mnx<-min(values)
  mxx<-max(values)
  f <- hist(values, maxpixels=length(values),breaks=50) # calculate from all cells of raster
  #abline(v=100, col="blue",lwd=3)
  dat <- data.frame(counts= ((f$counts/numcells)*100),breaks = f$mids)
  ggplot(dat, aes(x = breaks, y = counts, fill =counts)) + ## Note the new aes fill here
    geom_bar(stat = "identity",alpha = 0.8,fill="blue")+
    xlab(lab)+ ylab("%")+
    scale_x_continuous(breaks = seq(0,mxx,(roundUpNice(mxx/5))),
                       labels = seq(0,mxx,(roundUpNice(mxx/5))) )
}

# One matrix per variable row=cell number, col=year
dir_shinydata<-paste(root,"shinydata/",sep="")

# Uses getmap function in app
# for each variable
county<-"cornwall"
numcells<-ncell(raster(paste(dir_cornwall,"dem.tif",sep="")))
vars<-c("gdd10","gdd5","maxt","mint","meant","num20","num25","num30", "spfr","autfr","frfree")
labs<-c("Degree days","Degree days","degrees C","degrees C","degrees C","Number of days","Number of days","Number of days",
          "Day of Year (1:366)","Day of Year (1:366)","Number of days")

#y<-1; v<-1
for (v in 1: length(vars)){
  var<-vars[v]
  var.m<-matrix(data=NA,ncol=length(years),nrow=numcells)

for (y in 1:length(years)){
  r<-getmap(county,var,years[y])
  values<-getValues(r)
  plot(r)
  #plothist(labels[v],values)
  # save values to matrix
  var.m[,y]<-values
} # for y
  # write file of matrix - cell x year
  print(paste("Mean=",mean(var.m,na.rm=TRUE)))
  print(length(which(!is.na(var.m))))
  filename<-paste(dir_shinydata,county,"/",var,"-timeseries.R",sep="")
  print(filename)
  save(var.m,file=filename)
} # for v

###################################################################
### Prepare histogrms of whole rasters
library(ggplot2)
value<-100

# 1 Elevation histogram
# Number of non-sea cells 
numcells<-length(which(!is.na(getValues(dem.epsg))))
print(numcells)
f <- hist(dem.epsg, maxpixels=ncell(dem.epsg),breaks=50) # calculate from all cells of raster
#abline(v=100, col="blue",lwd=3)
dat <- data.frame(counts= ((f$counts/numcells)*100),breaks = f$mids)

dem.hist<-ggplot(dat, aes(x = breaks, y = counts, fill =counts)) + ## Note the new aes fill here
  geom_bar(stat = "identity",alpha = 0.8,fill="blue")+
  xlab("Elevation")+ ylab("%")+
  scale_x_continuous(breaks = seq(0,500,50),
                     labels = seq(0,500,50))

#+  scale_fill_gradient(low="blue", high="red")
#  + geom_vline(xintercept = value,color = "red", size=1)
#ggsave(paste(dir_shinydata,"demhist")
# Add vertical line marking individual (cell) value

# 2 Slope histogram
f <- hist(slope.epsg, maxpixels=ncell(slope.epsg),breaks=50) # calculate from all cells of raster
#abline(v=100, col="blue",lwd=3)
dat <- data.frame(counts= ((f$counts/numcells)*100),breaks = f$mids)

slope.hist<-ggplot(dat, aes(x = breaks, y = counts, fill =counts)) + ## Note the new aes fill here
  geom_bar(stat = "identity",alpha = 0.8,fill="blue")+
  xlab("Slope")+ ylab("%")+
  scale_x_continuous(breaks = seq(0,30,10),
                     labels = seq(0,30,10))

# 3 Aspect histogram
f <- hist(aspect.epsg, maxpixels=ncell(aspect.epsg),breaks=50,xlim=c(0,360)) # calculate from all cells of raster
#abline(v=100, col="blue",lwd=3)
dat <- data.frame(counts= ((f$counts/numcells)*100),breaks = f$mids)

aspect.hist<-ggplot(dat, aes(x = breaks, y = counts, fill =counts)) + ## Note the new aes fill here
  geom_bar(stat = "identity",alpha = 0.8,fill="blue")+
  xlab("Aspect")+ ylab("%")+
  scale_x_continuous(breaks = seq(0,360,50),
                     labels = seq(0,360,50))





# Cut Outs
  
  # Set extent of counties
  e.cornwall<-extent(cornwall)
  e.devon<-extent(devon)
  e.dorset<-extent(dorset)
  e.somerset<-extent(somerset)
  # ... or set to extent of whole of SouthWest
  # e.cornwall<-extent(r); e.devon<-extent(r); e.dorset<-extent(r); e.somerset<-extent(r)
  
for (n in 1:length(filelist)){
  r<-raster(paste(dir_input,filelist[n],sep=""))
  filename<-substr(filelist[n],1,nchar(filelist[n])-4)
  plot(r,main=filename)
  
  # mask and write county rasters
  devon.r <- mask(crop(r,e.devon), devon)
  cornwall.r<-mask(crop(r,e.cornwall), cornwall)
  dorset.r<-mask(crop(r,e.dorset), dorset)
  somerset.r<-mask(crop(r,e.somerset), somerset)
  
  # Reproject to epsg ???
  devon.r<-projectRaster(devon.r,crs=epsg,res=c(100,100))
  cornwall.r<-projectRaster(cornwall.r,crs=epsg,res=c(100,100))
  dorset.r<-projectRaster(dorset.r,crs=epsg,res=c(100,100))
  somerset.r<-projectRaster(somerset.r,crs=epsg,res=c(100,100))
  
  writeRaster(devon.r,file=paste(dir_devon,filename,".tif",sep=""),overwrite=TRUE)
  writeRaster(cornwall.r,file=paste(dir_cornwall,filename,".tif",sep=""),overwrite=TRUE)
  writeRaster(dorset.r,file=paste(dir_dorset,filename,".tif",sep=""),overwrite=TRUE)
  writeRaster(somerset.r,file=paste(dir_somerset,filename,".tif",sep=""),overwrite=TRUE)
}

######################################




###################################################################
# Density plot?
f <- density(dem.epsg,plot=TRUE)
dat <- data.frame(counts= f$counts,breaks = f$mids)
ggplot(dat, aes(x = breaks, y = counts, fill =counts)) + ## Note the new aes fill here
  geom_bar(stat = "identity",alpha = 0.8)+
  xlab("Elevation")+ ylab("Frequency")+
  scale_x_continuous(breaks = seq(-1,1,0.25),
                     labels = seq(-1,1,0.25))+
  scale_fill_gradient(low="blue", high="red") 