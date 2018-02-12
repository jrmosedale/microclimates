# Calculate key temperature parameters for each 100m cell
# Process by year and block
# Output by block x risk: xyz file of seasonal values with z as year 

# Carson job variables - none

print("Calling set up...")
source("/home/ISAD/jm622/rscripts/setup_carson.R") # loads & runs setup file

#dir_finalt<-"~/Documents/Exeter/Data2015/Temp100m/"
dir_gdd<-paste(root,"Outputs/gdd/",sep="")
dir_frost<-paste(root,"Outputs/frost/",sep="")
dir_flowering<-paste(root,"Outputs/flowering/",sep="")
dir_extremes<-paste(root,"Outputs/extremes/",sep="")

cells<-c(744,745,769,770)
#cells<-c(718,719,720,744,745,746,769,770,771,898,899,900,911,912,913,921,922,923)
years<-c(1991,1992,1993,1994,1995,1996,1997,1999,2000,2001,2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,2012)

print(cells)
print(years)
#par(mfrow=c(2,3))

for (ukcpcell in cells){
  # Get dem.block for cell
  print(paste("UKCP cell= ", ukcpcell,sep=""))
  print("Get Cell Coordinates")
  gridmask.r<-land5km.r # land5km.r defined in setup
  vals<-getValues(gridmask.r)
  xy<-xyFromCell(gridmask.r,1:ncell(gridmask.r))
  sel<-which(vals==1)
  landcells<-xy[sel,1:2] # = coordinates for middle of each ukcp09 cell
  print(dim(landcells))
  x<-landcells[ukcpcell,1]
  y<-landcells[ukcpcell,2]
  e.block<-extent(x-2500,x+2500,y-2500,y+2500)
  dem.block<-crop(demuk,e.block) 
  #plot(dem.block,main=paste("DEM ",ukcpcell,sep=""))
  
  for (year in years){
    print(paste("Analysing data for Year= ",year, " and for cell= ",ukcpcell,sep=""))
    # load this year's temperature data at 100m and hourly resolution
    infile<-paste(dir_finalt,"block-",ukcpcell,"-",year,".R",sep="")
    print(infile)
    load(infile) # = tmodel.r
    # calculate number of days in yr
    numdays<-dim(tmodel.r)[3]/24
    #print(numdays)
    
    ####################################################################################
    # Calculate daily stats for whole year- ie daily min/max etc
    ####################################################################################
    tmin.year<-stack()
    tmax.year<-stack()
    tmean.year<-stack()
    
    for (doy in 1:numdays){
      # Extract 24h of temperature data
      start<-1+(doy-1)*24
      end<-doy*24
      #print (paste(start,"  ",end))
      t.24h<-tmodel.r[,,start:end]
      # Calculate daily statistics
      tmin.24h<-apply(t.24h, c(1,2), min)
      tmax.24h<-apply(t.24h, c(1,2), max)
      tmean.24h<-apply(t.24h, c(1,2), mean)
      # hgdd10.24h<-apply(t.24h, c(1,2), function(x) sum(ifelse(x>10,(x-10)/24,0))) - daily - could be summed later
      # Add layers to stack of summary variables STACK or ARRAY?
      tmin.year<-addLayer(tmin.year,raster(as.matrix(tmin.24h),template=dem.block))
      tmax.year<-addLayer(tmax.year,raster(as.matrix(tmax.24h),template=dem.block))
      tmean.year<-addLayer(tmean.year,raster(as.matrix(tmean.24h),template=dem.block))
    }
    
    # Calculaate seasonal gdd using hourly data
    start.gdd<-1
    end.gdd<-244*24
    hgdd10.year<-apply(tmodel.r[,,start.gdd:end.gdd], c(1,2), function(x) sum(ifelse(x>10,(x-10)/24,0)))
    #hgdd5.year<-apply(tmodel.r[,,start.gdd:end.gdd], c(1,2), function(x) sum(ifelse(x>5,(x-5)/24,0)))
    #plot(raster(hgdd10.year,template=dem.block),main=paste("GDD(hr) Tbase=10C 1/4-30/10 ",ukcpcell," ",year,sep=""))
    
    # Write GDDhr file 
    outfile<-paste(dir_gdd,"GDD10hr-aproct-",ukcpcell,"-",year,".tif",sep="")
    print(outfile)
    writeRaster(raster(hgdd10.year,template=dem.block),file=outfile,overwrite=TRUE)
    
    ####################################################################################
    # Calculate Seasonal Rasters 
    ####################################################################################
    # Define Growing Season (or part of year of interest)
    start.gs<-90 # for 1st of March - ignoring leap years
    end.gs<-304 # end Oct or 334 for end of Nov
    
    # Calculate min temperature from 1 April to end May
    spfrdata.s<-subset(tmin.year,1:62) # to doy 150 if whole year data
    spring.tmin<-min(spfrdata.s)
    #plot(spring.tmin,main=paste("Min T 04-05 ",ukcpcell," ",year,sep=""))
    
    # 1. Calculate last spring and first fall frost - not limited to growing season
    # Calc last spring frost <= 1 C 
    # Assumes no frost between start of April and end of May & early September - end Nov 
    spfrdata.s<-subset(tmin.year,1:62) # to doy 150 if whole year data
    start.v<-rep(start.gs,(nlayers(spfrdata.s)))
    spfrost.r<-calc(spfrdata.s,fun=function(x){ifelse(length(which(x<=2))>0,tail(which(x<=2),1)+start.v,1)}) # extract layer of last frost day
    spfrost.r<-mask(spfrost.r,dem.block,maskvalue=NA)
    #plot(spfrost.r,main=paste("Last frost day ",ukcpcell," ",year,sep=""))
    
    # Calculate first autumn frost (after early sept) - needs tweaking
    #autfrdata.s<-subset(tmin.year,240:nlayers(tmin.year))
    #start.v<-rep(240,(nlayers(autfrdata.s)))
    #autfrost.r<-calc(autfrdata.s,fun=function(x){ifelse(length(which(x<=2))>0,head(which(x<=2),1)+start.v,nlayers(tmin.year))}) # extract layer of last frost day
    #autfrost.r<-mask(autfrost.r,dem.block,maskvalue=NA)
    #plot(autfrost.r,main=paste("First autumn frost day ",DMYjd(jd)$year,sep=""))
    
    # Calculate frost free period
    #frostfree.r<-overlay(spfrost.r,autfrost.r,fun=function(x,y){return(y-x)})
    #plot(frostfree.r,main=paste("Frost free period of year ",DMYjd(jd)$year,sep=""))
    
    # Calculate mean flowering temperature
    fl.start<-175-start.gs # doy
    fl.end<-188-start.gs # doy
    #  fl.start<-fl.start-start.gs; fl.end<-fl.end-start.gs   # adjust if not full year data
    fl.tmean<-calc(subset(tmean.year,fl.start:fl.end),fun=function(x){mean(x)})
    #plot(fl.tmean,main=paste("Mean flowering T  ",ukcpcell," ",year,sep=""))
    fl.numday<-calc(subset(tmean.year,fl.start:fl.end),fun=function(x){length(x[x>15])})
    fl.numday<-mask(fl.numday,dem.block,maskvalue=NA)
    #plot(fl.numday,main=paste("Flowering days Tmean>15C  ",ukcpcell," ",year,sep=""))
    
    # Write season rasters
    outfile<-paste(dir_frost,"SpringTmin-aprmay-",ukcpcell,"-",year,".tif",sep="")
    print(outfile)
    writeRaster(spring.tmin,file=outfile,overwrite=TRUE)
    
    outfile<-paste(dir_frost,"LastFrost-aprmay-",ukcpcell,"-",year,".tif",sep="")
    print(outfile)
    writeRaster(spfrost.r,file=outfile,overwrite=TRUE)
    
    outfile<-paste(dir_flowering,"Flowering-tmean-",ukcpcell,"-",year,".tif",sep="")
    print(outfile)
    writeRaster(fl.tmean,file=outfile,overwrite=TRUE)
    
    outfile<-paste(dir_flowering,"Flowering-numday-",ukcpcell,"-",year,".tif",sep="")
    print(outfile)
    writeRaster(fl.numday,file=outfile,overwrite=TRUE)
    
  } # end years
} # end ukcpcell


