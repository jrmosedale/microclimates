# Calculate key temperature parameters for each 100m cell
# Process by year and block
# Output by block x risk: xyz file of seasonal values with z as year 
## ASSUME start end end correspond to start and end of single YEAR or GROWING SEASON

# Carson job variables
args <-commandArgs(trailingOnly = TRUE)
print(args)
year <- as.integer(args[1])

# For testing only: 
#  root<-"~/Documents/Exeter/Data2015/"; in.root<-"~/Documents/Exeter/Data2015/" 

# Source data and output data
source("/home/ISAD/jm622/rscripts/setup_v4_carson.R") 
dir_results<-paste(root,"Outputs/",sep="")

# Set cells
cells<-c(1:935)
numcells<-length(cells)
print(paste("Analysing ",numcells," grid cells",sep=""))

# Define growing season
gsstart.day<-1
gsstart.month<-4
gsend.day<-31
gsend.month<-10

# calculate number of days in yr 
numdays<-daydoy(year,12,31)
print(paste("Number of Days in ",year," is ", numdays,sep=""))

# Calculate days of year for season start/end
gsstart.doy<-daydoy(year,gsstart.month,gsstart.day)
gsend.doy<-daydoy(year,gsend.month,gsend.day)

par(mfrow=c(2,3))

####################################################################################
# Create lists to merge blocks into single raster

gdd10_gs<-vector("list",numcells)
gdd5_gs<-vector("list",numcells)
tmean_gs<-vector("list",numcells)
tmin_year<-vector("list",numcells)
tmax_year<-vector("list",numcells)
t20_gsdays<-vector("list",numcells)
t25_gsdays<-vector("list",numcells)
t30_gsdays<-vector("list",numcells)

lastspfr_doy<-vector("list",numcells)
firstautfr_doy<-vector("list",numcells)
frostfree_days<-vector("list",numcells)

fl_tmean<-vector("list",numcells)
fl_numday<-vector("list",numcells)

i<-1 # index for blocks

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
  print(paste("Analysing data for Year= ",year, " and Cell= ",ukcpcell,sep=""))
  
  # load this year's temperature data at 100m and hourly resolution
  infile<-paste(dir_finalt,"block-",sprintf("%03d",ukcpcell,sep=""),"-",year,".R",sep="")
  print(infile)
  load(infile) # = tmodel.r
  if (dim(tmodel.r)[3]!=numdays*24) warning("Number of days in tmodel.r NOT equal to number of hours in year")

  ####################################################################################
  # Calculate Seasonal Risks 
  ####################################################################################
  start.gs<-(gsstart.doy*24)-23 # layer coresponding to 0h00 on 1st day
  end.gs<-(gsend.doy*24) # later for 23h00 for last day
  
  # 1. Calculate GS tmean and GDD using hourly data
  tmean.gs<-apply(tmodel.r[,,start.gs:end.gs], c(1,2), function(x) mean(x))
  
  hgdd10.gs<-apply(tmodel.r[,,start.gs:end.gs], c(1,2), function(x) sum(ifelse(x>10,(x-10)/24,0)))
  hgdd5.gs<-apply(tmodel.r[,,start.gs:end.gs], c(1,2), function(x) sum(ifelse(x>5,(x-5)/24,0)))
  #plot(raster(hgdd10.gs,template=dem.block),main=paste("GDD(hr) Tbase=10C 1/4-30/10 ",ukcpcell," ",year,sep=""))
  
  # 2a. Yearly min/max  temperature 
  tmin.year<-apply(tmodel.r[,,1:(numdays*24)], c(1,2), function(x) min(x))
  tmax.year<-apply(tmodel.r[,,1:(numdays*24)], c(1,2), function(x) max(x))
  
  # 2b Number of GS days with max t> t
  start.hr<-start.gs
  t.day<-array(rep(0,ncell(dem.block)),dim=c(nrow(dem.block),ncol(dem.block),(1+end.gs-start.gs)/24) )
  n<-1
  t20.numday<-array(rep(0,ncell(dem.block)),dim=c(nrow(dem.block),ncol(dem.block)) )
  t25.numday<-array(rep(0,ncell(dem.block)),dim=c(nrow(dem.block),ncol(dem.block)) )
  t30.numday<-array(rep(0,ncell(dem.block)),dim=c(nrow(dem.block),ncol(dem.block)) )
  
  while (start.hr<=end.gs){ # calculate daily mean temperatures
    #print(flhr)
    t.day[,,n]<-apply(tmodel.r[,,start.hr:(start.hr+23)], c(1,2), function(x) max(x))
    start.hr<-start.hr+24
    n<-n+1
  } # end while
  
  t20.numday<-apply(t.day, c(1,2), function(x) {length(which(x>=20))} )
  t25.numday<-apply(t.day, c(1,2), function(x) {length(which(x>=25))} )
  t30.numday<-apply(t.day, c(1,2), function(x) {length(which(x>=30))} )
  #plot(raster(t20.numday,template=dem.block),main=paste("Days Tmean>20 " ,ukcpcell," ",year,sep=""))
  #plot(raster(t25.numday,template=dem.block),main=paste("Days Tmean>25 " ,ukcpcell," ",year,sep=""))
  #plot(raster(t30.numday,template=dem.block),main=paste("Days Tmean>30 " ,ukcpcell," ",year,sep=""))
  
  # 3. Calculate last spring and first fall frost - returns Days of Year and not limited to growing season
  # Define max temperature which is considered a frost event = frostrisk.t in degrees C 
  t.frost<-1
   
  # NB Assumes no frost between 31 May & 1 September. REQUIRE LATER MASKING WHEN CONVERTED TO RASTER TO ENSURE SEA=NA
  end.spfrostrisk<-daydoy(year,5,31)*24 # last hr to be checked for sp frost even (end May)
  lastspfr.doy<-apply(tmodel.r[,,1:end.spfrostrisk], c(1,2), function(x) {ifelse(length(which(x<=t.frost))>0,tail(which(x<=t.frost),1)+1,1)})
  lastspfr.doy<-ceiling(lastspfr.doy/24)
  #plot(raster(lastspfr.doy,template=dem.block),main=paste("Last Spring frost ",ukcpcell," ",year,sep=""))
  
  start.autfrostrisk<-(daydoy(year,9,1)*24)-23
  firstautfr.doy<-apply(tmodel.r[,,start.autfrostrisk:(numdays*24)], c(1,2), function(x) {ifelse(length(which(x<=t.frost))>0, head(which(x<=t.frost),1)+start.autfrostrisk, numdays*24 )})
  firstautfr.doy<-ceiling(firstautfr.doy/24)
  #plot(raster(firstautfr.doy,template=dem.block),main=paste("First Autumn Frost ",ukcpcell," ",year,sep=""))
  
  # Calculate frost free period
  frostfree.days<-firstautfr.doy - lastspfr.doy
  #plot(raster(frostfree.days,template=dem.block),main=paste("Frost free period of year ",DMYjd(jd)$year,sep=""))
  
  # 4. Calculate flowering risks
  
  # Calculate days of flowering
  fl.start<-daydoy(year,6,22) # doy
  fl.end<-daydoy(year,7,5) # doy
  # Convert to hours
  fl.start<-(fl.start*24)-23
  fl.end<-fl.end*24 
  
  # define critical mean t
  t.fl<-15
  
  # Mean flowering temperature - all hours
  fl.tmean<- apply(tmodel.r[,,fl.start:fl.end], c(1,2), function(x) mean(x))
  #plot(raster(fl.tmean,template=dem.block),main=paste("Mean flowering T  ",ukcpcell," ",year,sep=""))

  # Number of flowering days where mean T>15C
  flhr<-fl.start
  fl.numday<-array(rep(0,ncell(dem.block)),dim=c(nrow(dem.block),ncol(dem.block)) )
  while (flhr<=fl.end){
    #print(flhr)
    good.day<-apply(tmodel.r[,,flhr:(flhr+23)], c(1,2), function(x) {ifelse(mean(x)>t.fl,1,0)} )
    fl.numday<-fl.numday+good.day
    #plot(raster(apply(tmodel.r[,,flhr:(flhr+23)], c(1,2), function(x) {mean(x)} ),template=dem.block),main="Good day?")
    flhr<-flhr+24
  } # end while
  #plot(raster(fl.numday,template=dem.block),main=paste("Flowering days Tmean>",t.fl," ",ukcpcell," ",year,sep=""))
  
  # Add block raster to lists
  gdd10_gs[[i]]<-raster(as.matrix(hgdd10.gs),template=dem.block)
  gdd5_gs[[i]]<-raster(as.matrix(hgdd5.gs),template=dem.block)
  tmean_gs[[i]]<-raster(as.matrix(tmean.gs),template=dem.block)
  tmin_year[[i]]<-raster(as.matrix(tmin.year),template=dem.block)
  tmax_year[[i]]<-raster(as.matrix(tmax.year),template=dem.block)
  
  t20_gsdays[[i]]<-raster(as.matrix(t20.numday),template=dem.block)
  t25_gsdays[[i]]<-raster(as.matrix(t25.numday),template=dem.block)
  t30_gsdays[[i]]<-raster(as.matrix(t30.numday),template=dem.block)
  
  lastspfr_doy[[i]]<-raster(as.matrix(lastspfr.doy),template=dem.block)
  firstautfr_doy[[i]]<-raster(as.matrix(firstautfr.doy),template=dem.block)
  frostfree_days[[i]]<-raster(as.matrix(frostfree.days),template=dem.block)
  
  fl_tmean[[i]]<-raster(as.matrix(fl.tmean),template=dem.block)
  fl_numday[[i]]<-raster(as.matrix(fl.numday),template=dem.block)
  
  i<-i +1
} # end cell

# Merge blocks into single raster for each stat, plot and write raster

# List of statistics
statistics<-list(gdd10_gs,gdd5_gs,tmean_gs,tmin_year,tmax_year,
                 t20_gsdays,t25_gsdays,t30_gsdays, 
                 lastspfr_doy, firstautfr_doy, frostfree_days,
                 fl_tmean, fl_numday)
statnames<-c("gdd10_gs","gdd5_gs","tmean_gs","tmin_year","tmax_year",
                 "t20_gsdays","t25_gsdays","t30_gsdays", 
                "lastspfr_doy", "firstautfr_doy", "frostfree_days",
                "fl_tmean", "fl_numday")

for(n in 1:length(statistics)){
  
  map.r<-do.call(merge,statistics[[n]] )
  map.r<-mask(map.r,crop(dem,map.r))
  fileout<-paste(dir_results,statnames[n],"_",year,".tif",sep="")
  print(fileout)
  plot(map.r,main=fileout)
  writeRaster(map.r,file=fileout,format="GTiff",overwrite=TRUE)
}


