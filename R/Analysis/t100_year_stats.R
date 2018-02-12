# Calculate key temperature parameters for each 100m cell
# Process by year and block
# Output by block x risk: xyz file of seasonal values with z as year 

# Carson job variables
args <-commandArgs(trailingOnly = TRUE)
print(args)
start.day <- as.integer(args[1])
start.month<-as.integer(args[2])
start.year<-as.integer(args[3] )
end.day<-as.integer(args[4] )
end.month<-as.integer(args[5] )
end.year<-as.integer(args[6] )
ukcpcell<-as.integer(args[7])

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
plot(dem.block)


#source("/home/ISAD/jm622/rscripts/setup_carson.R") # loads & runs setup file (jd functions, dem etc)

## ASSUME start end end correspond to start and end of single YEAR or GROWING SEASON
#start.day<-1; start.month<-7; start.year<-1992
#end.day<-2; end.month<-7; end.year<-1992
# dir_temp<-"C:/Data2015/Temp5km/extract/"
# dir_results<-"C:/Results2015/year_stats_5km/"
# dir_allyr<-"C:/Results2015/allyear_stats_5km"
# Uses 5km grid mask: grid5km.r
# dir_grids<-"C:/Data2015/Templates/"
# dir_finalt<-"C:/Data2015/Temp5km/hourly/"

dir_finalt<-"~/Documents/Exeter/Data2015/Temp100m/"
start.year<-1991
end.year<-1999
ukcpcell<-900
#year<-start.year

start.jd<-JDdmy(1,1,start.year)
end.jd<-JDdmy(31,12,end.year)
print(start.jd)
print(end.jd)
#jd<-start.jd

####################################################################################
# for (ukcpcell in cells[1]:cells[length(cells)]){}

# for (year in start.year:end.year){}
print(paste("Analysing data for Year= ",year,sep=""))
# load this year's temperature data at 100m and hourly resolution
infile<-paste(dir_finalt,"block-",ukcpcell,"-",year,".R",sep="")
print(infile)
load(infile) # = tmodel.r
# calculate number of days in yr
days.in.yr<-dim(tmodel.r)[3]/24
print(days.in.yr)

####################################################################################
# Calculate daily stats for whole year- ie daily min/max etc
####################################################################################
tmin.year<-stack()
tmax.year<-stack()
tmean.year<-stack()

for (doy in 1:days.in.yr){
  # Extract 24h of temperature data
  start<-1+(doy-1)*24
  end<-doy*24
  #print (paste(start,"  ",end))
  t.24h<-tmodel.r[,,start:end]
  # Calculate daily statistics
  tmin.24h<-apply(t.24h, c(1,2), min)
  tmax.24h<-apply(t.24h, c(1,2), max)
  tmean.24h<-apply(t.24h, c(1,2), mean)
  
  # Add layers to stack of summary variables STACK or ARRAY?
  tmin.year<-addLayer(tmin.year,raster(as.matrix(tmin.24h),template=dem.block))
  tmax.year<-addLayer(tmax.year,raster(as.matrix(tmax.24h),template=dem.block))
  tmean.year<-addLayer(tmean.year,raster(as.matrix(tmean.24h),template=dem.block))
  
}

####################################################################################
# Calculate Seasonal Rasters 
####################################################################################
# Define Growing Season (or part of year of interest)
start.gs<-1; end.gs<-nlayers(tmean.year)

# 1. Calculate last spring and first fall frost - not limited to growing season
# Calc last spring frost <= 1 C 
# Assumes no frost between end of May & early September - so reduces vector to 1-150 doy
spfrdata.s<-subset(tmin.year,1:150)
start.v<-rep(start.gs,(nlayers(spfrdata.s)))
spfrost.r<-calc(spfrdata.s,fun=function(x){ifelse(length(which(x<=2))>0,tail(which(x<=2),1)+start.v,1)}) # extract layer of last frost day
spfrost.r<-mask(spfrost.r,dem.block,maskvalue=NA)
plot(spfrost.r,main=paste("Last spring frost day ",DMYjd(jd)$year,sep=""))

# Calculate first autumn frost (after early sept)
autfrdata.s<-subset(tmin.year,240:nlayers(tmin.year))
start.v<-rep(240,(nlayers(autfrdata.s)))
autfrost.r<-calc(autfrdata.s,fun=function(x){ifelse(length(which(x<=2))>0,head(which(x<=2),1)+start.v,nlayers(tmin.year))}) # extract layer of last frost day
autfrost.r<-mask(autfrost.r,dem.block,maskvalue=NA)
plot(autfrost.r,main=paste("First autumn frost day ",DMYjd(jd)$year,sep=""))

# Calculate frost free period
frostfree.r<-overlay(spfrost.r,autfrost.r,fun=function(x,y){return(y-x)})
plot(frostfree.r,main=paste("Frost free period of year ",DMYjd(jd)$year,sep=""))


# 2. Calculate growing season stats - temperature extremes, gdd etc
#start.gs<-90; end.gs=305
# correct start and end dates to reflect zone of interest
#end.jd<-start.jd+end.gs-1
#start.jd<-start.jd+start.gs-1
# OR SET GS to first and last frosts??

# Calc gdd 
Tbase<-10;  tbase.v<-rep(Tbase,(end.gs-start.gs+1))
gdd10.r<-calc(subset(tmean.year,start.gs:end.gs),fun=function(x){sum(x-tbase.v)})

gdd10.r<-mask(gdd10.r,dem.block,maskvalue=NA)
plot(gdd10.r,main=paste("GDD10 ",year," Tbase= ",Tbase,sep=""))

Tbase<-5;  tbase.v<-rep(Tbase,(end.gs-start.gs+1))
gdd5.r<-calc(subset(tmean.year,start.gs:end.gs),fun=function(x){sum(x-tbase.v)})
plot(gdd5.r,main=paste("GDD5 ",year," Tbase= ",Tbase,sep=""))

# Calc mean T
meangst.r<-calc(subset(tmean.year,start.gs:end.gs),fun=function(x){mean(x)})
plot(meangst.r,main=paste("MeanT ",year,sep=""))

# Calc max T
maxgst.r<-calc(subset(tmax.year,start.gs:end.gs),fun=function(x){max(x)})
plot(maxgst.r,main=paste("Max T ",year,sep=""))

# Calc #days where max temp>  20C, 25C, 30C from April-Oct
days20.r<-calc(subset(tmax.s,start.gs:end.gs),fun=function(x){ifelse(length(which(x>=20))>0,length(which(x>=20)),0)} )
days20.r<-mask(days20.r,grid5km.r,maskvalue=NA)

days25.r<-calc(subset(tmax.s,start.gs:end.gs),fun=function(x){ifelse(length(which(x>=25))>0,length(which(x>=25)),0)} )
days25.r<-mask(days25.r,grid5km.r,maskvalue=NA)

days30.r<-calc(subset(tmax.s,start.gs:end.gs),fun=function(x){ifelse(length(which(x>=30))>0,length(which(x>=30)),0)} )
days30.r<-mask(days30.r,grid5km.r,maskvalue=NA)

plot(days20.r,main=paste("Days=>20 ",DMYjd(jd)$year,sep=""))
plot(days25.r,main=paste("Days=>25 ",DMYjd(jd)$year,sep=""))
plot(days30.r,main=paste("Days=>30 ",DMYjd(jd)$year,sep=""))

# Calc min T in growing season
mingst.r<-calc(subset(tmin.s,start.gs:end.gs),fun=function(x){min(x)})
plot(mingst.r,main=paste("Min T ",DMYjd(jd)$year,sep=""))



####################################################################################
# Output raster files for year - seasonal characteristics by block
####################################################################################
gdd10.fileout<-paste(dir_results,"block-",ukcpcell,"-gdd10-",year,".tif" ,sep="")
gdd5.fileout<-paste(dir_results,"block-",ukcpcell,"-gdd5-",year,".tif" ,sep="")
meant.fileout<-paste(dir_results,"block-",ukcpcell,"-meant-",year ,".tif" ,sep="")
maxt.fileout<-paste(dir_results,"block-",ukcpcell,"-maxt-",year ,".tif" ,sep="")
days20.fileout<-paste(dir_results,"block-",ukcpcell,"-days20-",year ,".tif" ,sep="")
days25.fileout<-paste(dir_results,"block-",ukcpcell,"-days25-",year ,".tif" ,sep="")
days30.fileout<-paste(dir_results,"block-",ukcpcell,"-days30-",year ,".tif" ,sep="")
mint.fileout<-paste(dir_results,"block-",ukcpcell,"-mint-",year ,".tif" ,sep="")
spfrost.fileout<-paste(dir_results,"block-",ukcpcell,"-spfrost-",year ,".tif" ,sep="")
autfrost.fileout<-paste(dir_results,"block-",ukcpcell,"-autfrost-",year ,".tif" ,sep="")
frostfree.fileout<-paste(dir_results,"block-",ukcpcell,"-frostfree-",year ,".tif" ,sep="")

writeRaster(gdd10.r,file=gdd10.fileout,format="GTiff")
writeRaster(gdd5.r,file=gdd5.fileout,format="GTiff")
writeRaster(maxt.r,file=maxt.fileout,format="GTiff")
writeRaster(days20.r,file=days20.fileout,format="GTiff")
writeRaster(days25.r,file=days25.fileout,format="GTiff")
writeRaster(days30.r,file=days30.fileout,format="GTiff")
writeRaster(mint.r,file=mint.fileout,format="GTiff")
writeRaster(spfrost.r,file=spfrost.fileout,format="GTiff")
writeRaster(autfrost.r,file=autfrost.fileout,format="GTiff")
writeRaster(frostfree.r,file=frostfree.fileout,format="GTiff")




# CUT OUTS
# TESTING ON INTERPOLATED 5km DATA
#day.file<-paste(dir_hrtemp,"HrTemp_", DMYjd(jd)$year[1], "-",sprintf("%02d",DMYjd(jd)$month[1],sep=""),"-", sprintf("%02d",DMYjd(jd)$day[1],sep=""),"_100m.r", sep="") 
#print(day.file)
#load( file=file.out) # loads: t100m.day

# Calculate daily statistics
day.tmin<-apply(t100m.day, c(1,2), min)
day.tmax<-apply(t100m.day, c(1,2), max)
day.tmean<-apply(t100m.day, c(1,2), mean)
day.GDH10<-
  day.GDH0<-
  day.GDH15<-
  
  day.GDH20<-
  day.GDH25<-
  day.GDH30<-
  
  day.FDH0<-
  day.FDH2<-
  
  day.H25<-
  day.H30<-
  
  # Daily summary stats= max, min, mean T, GDhr (2 ver - 0C, 10C, 15C??), ~hrs>maxT, #hrs<minT, frost degree hrs = FDhr, 
  # ...
  
  # Add layers to stack of summary variables
  tmin.year<-addLayer(tmin.year,tmin.day)


# Calculate by cell - vector analysis after extracting from stack??
for (cell in 1:ncell(t100.24[[1]])){    
}

# Or for every day calculate summary stats and create yearly stack with 365 layers - apply functions to bricks
