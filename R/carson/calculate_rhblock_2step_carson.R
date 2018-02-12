##########################################################################################
# CALCULATE HOURLY 5KM RELATIVE HUMIDITY
# Input:
# Output: Daily files of hourly interpolated data
# INterpolation uses: Raster approx NA function
##########################################################################################
source("/home/ISAD/jm622/rscripts/setup_carson.R") # loads & runs setup file
args <-commandArgs(trailingOnly = TRUE)
print(args)
start.day <- as.integer(args[1])
start.month<-as.integer(args[2])
start.year<-as.integer(args[3] )
end.day<-as.integer(args[4] )
end.month<-as.integer(args[5] )
end.year<-as.integer(args[6] )

#start.day<-1; start.month<-7; start.year<-1992
#end.day<-4; end.month<-7; end.year<-1992
start.jd<-JDdmy(start.day,start.month,start.year)
end.jd<-JDdmy(end.day,end.month,end.year)

  

#######################################################################################
# FUNCTIONS TO WRITE DAILY FILES OF RELATIVE HUMIDITY
#######################################################################################
#dir_rh5km<-"~/Documents/Exeter/Data2015/RelHumidity/rh5km/"
add.zero<-function(x)
{
  y<-x
  if (y<9) y<-paste("0",x,sep="")
  y
}

rel.to.abs<-function(rh,t)
{
  e0.1<-0.6108*exp(17.27*t/(t+237.3))
  e<-e0.1*(rh/100)
  abso<-(2165*e)/(t+273.16) # grams per metre cubed
  abso
}

abs.to.rel<-function(ab,t)
{
  e<-ab*(t+273.16)/2165
  e0.1<-0.6108*exp(17.27*t/(t+237.3))
  rh<-100*(e/e0.1)
  rh
}

# Centres raster on GML 
reslice.raster<-function(r)
{
  arr<-array(0,dim=c(73,144))
  e1<-extent(1.25,181.25,-91.25,91.25)
  e2<-extent(181.25,358.75,-91.25,91.25)
  e3<-extent(-1.25,1.25,-91.25,91.25)
  arr[,1:71]<-getValues(crop(r,e2),format="matrix")
  arr[,72]<-getValues(crop(r,e3),format="matrix")
  arr[,73:144]<-getValues(crop(r,e1),format="matrix")
  r2<-raster(arr,xmn=-178.25,xmx=181.25,ymn=-91.25,ymx=91.25)
  res(r2)<-2.5
  r2
}

reproject.r<-function(a,e,dem)
{
  r<-raster(a,xmn=xmin(e),xmx=xmax(e),ymn=ymin(e),ymx=ymax(e))
  projection(r)<-"+init=epsg:4326"
  r2<-projectRaster(r,crs="+init=epsg:27700")
  r3<-resample(r2,dem)
  #r4<-crop(r3,e2)
  r3
}

rh.plot<-function(rh,hr,day,month,year) 
{
  t<-paste("RH on ",day,"-",month,"-",year," ",hr,":00",sep="")
  brk<-c(0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,160,170,180,190,200)
  col<-rev(heat.colors(21))
  plot(rh,main=t,col=col,breaks=brk)
}
###########################################################################
# Main Function: Write daily files of hourly 5km Rel Humidity 
# RH data in = 4xdaily global data
###########################################################################
calc_abs_day<-function(jd,dem,hourplot=FALSE,writefile=TRUE)
{
  e<-extent(-5.75,-0.75,48.75,53.75)  # lon/lat extent of cropped RH rasters
  #e<-extent(353.75,359.25,48.75,53.75)
  ptm<-proc.time()
  # Define output matrix to hold one day of hourly absolute humidity values at original resolution 
  abs.day<-array(0,dim=c(nrow(dem.block),ncol(dem.block),24))
  
  # Reads in Relative Humidity data from yearly file if new year
  if (jd==start.jd | (DMYjd(jd)$day==1 & DMYjd(jd)$month==1)){
    infile.rh<-paste(dir_rh,"rhum.sig995.",DMYjd(jd)$year,".nc",sep="") # reads yearly data file
    nextfile.rh<-paste(dir_rh,"rhum.sig995.",DMYjd(jd)$year+1,".nc",sep="") # reads yearly data file
    print(paste("This year's data file= ",infile.rh,sep=""))
    print(paste("Next year's data file= ",nextfile.rh,sep=""))
    rh.s<-stack(infile.rh) # all times
    nextrh.s<-stack(nextfile.rh)
    rh.s<-addLayer(rh.s,nextrh.s[[1]]) # Adds first layer of next year to currnet year
    remove(nextrh.s)
    # get times
    netRH<-nc_open(infile.rh)
    tm<-ncvar_get(netRH,"time")
    tm<-c(tm,tm[length(tm)]+6) # =adds time of first reading of next year
  }
  
  # Read in 100m daily Temperature file and resample to low res lat/lon
  filein1<-paste(dir_hrtemp,"HrTemp_", DMYjd(jd)$year, "-",sprintf("%02d",DMYjd(jd)$month,sep=""),"-", sprintf("%02d",DMYjd(jd)$day,sep=""),"_100m.r", sep="") # define file name from year,month,day,hr
  print(filein1)
  load(filein1)
  hr.temp1<-t100m.day
  filein2<-paste(dir_hrtemp,"HrTemp_", DMYjd(jd+1)$year, "-",sprintf("%02d",DMYjd(jd+1)$month,sep=""),"-", sprintf("%02d",DMYjd(jd+1)$day,sep=""),"_100m.r", sep="") # define file name from year,month,day,hr
  print(filein2)
  load(filein2)
  hr.temp2<-t100m.day
  
  for (period in 0:3)
  { print(paste("Period= ",period,sep=""))
    hr<-(period*6)
    
    # Identify rh layers and crop
    Jul.base<-JDdmy(1,1,1800)
    #Jul.actual<-JD(ISOdate(year,month,day))
    hr.val<-(jd-Jul.base)*24+hr
    sel<-which(tm==hr.val)
    rh.all0<-subset(rh.s,sel)
    rh.all6<-subset(rh.s,(sel+1))
    rh.all0<-reslice.raster(rh.all0)
    rh.all6<-reslice.raster(rh.all6)
    projection(rh.all0)<-"+init=epsg:4326"
    projection(rh.all6)<-"+init=epsg:4326"
    rh0<-crop(rh.all0,e)
    rh6<-crop(rh.all6,e) 
    
    # read in temperature rasters
    t0.100m<-raster(hr.temp1[,,hr+1],template=dem)
    if (period<3) t6.100m<-raster(hr.temp1[,,hr+7],template=dem)
    if (period==3) t6.100m<-raster(hr.temp2[,,1],template=dem) # 0hr00 for following day
    
    # reproject and calculate aggregate temperatures to match humidity data resolution
    projection(t0.100m)<-"+init=epsg:27700"
    projection(t6.100m)<-"+init=epsg:27700"
    t0.ll<-projectRaster(t0.100m,crs="+init=epsg:4326")
    t6.ll<-projectRaster(t6.100m,crs="+init=epsg:4326")
    t0.ll<-raster::extend(t0.ll,e) # required to allow aggregate values of rh cells only partly covered by temp data
    t6.ll<-raster::extend(t6.ll,e)
    tem<-raster(array(0,dim=c(2,2)),xmn=xmin(e),xmx=xmax(e),ymn=ymin(e),ymx=ymax(e))
    t0<-resample(t0.ll,rh6)
    t6<-resample(t6.ll,rh6)
    
    # Convert to absolute humidity
    abs.hum0<-rel.to.abs(getValues(rh0,format="matrix"),getValues(t0,format="matrix"))
    abs.hum6<-rel.to.abs(getValues(rh6,format="matrix"),getValues(t6,format="matrix"))
    # interpolate for missing hours
    abs.hum1<-(abs.hum0*5+abs.hum6*1)/6
    abs.hum2<-(abs.hum0*4+abs.hum6*2)/6
    abs.hum3<-(abs.hum0*3+abs.hum6*3)/6
    abs.hum4<-(abs.hum0*2+abs.hum6*4)/6
    abs.hum5<-(abs.hum0*1+abs.hum6*5)/6
    
    # Save output matrix
    abs.day[,,hr+1]<-abs.hum0
    abs.day[,,hr+2]<-abs.hum1
    abs.day[,,hr+3]<-abs.hum2
    abs.day[,,hr+4]<-abs.hum3
    abs.day[,,hr+5]<-abs.hum4
    abs.day[,,hr+6]<-abs.hum5
  } # for period
  
  # WRITE DAILY files of ABSOLUTE humidity in lon/lat at resolution of origina RH data
  if (writefile==TRUE){
    fileout<-paste(dir_rh5km,"AbsHum_",DMYjd(jd)$year,"_",DMYjd(jd)$month,"_",DMYjd(jd)$day,".R",sep="")
    print(fileout)
    save(abs.day,file=fileout)
  } # end if
} 
    


calc_RH_day<-function(jd,t100m.day,dem,dem.block,hourplot=FALSE)
{
  # tref.day = matrix of dem.block tref data for all 24 hrs
  filein<-paste(dir_rh5km,"AbsHum_",DMYjd(jd)$year,"_",DMYjd(jd)$month,"_",DMYjd(jd)$day,".R",sep="")
  print(filein)
  load(filein) # abs.day
  e<-extent(-5.75,-0.75,48.75,53.75)  # lon/lat extent of cropped RH rasters
  abshum.day<-array(0,dim=c(nrow(dem.block),ncol(dem.block),24))
  rh.day<-array(0,dim=c(nrow(dem.block),ncol(dem.block),24))
  for (hr in 0:23){
    # Resample to dem.block
    abshum.block<-reproject.r(abs.day[,,hr+1],e,dem.block)
    # Crop t data to dem.block
    tref.block<-crop(raster(t100m.day[,,hr+1],template=dem),dem.block)
    # Convert to relative humidity
    rh.day[,,hr+1]<-abs.to.rel(getValues(abshum.block),getValues(tref.block))
  }
  print(proc.time()-ptm)
  if (hourplot==TRUE){
    for (t in 1:24) {
      rh.plot(raster(rh.day[,,t],template=dem.block),t-1,DMYjd(jd)$day,DMYjd(jd)$month,DMYjd(jd)$year) 
    } }
} # end function


#######################################################################################
# FUNCTION CALLS
for (jd in start.jd:end.jd){
  calc_abs_day(jd,dem)
}

calc_RH_day(jd,t100m.day,dem,dem.block,hourplot=FALSE)
  
