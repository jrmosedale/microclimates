library(ncdf4)
library(mgcv) # require package
library(raster)
# # # # # # # # # # # # # # # #
# Imputes missing values of CAL
# # # # # # # # # # # # # # # #

# FUNCTION - imputation
# input: dfi = dataset with missing values (effective Cloud albedo with missing night time values)
# output: dataset with missing values imputed (effective Cloud albedo with no missing night time values)
imputation<-function(dfi)
{
  # Ensure values lie between zero and one
  sel<-which(dfi<0.001); dfi[sel]<-0.001
  sel<-which(dfi>0.999); dfi[sel]<-0.999
  dfi<-log(dfi/(1-dfi))
  x<-c(1:length(dfi))
  # use spline to interpolate
  sp<-spline(x,dfi,n=length(x))
  dfi<-1/(1+exp(-1*sp$y))
  dfi
}

# **********************

# FUNCTION - cal.5km.impute - writes 5km matrix with imputed nightime values
# Writes files of imputed CAL values from files missing nightime values
# Ensure CAL input files already cropped to correct extent

cal.5km.impute<-function(dir_cal,dir_calimp,start.jd,end.jd,grid5km,plotcal=FALSE){
  files<-list.files(dir_cal) 
  
  # Define initial cropping in lat/lon and resulting number of cols/rows to be stored
  e<-extent(c(-7,-2,49,52)) # extent for initial cropping 
  e.5km<-extent(grid5km.r)
  numrows<-nrow(grid5km.r)
  numcols<-ncol(grid5km.r)
  max.i<-(1+end.jd-start.jd)*24 # = total number of hrs to be stored
  CAL.vals<-array(NA,c(numrows,numcols,max.i))
  
  # load and store CAL files 
  i<-1
  for (t in start.jd:end.jd){
    for (hr in 0:23){ 
      year<-DMYjd(t)$year; month<-DMYjd(t)$month ; day<-DMYjd(t)$day 
      # check if file exists and print warning and record missing if not
      infile.cal<-paste(dir_cal,"CALhm",year,sprintf("%02d",month,sep=""),sprintf("%02d",day,sep=""),sprintf("%02d",hr,sep=""),"00002UD1000101UD.nc",sep="")
      print(infile.cal)
      #if (infile.cal not in files) print("File not found") 
      
      # Resample CAL data to 5km grid
      cal.r<-raster(infile.cal)
      projection(cal.r)<-"+init=epsg:4326"
      calosgb.r<-projectRaster(cal.r,crs="+init=epsg:27700")
      #plot(calosgb.r)
      cal5km.r<-raster::resample(crop(calosgb.r,e.5km),grid5km.r)
      #plot(cal5km.r)
      #print(dim(m))
      #print (i)
      m<-getValues(cal5km.r,format="matrix")
      CAL.vals[,,i]<-m
      i<-i+1
      # else{warning(paste("Files not found: ",infile.cal,sep=""))}
      }# end hr loop
    }# end day loop
  
  # Impute missing values for each grid cell in turn
  CAL.imp<-array(NA,c(numrows,numcols,max.i))
  for (r in 1:numrows){
    for (c in 1:numcols){
      v<-CAL.vals[r,c,]
      CAL.imp[r,c,]<-imputation(v)
      print(paste("row: ",r," col: ",c,sep=""))
    } # rows
  } # cols
  
  
  
  # Write DAILY files of imputed CAL values at 5km resolution 
  for (t in start.jd:end.jd){
    CALimp.day<-array(NA,c(numrows,numcols,24))
    first<-(t-start.jd)*24+1
    last<-first + 23
    print(paste("Day: ",t," first=",first," last=",last,sep=""))
    CALimp.day[,,1:24]<-CAL.imp[,,first:last]
    
    # plot each hour of each day ?
    if (plotcal==TRUE){
      par(mfrow=c(4,3))
      for (i in 1:24){
        calsw.r<-raster(CALimp.day[,,i],template=grid5km.r)
        titletext<-paste(DMYjd(t)$day,"/",DMYjd(t)$month,"/",DMYjd(t)$year," ",i,"hr00",sep="")
        plot(calsw.r,main=titletext)
      } }
    
    fileout<-paste(dir_calimp,"CALimp_5km_",DMYjd(t)$year,"_",DMYjd(t)$month,"_",DMYjd(t)$day,".R",sep="")
    print (fileout)
    save(CALimp.day,file=fileout) 
  }# end for write files
  
} # end function


#######################################################
# TESTING OF FUNCTIONS

start.jd<-JDdmy(1,7,1992) 
end.jd<-JDdmy(10,7,1992)

dir_cal<-"C:/Data2015/CMSAF-CAL/extract/"
dir_calimp<-"C:/Data2015/CMSAF-CAL/imputed/"


# 1st function
cal.impute(dir_cal,dir_calimp,start.jd,end.jd)


projection(cal.r)<-"+init=espg:4326"
calosgb.r<-projectRaster(cal.r,crs="+init=epsg:27700")
plot(calosgb.r)
cal5km.r<-resample(crop(calosgb.r,e.5km),grid5km.r)
plot(cal5km.r)

# 2nd function
cal.5km.impute(dir_cal,dir_calimp,start.jd,end.jd,grid5km)


load("~/Documents/Exeter/Data2015/CMSAF-CAL/imputed/CALimp_5km_1992_7.R")


#######################################################


# ******** CUT OUTS **********
# FUNCTION - cal.impute - writes lat/lon matrix at res 0.05 with imputed nightime values
# Writes files of imputed CAL values from files missing nightime values
# Ensure CAL input files already cropped to correct extent
cal.impute<-function(dir_cal,dir_calimp,start.jd,end.jd){
  files<-list.files(dir_cal) 
  
  # Define initial cropping in lat/lon and resulting number of cols/rows to be stored
  #e<-extent(c(-6.725,-2.675,49.775,51.375)) # for testing
  e<-extent(c(-7,-2,49,52)) # gives 60 rows and 100 cols
  numcols<-100; numrows<-60
  max.i<-(1+end.jd-start.jd)*24 # = total number of hrs to be stored
  CAL.vals<-array(NA,c(numrows,numcols,max.i))
  
  # load and store CAL files 
  # t<-start.jd; hr<-0
  i<-1
  for (t in start.jd:end.jd){
    for (hr in 0:23){ 
      year<-DMYjd(t)$year; month<-DMYjd(t)$month ; day<-DMYjd(t)$day 
      # check if file exists and print warning and record missing if not
      infile.cal<-paste(dir_cal,"CALhm",year,sprintf("%02d",month,sep=""),sprintf("%02d",day,sep=""),sprintf("%02d",hr,sep=""),"00002UD1000101UD.nc",sep="")
      print(infile.cal)
      #if (infile.cal in files) 
      cal.r<-raster(infile.cal)
      cal.r<-crop(cal.r,e)
      m<-getValues(cal.r,format="matrix")
      CAL.vals[,,i]<-m
      i<-i+1
      # else{warning(paste("Files not found: ",infile.cal,sep=""))}
      #print(paste(day,"/",month,"/",year," ",hr,"hrs",sep=""))
      
    }# end hr loop
  }# end day loop
  
  # Impute missing values for each grid cell in turn
  CAL.imp<-array(NA,c(numrows,numcols,max.i))
  for (r in 1:numrows){
    for (c in 1:numcols){
      v<-CAL.vals[r,c,]
      CAL.imp[r,c,]<-imputation(v)
      print(paste("row: ",r," col: ",c,sep=""))
      
    } # for r
  } # for c
  fileout<-paste(dir_calimp,"CALimp_",year,"_",month,".R",sep="")
  print (fileout)
  save(CAL.imp,file=fileout)
} # end function

# MOnthly files - CUT OUT
# Assumes start.jd = 1st of a month
for (jd in start.jd:end.jd){
  lastofmonth<-days.in.month(DMYjd(jd)$day,DMYjd(jd)$month,DMYjd(jd)$year)
  if (DMYjd(jd)$day==lastofmonth) {
    firstofmonth<-(jd-start.jd)-lastofmonth+2
    first<-(1+jd-start.jd-lastofmonth)*24
    last<-(1+jd-start.jd)*24
    CAL.imp.month<-CAL.imp[r,c,first:last] # for whole month of data
    fileout<-paste(dir_calimp,"CALimp_5km_",DMYjd(jd)$year,"_",DMYjd(jd)$month,".R",sep="")
    print (fileout)
    save(CAL.imp.month,file=fileout) }
  
  if(jd==end.jd & DMYjd(jd)$day!=lastofmonth){
    firstofmonth<-(jd-start.jd)-DMYjd(jd)$day+2
    end<-jd-start.jd+1
    CAL.imp.month<-CAL.imp[r,c,firstofmonth:end] # for whole month of data
    fileout<-paste(dir_calimp,"CALimp_5km_",DMYjd(jd)$year,"_",DMYjd(jd)$month,".R",sep="")
    print (fileout)
    #save(CAL.imp.month,file=fileout) 
  }
  