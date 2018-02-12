#################################################################################
# 2. Means across seasons
# load seasonal files into single stack of 30 years of data 1983-2014
# create empty stack
tmax.s<-stack()
tmin.s<-stack()
tmean.s<-stack()
gdd10.s<-stack()
gdd5.s<-stack()
days20.s<-stack()
days25.s<-stack()
days30.s<-stack()
spfrost.s<-stack()

# Load 5km data for entire time period
for (jd in 1983:2014){
  # Read day temperature data
  tmax.s<-stack(tmax.s,raster(file=paste(dir_results,"maxt_5km_",year ,".grd" ,sep="") )
  tmin.s<-stack(tmin.s,raster(file=paste(dir_results,"mint_5km_",year ,".grd" ,sep="") )
  tmean.s<-stack(tmean.s,raster(file=paste(dir_results,"meant_5km_",year ,".grd" ,sep="") )
  days20.s<-stack(days20.s,raster(file=paste(dir_results,"days20_5km_",year ,".grd" ,sep="") )
  days25.s<-stack(days25.s,raster(file=paste(dir_results,"days25_5km_",year ,".grd" ,sep="") )
  days30.s<-stack(days30.s,raster(file=paste(dir_results,"days30_5km_",year ,".grd" ,sep="") )
  gdd10.s<-stack(gdd10.s,raster(file=paste(dir_results,"gdd10_5km_",year ,".grd" ,sep="") )
  gdd5.s<-stack(gdd5.s,raster(file=paste(dir_results,"gdd5_5km_",year ,".grd" ,sep="") )                        
                                                                                                                          spfrost.s <-stack(spfrost.s,raster(file=paste(dir_results,"maxt_5km_",year ,".grd" ,sep="") ) 
}

# Calculate summaries of whole seasons 
# For specific 5km site
site<-  # raster cell number

# Mean of each variable for time periods
yr1<-1983; yr2<-
  
  
# Time series plots
  
# Calculate difference from norm (gdd-mean gdd )
  
