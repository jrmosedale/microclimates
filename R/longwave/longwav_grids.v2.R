# CALCULATE LW RADIATION from CAL , Temp, RH
# Units of output = MJ/m2
# Historic data: run rel.hum.v2, cal.cal.hrly & t5km_to_hrmatrix
# Requires jd functions

# FUNCTION - calculate long wave radiation from Temp, RH & CAL
lwr<-function(Temp,RH,CAL)
{
  e0<-0.6108*exp(17.27*Temp/(Temp+237.3)) # saturated vapour pressure
  ea<-e0*(RH/100) # actual vapour pressure
  moisture.absorb<-0.34-0.14*sqrt(ea)
  cloud.absorb<-1.35*(1-CAL)-0.35
  rnl<-2.043*10^-10*(Temp+273.16)^4*moisture.absorb*cloud.absorb
  rnl
}

####################################

write_lwr_dayfies<-function(start.jd,end.jd,grid5km.r,plotlwr=FALSE)
{
  # Define output file - one day of hourly 5km data
  lwr.day<-array(0,dim=c(nrow(grid5km.r),ncol(grid5km.r),24))
  for (t in start.jd:end.jd)
  {
    year<-DMYjd(t)$year
    month<-DMYjd(t)$month
    day<-DMYjd(t)$day
    
    # Reads in daily files - all assumed to be 5km OSGB gridcells with identical extents of grid5km.r
    rh.filein<-paste(dir_rh5km,"RH_5km_",year,"_",month,"_",day,".R",sep="")
    print(rh.filein)
    load(rh.filein) # rh.day
    
    t.filein<-paste(dir_hrtemp,"HrTemp_", year, "-",sprintf("%02d",month,sep=""),"-", sprintf("%02d",day,sep=""),".r", sep="")
    print(t.filein)
    load(t.filein) # t5km.day
    
    cal.filein<-paste(dir_calimp,"CALimp_5km_",year,"_",month,"_",day,".R",sep="")
    print(cal.filein)
    load(cal.filein) # CALimp.day
  
     for (hr in 1:24)
      {
         # rel hum
         m.rh<-rh.day[,,hr]
         r.rh<-raster(m.rh,template=grid5km.r)
         # temperature
         m.temp<-t5km.day[,,hr]
         r.temp<-raster(m.temp,template=grid5km.r)
         # CAL 
         m.cal<-CALimp.day[,,hr]
         r.cal<-raster(m.cal,template=grid5km.r)
         
         # Calculate longwave radiation
         lwr.day[,,hr]<-lwr(getValues(r.temp,format="matrix"),
                             getValues(r.rh,format="matrix"),
                             getValues(r.cal,format="matrix"))
         
         if (plotlwr==TRUE) {
             par(mfrow=c(2,2))
             dayhr<-paste(day,"/",month,"/",year," ",hr-1,"h00",sep="")
             plot(r.cal,main=paste("Effective cloud albedo ",dayhr,sep=""))
             plot(r.rh,main=paste("Relative humdidity ",dayhr,sep=""))
             plot(r.temp,main=paste("Temperature ",dayhr,sep=""))
             plot(raster(lwr.day[,,hr],template=r.rh),main=paste("Long-wave radiation ",dayhr,sep=""))
         }
       } # for hr
    # WRITE DAILY lwr files
    file.out<-paste(dir_lwr,"lwr_",year,"_",month,"_",day,".R",sep="")
    print(file.out)
    save(lwr.day,file=file.out)
         
    } # for day
} # end function































