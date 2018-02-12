##########################################################################################
# CALCULATE HOURLY 5KM LW RADIATION 
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

## Downscale daily TEMPERATURE data to hourly  Prog: t5km_to_hrmatrix ###
#start.day<-1; start.month<-1; start.year<-1983
#end.day<-31; end.month<-12; end.year<-2013
start.jd<-JDdmy(start.day,start.month,start.year)
end.jd<-JDdmy(end.day,end.month,end.year)

#######################################################################################
# FUNCTIONS TO CALCULATE LWR from CAL, RH and T
#######################################################################################
# CALCULATE LW RADIATION from CAL , Temp, RH
# Units of output = MJ/m2
# Inputs: rh cropped to block, tref requires cropping, CAL - requires downscaling and cropping
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

calc_lwr_block<-function(jd,cal.block,rhref.block,tref.block,plotlwr=FALSE,writefile=FALSE)
{
  # Define output file - one day of hourly 100m data
    #proc.time()->ptm
    year<-DMYjd(jd)$year
    month<-DMYjd(jd)$month
    day<-DMYjd(jd)$day
    #print(paste("Date: ",day,"/",month,"/",year,sep=""))
    compareRaster(cal.block,rhref.block,tref.block)
    
    # Calculate longwave radiation
    lwr.m<-lwr(getValues(tref.block,format="matrix"),
                       getValues(rhref.block,format="matrix"),
                       getValues(cal.block,format="matrix"))
    lwr.block<-raster(lwr.m,template=tref.block)
    
    if (plotlwr==TRUE) {
      par(mfrow=c(2,2))
      dayhr<-paste(day,"/",month,"/",year," ",hr+1,"h00",sep="")
      plot(cal.block,main=paste("Effective cloud albedo ",dayhr,sep=""))
      plot(rhref.block,main=paste("Relative humdidity ",dayhr,sep=""))
      plot(tref.block,main=paste("Temperature ",dayhr,sep=""))
      plot(lwr.block,main=paste("Long-wave radiation ",dayhr,sep=""))
    }
    #print(proc.time()-ptm)
    if (writefile==TRUE){
    file.out<-paste(dir_lwr,"lwr_",year,"_",month,"_",day,"_100m.R",sep="")
    print(paste("File out: ",file.out,sep=""))
    save(lwr.day,file=file.out)
    print(proc.time()-ptm)
    } # if writefile
} # end function


#######################################################################################
# FUNCTION CALL
#######################################################################################
# Calculate LONG WAVE RADITAION at 5km hourly res ffrom CAL, RH T  Prog: longwav_grids
#rhref.block<-crop(raster(rh.day[,,hr+1],template=dem),dem.block)
#tref.block<-crop(raster(t100m.day[,,hr+1],template=dem),dem.block)
#calc_lwr_block(jd,cal.block,rhref.block,tref.block,plotlwr=TRUE,writefile=FALSE)

