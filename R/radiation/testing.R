par(mfrow=c(2,2))
jd<-start.jd
for (i in 0:23){
  print(i)
  #if (i==0) n<-0 else n<-i/24
  #print(sza.angle(jd+n,-2,50))
  # Calculate jd including hour fraction
  h<--12+i; #print(h)
  h.jd<-h/24;#print(h.jd)
  jd.h<-jd+h.jd;#print(jd.h)
  
  # Using insol package
  sv<-(sunvector(jd.h,lat,long,0));print(paste("Sunvector=",sv,sep="")) # z angle >0 in daylight
  print(paste("Sunpos=",sunpos(sv),sep="")) # daylight solar zenith<=90; prints azimuth and zenith
  print(paste("sza.angle=",sza.angle(jd.h,lat,long,tz=0)* (pi/180),sep=""))
  
  print(paste("Solalt=",solalt(i,lat,long,jd),sep=""))#saltitude
  print(paste("Solalt=",solalt(i,lat,long,jd)*(pi/180),sep=""))
  zen<-pi/2 -(solalt(i,50,-5,jd)*(pi/180)); print(paste("Zen=",zen,sep=""))
  print(paste("Sazimuth=",solazi(i,lat,long,jd),sep=""))
  azi<-solazi(i,lat,long,jd)* (pi/180); print(paste("Azi=",azi,sep=""))
  #horangle<-horizonangle(dtm,solazi(i,lat,long,jd))
  horangle<-horizonangle(dtm,sunpos(sv)[1])
  #plot(raster(horangle,template=dem.buffer),main=i)
  shadowmask <- array(1,dim(dtm))
  shadowmask[horizonangle(dtm,sazimuth)>tan(solalt(i,lat,long,jd)*(pi/180))] <- 0
  print(paste("Tan(alt)= ",tan(solalt(i,lat,long,jd)*(pi/180)),sep=""))
  #plot(raster(shadowmask,template=dem.buffer))
  shade.r<-doshade(dem.buffer,sv)
  #plot(mask(shade.r,dem.buffer),main="shade.r")
 
  asp<-360* (pi/180); sl<-15* (pi/180)
  index <- cos(zen) * cos(sl) + sin(zen) * sin(sl) * cos(azi - asp)
  print(paste("Index= ",index,sep=""))
 
  ins<-insolation(sunpos(sv)[2], jd.h, 20,100,90,288,0.01,0.2)
}
sazimuth<-solazi(i,lat,long,jd)
dtm=array(0,dim=c(1,1))
si<-solarindex(slope=m.slope,aspect=m.aspect,localtime=hr,
               Lat=lat,Long=long,Julian=jd,dtm=m.dem)
solarindex <- function(slope,aspect,localtime,Lat,Long,Julian,dtm=array(0,dim=c(1,1)),res=100,merid=0,dst=0,shadow=TRUE)
{
  saltitude<-solalt(localtime,Lat,Long,Julian,merid,dst)
  alt <- saltitude * (pi/180)
  zen <- pi/2 - alt
  sazimuth<-solazi(localtime,Lat,Long,Julian,merid,dst)
  azi <- sazimuth * (pi/180)
  sl <- slope * (pi/180)
  asp <- aspect * (pi/180)
  shadowmask <- array(1,dim(dtm))
  horangle<-horizonangle(dtm,sazimuth)
  if(shadow) {
    shadowmask[horizonangle(dtm,sazimuth)>tan(alt)] <- 0
  }
  index <- array(0,dim(dtm))
  index <- cos(zen) * cos(sl) + sin(zen) * sin(sl) * cos(azi - asp)
  index[index<0] <- 0
  index <- index * shadowmask
  index
}

h.r<-raster(horizonangle(dtm,solazi(i,lat,long,jd)),template=dem.buffer)
plot(h.r,main="h.r")
h.nxt<- calc(h.r,fun=function(x){ifelse(x>0.05,x<-0,x<-1)} )
plot(h.nxt,main="h.nxt")


# Calculate SZA and create raster of SZA values
sza<-sza.angle(jdays,lt,ln)/360  # in decial degrees!!!
#csza<-cos(sza)

sza.m<-matrix(sza,nrow=nrow(dni.r),byrow=TRUE) # convert array to matrix
#csza.m<-matrix(csza,nrow=nrow(dni.r),byrow=TRUE) # convert array to matrix

sza.r<-raster(sza.m, template=dni.r) # convert matrix to raster
#csza.r<-raster(csza.m, template=dni.r) # convert matrix to raster

# Calculate SID as dni * cosine(sza)
sid.r<-overlay(x = dni.r,y = sza.r,fun = function(x,y) {
  z<-x*cos(y) 
  return(z)
} ) 

# Calculate indirect from SIS and SID
ind.r<-overlay(x = sis.r,y = sid.r,fun = function(x,y) {
  z<-x-y 
  return(z)
} ) 

dif.r<-overlay(x = dni.r,y = sid.r,fun = function(x,y) {
  z<-x-y 
  return(z)
} ) 
dif2.r<-overlay(x = sis.r,y = dni.r,fun = function(x,y) {
  z<-x-y 
  return(z)
} )  



# Insol package
sunvector(jd,lat,long,0)



# MODIFICATIONS
solarindex <- function(slope,aspect,localtime,Lat,Long,Julian,dtm=array(0,dim=c(1,1)),res=100,merid=0,dst=0,shadow=TRUE)
{
  saltitude<-array(solalt(localtime,Lat,Long,Julian,merid,dst),dim(dtm))
  alt <- saltitude * (pi/180)
  zen <- pi/2 - alt
  sazimuth<-array(solazi(localtime,Lat,Long,Julian,merid,dst),dim(dtm))
  azi <- sazimuth * (pi/180)
  sl <- slope * (pi/180)
  asp <- aspect * (pi/180)
  shadowmask <- array(1,dim(dtm))
  horangle<-horizonangle(dtm,sazimuth)
  #plot(raster(horangle,template=dem.buffer))
  if(shadow) {
    shadowmask[horizonangle(dtm,sazimuth)>tan(alt)] <- 0
  }
  index <- array(0,dim(dtm))
  index <- cos(zen) * cos(sl) + sin(zen) * sin(sl) * cos(azi - asp)
  index[index<0] <- 0
  index <- index * shadowmask
  plot(raster(index,template=dem.buffer))
  index
}

par(mfrow=c(2,2))
for (hr in 6:20){
  # Construct time label for plots 
  timelabel<-paste(day,"/",month,"/",year," ",hr,"h00",sep="")
  print(timelabel)
  
  print("Radiation")
  # Extract hourly data from day stack
  dnr.r<-raster(dnr.24h,layer=hr+1); #plot(dnr.r)
  sis.r<-raster(sis.24h,layer=hr+1); #plot(sis.r)
  rad.results<-radiation_downscale_stack(day,month,year,hr,sis.r,dnr.r,dem.buffer,dem.block,slope.buffer,aspect.buffer,print.results=FALSE)
  direct.block<-rad.results[[1]]*W.to.MJhr
  diffuse.block<-rad.results[[2]]*W.to.MJhr
  total.block<-rad.results[[3]]*W.to.MJhr
  #plot(total.block,main=paste("Total MJ ",timelabel,sep=""))
  ### ALTITUDE EFFECT ON TEMPERATURE ###
  elev.tdif<- (elevdif.block/100)*0.66  # altitude effect=1.98 per 300m - see Wikipedia!
  
  
  anom.block<-(params$estimates[1]+elev.tdif+   # difference in temperature due to altitude
                 2*total.block+      # net short wave radiation
                 params$estimates[5]*albedo+        # albedo
                 params$estimates[14]*total.block*albedo)
  plot(total.block,main=paste("Total Rad ",timelabel,sep=""))
  plot(direct.block,main=paste("Direct ",sep=""))
  plot(diffuse.block,main=paste("Diffuse ",sep=""))               
  plot(anom.block,main=paste("ANOMALY ",sep=""))               
  
  #plot(params$estimates[3]*total.block,main="Temperature effect")
}
hr<-11

day<-15
filein1<-paste(dir_dniday,"DNIhm",year,sprintf("%02d",month,sep=""),sprintf("%02d",day,sep=""),sep="")
filein2<-paste(dir_sisday,"SIShm",year,sprintf("%02d",month,sep=""),sprintf("%02d",day,sep=""),sep="")
dnr.24h<-stack(filein1)
sis.24h<-stack(filein2)
for (hr in 5:20){
plot(crop(raster(dnr.24h,layer=hr),dem.buffer),main="dnr")
plot(crop(raster(sis.24h,layer=hr),dem.buffer),main="sis")
}

for (hr in 6:20){
  dnr.r<-raster(dnr.24h,layer=hr+1); #plot(dnr.r)
  sis.r<-raster(sis.24h,layer=hr+1); #plot(sis.r)
  rad.results<-radiation_downscale_stack(day,month,year,hr,sis.r,dnr.r,dem.buffer,dem.block,slope.buffer,aspect.buffer,print.results=TRUE)
  
}



# Calculate Rdir and Rdif from SI/SV to give Rtotal
sza<-sza.angle(jd.h,lat,long,tz=0)* (pi/180)  # in decial degrees!!!

sza.m<-matrix(sza,nrow=nrow(dni.r),byrow=TRUE) # convert array to matrix
#csza.m<-matrix(csza,nrow=nrow(dni.r),byrow=TRUE) # convert array to matrix

sza.r<-raster(sza.m, template=dni.r) # convert matrix to raster
#csza.r<-raster(csza.m, template=dni.r) # convert matrix to raster



# TEST AND TRIALS ON DIFFERENT APPROACHES
# Trying out eq from Bennie et al
a<-ifelse(wstr>1,0.013,0.022)

T<-Tsh+Ta+Tb*
  
  Tsh<-tref.block+a*(lwr.block+diffuse.block)
Ta<-a*direct.block*(1-albedo)
Tb<-a*direct.block*(sin(sza)*sin(slope)*(1-albedo)
                    
                    #simplified
                    t<-tref.block+a*total.block
                    plot(t)
                    
                    
                    
                    # Test Rdir equation 1
                    dnr.block<-crop(dnr.buffer,dem.block)
                    plot(dnr.block)
                    
                    slope.block<-crop(slope.buffer,dem.block)
                    aspect.block<-crop(aspect.buffer,dem.block)* (pi/180)
                    slope.cos<-cos(slope.block* (pi/180))
                    slope.sin<-sin(slope.block* (pi/180))
                    
                    h<--12+hr; #print(h)
                    h.jd<-h/24;#print(h.jd)
                    jd.h<-jd+h.jd;#print(jd.h)
                    
                    sv<-(sunvector(jd.h,lat,long,0));print(paste("Sunvector=",sv,sep="")) # z angle >0 in daylight
                    s.azi<-sunpos(sv)[1]* (pi/180)
                    s.zen<-sunpos(sv)[2]* (pi/180)
                    print(s.azi);print(s.zen)
                    
                    # Convert to rasters for calculation
                    block.m<-getValues(dem.block,format="matrix")
                    
                    s.zen.r<-raster(array(s.zen,dim(block.m)),template=dem.block)
                    s.azi.r<-raster(array(s.azi,dim(block.m)),template=dem.block)
                    
                    
                    Rdir<-dnr.block*( cos(s.zen.r)*slope.cos + sin(s.zen.r)*slope.sin*cos(s.azi.r-aspect.block) )
                    
                    # OR
                    
                    si.block<-crop(raster(si,template=dem.buffer),dem.block)
                    plot(tref.block*(si.block/si.flat))
                    
                    
                    
                    # Calculate mean Rref on flat surface across whole 5km grid cell
                    #   =SIS ?? (ie Rdir(=SID)+Rdif)
                    #Calculate SI/SIflat - correct for shadow areas etc - set to 0 for these areas
                    
                    #Calculate SV
                    
                    # Calculate SID as dni * cosine(sza)
                    sid.r<-overlay(x = dni.r,y = sza.r,fun = function(x,y) {
                      z<-x*cos(y) 
                      return(z)
                    } ) 
                    # Rtotal/Rref * T