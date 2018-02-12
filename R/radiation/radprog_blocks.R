####################################################################################
# FUNCTIONS USED
####################################################################################

# Function for printing multiple maps of 1+ layers in same stack. 
# Argumennts: stack, vector of layer names, hour (for label) 
plot.stack<-function(stk,lyrs,hour) {
  # Define common colour scheme and scale for plots
  #par(mfrow=c(2,2))
  brk<-c(-50,0,100,200,300,400,500,600,700,800,900,1000)
  col<-rev(rainbow(11,start=1/6,end=4/6))
  for (map in 1:length(lyrs)){
    text<-paste(lyrs[map]," at ", hour, ":00",sep="")
    plot(x=stk,lyrs[map],col=col,breaks=brk, main="")
    title(main=text)
  } # for loop
} # function

# USES FUNCTION JDdoy(doy,year) to give JD

# Converts a raster to a matrix for use with solar index functions NOT USED - USE getValues(r,format="matrix")
use.raster<-function(r)
{
     xr<-dim(r)[1]
     xc<-dim(r)[2]
     m<-array(getValues(r),dim=c(xc,xr))
     m<-t(m) # transpose
     m
}

# OS GB grid ref to Lat and Long
OSGBtolatlong<-function(x,y)
{
  pt = data.frame(x,y)
  coordinates(pt)=~x+y
  proj4string(pt)=CRS("+init=epsg:27700")
  latlong<-spTransform(pt,CRS("+init=epsg:4326"))
  ll<-as.data.frame(latlong)
  ll
}

# Needed for solar index function
solartime <- function(localtime,Long,Julian,merid=0,dst=0)
{
	Bn <- 2 * 3.141 * (Julian - 81) / 364
	eot <- 9.87 * sin(2 * Bn) - 7.53 * cos(Bn) - 1.5 * sin(Bn)
	solartime <- localtime + (4 / 60) * (2 * 3.141 * (merid - Long) / 360) + (eot / 60) - dst
	solartime
}

solalt <- function(localtime,Lat,Long,Julian,merid=0,dst=0)
{
  stime<-solartime(localtime,Long,Julian,merid,dst)
	tt <- 0.261799 * (stime - 12)
	declin <- (pi * 23.5 / 180) * cos(2 * pi * ((Julian - 171) / 365.25))
	Sinh = sin(declin) * sin(Lat * pi / 180) + cos(declin) * cos(Lat * 3.141 / 180) * cos(tt)
	solalt = (180 * atan(Sinh / sqrt(1 - Sinh * Sinh))) / pi
	solalt
}

solazi <- function(localtime,Lat,Long,Julian,merid=0,dst=0)
{
  stime<-solartime(localtime,Long,Julian,merid,dst)
	tt = 0.261799 * (stime - 12)
	declin = (pi * 23.5 / 180) * cos(2 * pi * ((Julian - 171) / 365.25))
	Sinh = sin(declin) * sin(Lat * pi / 180) + cos(declin) * cos(Lat * pi / 180) * cos(tt)
	hh = (atan(Sinh / sqrt(1 - Sinh * Sinh)))
	Sinazi = cos(declin) * sin(tt) / cos(hh)
	cosazi = (sin(Lat * pi / 180) * cos(declin) * cos(tt) - cos(pi * Lat / 180) * sin(declin)) / sqrt((cos(declin) *
            sin(tt)) ^ 2 + (sin(pi * Lat / 180) * cos(declin) * cos(tt) - cos(pi * Lat / 180) * sin(declin)) ^ 2)
	solazi = 180 + (180 * atan(Sinazi / sqrt(1 - Sinazi * Sinazi))) / pi
	if (cosazi < 0) {
		if (Sinazi < 0) {
		solazi = 180 - solazi
		} else {
		solazi = 540 - solazi
		}
	}
	solazi
}

horizonangle <- function(dtm,azimuth,res=100)
{
  dtm<-(dtm*5)/res
  azimuth<-azimuth-90
  azi <- azimuth * (pi/180)
	horizon <- array(0,dim(dtm))
	dtm3 <- array(0,dim(dtm)+200)
	x <- dim(dtm)[1]
	y <- dim(dtm)[2]
	dtm3[101:(x+100),101:(y+100)] <- dtm
	for (step in 1:10) {
		horizon[1:x,1:y] <- pmax(horizon[1:x,1:y], (dtm3[(101+sin(azi)*step^2):(x+100+sin(azi)*step^2),(101+cos(azi)*step^2):(y+100+cos(azi)*step^2)]-dtm3[101:(x+100),101:(y+100)])/(5*step^2))
	}
horizon
}

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

skyview <-function(dtm,steps=36)
{
	sky <- array(1,dim(dtm))
	for (s in 1:steps) {
		sky <- sky-atan(horizonangle(dtm,s*360/steps))/((pi/2)*steps)
	}
sky
}

####################################################################################
# Direct normal radiation: to downscaled direct radiation
####################################################################################
    
radiation_downscale<-function(day,month,year,hr,demuk,dem.buffer,dem.block,slope.buffer,aspect.buffer,dir_rad,print.results=TRUE,write.files=TRUE)
{
    e.demuk<-extent(demuk)
    mn<-0
    datetime<-paste(year,"/",sprintf("%02d",month,sep=""),"/",sprintf("%02d",day,sep="")," ",sprintf("%02d",hr,sep=""),":",sprintf("%02d",mn,sep=""),sep="")
    print(datetime)
    dir_direct<-paste(dir_rad,"rasters/direct/",sep="")
    dir_diffuse<-paste(dir_rad,"rasters/diffuse/",sep="")
    dir_total<-paste(dir_rad,"rasters/total/",sep="")

    # ensure SIS file read in matches DNR file
    infile.dnr<-paste(dir_dni,"DNIhm",year,sprintf("%02d",month,sep=""),sprintf("%02d",day,sep=""),sprintf("%02d",hr,sep=""),sprintf("%02d",mn,sep=""),"002UD1000101UD.nc",sep="")
    infile.sis<-paste(dir_sis,"SIShm",year,sprintf("%02d",month,sep=""),sprintf("%02d",day,sep=""),sprintf("%02d",hr,sep=""),sprintf("%02d",mn,sep=""),"002UD1000101UD.nc",sep="")
    
    # Read in data from ncdf file
    ncdf_dnr<-nc_open(infile.dnr)
    ncdf_sis<-nc_open(infile.sis)
    dnr.r<-raster(infile.dnr)
    sis.r<-raster(infile.sis)
    #e<-extent(c(-7,-2,49,52)) # IMPORTANT: initial cropping to area of interest - CHECK & cf with demuk cropping later ?? !!!!
    #dnr.r<-crop(dnr.r,e)
    #sis.r<-crop(sis.r,e)
    projection(dnr.r)<-"+init=epsg:4326"
    projection(sis.r)<-"+init=epsg:4326"
    
    # reproject to OSGB and set extent to same as DEM
    dnr.rprj<-projectRaster(dnr.r,crs="+init=epsg:27700")
    sis.rprj<-projectRaster(sis.r,crs="+init=epsg:27700")
    dnr.rprj<-crop(dnr.rprj,e.demuk)
    sis.rprj<-crop(sis.rprj,e.demuk)
    
    if (cellStats(sis.rprj,max)>0) # if DAYTIME
    {
      # resample to 100m
      dnr.buffer<-raster::resample(dnr.rprj,dem.buffer)
      sis.buffer<-raster::resample(sis.rprj,dem.buffer)
      #par(mfrow=c(2,2))
      #plot(dnr.buffer)
      #plot(sis.buffer)
      
      # work out Julian day and time - ncdf time in hrs since 1/1/1983 0:00
      jul.base<-JDdoy(1,1983)
      hrs<-ncvar_get(ncdf_dnr,"time")
      days<-floor(hrs/24)
      jul.day<-as.numeric(jul.base+days) #JD at 12:00 on ncdf day
      h<-as.numeric(hrs%%24)
      if (h!=hr) {print("WARNING h^= hr!!!")}     
      
      # creates raster template for storing direct and diffuse radiation values with cell values of -9
      direct.r<-dem.buffer*0-9
      diffuse.r<-dem.buffer*0-9
      total.r<-dem.buffer*0-9

      #plot(dem.buffer)
      #plot(slope.buffer)
      #plot(aspect.buffer)
      
      # convert values to matrices for use with solar index function
      m.dem<-getValues(dem.buffer,format="matrix") # use.raster function above
      m.slope<-getValues(slope.buffer,format="matrix")
      m.aspect<-getValues(aspect.buffer,format="matrix")
      
      # converts NAs to zeros - NB: buffer of boundary blocks = NA->0
      sel<-which(is.na(m.dem)==T); m.dem[sel]<-0
      sel<-which(is.na(m.slope)==T); m.slope[sel]<-0
      sel<-which(is.na(m.aspect)==T); m.aspect[sel]<-0
      
      # Calculate lat and long of centre of grid
      ll<-OSGBtolatlong(xmin(dem.buffer)+(0.5*(xmax(dem.buffer)-xmin(dem.buffer))) , ymin(dem.buffer)+(0.5*(ymax(dem.buffer)-ymin(dem.buffer))) )
      lat<-as.numeric(ll[2])
      long<-as.numeric(ll[1])
      si<-solarindex(slope=m.slope,aspect=m.aspect,localtime=hr,
            Lat=lat,Long=long,Julian=jul.day,dtm=m.dem)
      
      si.flat<-solarindex(slope=0,aspect=0,localtime=hr,
              Lat=lat,Long=long,Julian=jul.day,shadow=F)[1,1]
      
      # Direct normal radiation:
      dnr.m<-getValues(dnr.buffer,format="matrix")
      # Direct radiation: flat
      dir.flat<-dnr.m*si.flat
      # Diffuse radiation flat
      sis.m<-getValues(sis.buffer,format="matrix")
      dif.m<-sis.m-dir.flat
      
      # Calculate matrix coords that define central 'block'
      xmn<-1+(xmin(dem.block)-xmin(dem.buffer))/res(dem.block)[1]
      xmx<-dim(dem.buffer)[1]-(xmax(dem.buffer)-xmax(dem.block))/res(dem.block)[1]
      ymn<-1+(ymin(dem.block)-ymin(dem.buffer))/res(dem.block)[1]
      ymx<-dim(dem.buffer)[2]-(ymax(dem.buffer)-ymax(dem.block))/res(dem.block)[1]
      
      # downscaled direct radiation
      direct<-dnr.m*si
      direct<-direct[xmn:xmx,ymn:ymx] # Extract block
      # downscaled diffuse radiation
      sv<-skyview(m.dem)
      diffuse<-dif.m*sv
      diffuse<-diffuse[xmn:xmx,ymn:ymx] # Extract block
      # add mask back in so that SIs are only produced for land
      mask<-getValues(dem.block,format="matrix")*0
      direct<-direct+mask
      diffuse<-diffuse+mask
      total<-direct+diffuse
      
      # Extract central Block and Convert to rasters
      direct.r<-raster(direct,template=dem.block)
      diffuse.r<-raster(diffuse,template=dem.block)
      total.r<-raster(total,template=dem.block) 

    } else {         # IF NIGHTIME         
      direct.r<-dem.block*0
      diffuse.r<-dem.block*0
      total.r<-dem.block*0
    }
    
    # Write files
    if (write.files){
      fileout1<-paste(dir_direct,"direct_",year,"_",sprintf("%02d",month,sep=""),"_",
                      sprintf("%02d",day,sep=""),"_",hr,"h.tif",sep="")
      fileout2<-paste(dir_diffuse,"diffuse_",year,"_",sprintf("%02d",month,sep=""),"_",
                      sprintf("%02d",day,sep=""),"_",hr,"h.tif",sep="")
      fileout3<-paste(dir_total,"total_",year,"_",sprintf("%02d",month,sep=""),"_",
                      sprintf("%02d",day,sep=""),"_",hr,"h.tif",sep="")
      writeRaster(direct.r,file=fileout1,overwrite=T)
      writeRaster(diffuse.r,file=fileout2,overwrite=T)
      writeRaster(total.r,file=fileout3,overwrite=T)
    }
    #output<-list(direct.r,diffuse.r,total.r)
    
    #direct.r;diffuse.r;total.r;dem.block;dem.buffer
    #Create raster stack and print if requested
    if (print.results==TRUE){
      result.stack<-stack(direct.r,diffuse.r,total.r,dem.block)
      names(result.stack)<-c("direct","diffuse","total","dem block")
      par(mfrow=c(2,2))
      plot(result.stack)  
    } # end if
    
    results<-c(direct.r,diffuse.r,total.r)
    return(results)
} # end function


# CUT OUTS
#direct.raster<-mosaic(direct.raster,dir.r,fun=max)
#diffuse.raster<-mosaic(diffuse.raster,dif.r,fun=max)
#total.raster<-mosaic(total.raster,tot.r,fun=max)
