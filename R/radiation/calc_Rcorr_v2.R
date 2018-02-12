# calculate_Rcorr(jd,hr,sis.buffer,dnr.buffer,dem.buffer,slope.buffer,aspect.buffer)

calculate_Rcorr<-function(jd,hr,sis.buffer,dnr.buffer,dem.buffer,slope.buffer,aspect.buffer,dem.block)
{
  # TEST LINES
  #for(hr in 6:20){
    # load rad files
    #dnr.buffer<-raster(rad.24h[[1]],hr+1)
    #sis.buffer<-raster(rad.24h[[2]],hr+1)
    
    # Calculate jd.h - Julianday incl hour decimal
    h<--12+hr; 
    h.jd<-h/24;
    jd.h<-jd+h.jd;#print(jd.h)
    
    # Calculate lat and long of centre of grid
    ll<-OSGBtolatlong(xmin(dem.buffer)+(0.5*(xmax(dem.buffer)-xmin(dem.buffer))) , ymin(dem.buffer)+(0.5*(ymax(dem.buffer)-ymin(dem.buffer))) )
    lat<-as.numeric(ll[2])
    long<-as.numeric(ll[1])
    
    # Calculate m.dem,slope & aspect
    # convert values to matrices for use with solar index function
    m.dem<-getValues(dem.buffer,format="matrix") # use.raster function above
    m.slope<-getValues(slope.buffer,format="matrix")
    m.aspect<-getValues(aspect.buffer,format="matrix")
    
    # converts NAs to zeros - NB: buffer of boundary blocks = NA->0
    sel<-which(is.na(m.dem)==T); m.dem[sel]<-0
    sel<-which(is.na(m.slope)==T); m.slope[sel]<-0
    sel<-which(is.na(m.aspect)==T); m.aspect[sel]<-0
    
    if (sunvector(jd.h,lat,long,0)[3]>=0) # if DAYTIME 
    {
      # Calculate indexes to adjust R
      # Calculate solar index for direct radiation 
      si<-solarindex(slope=m.slope,aspect=m.aspect,localtime=hr,
                     Lat=lat,Long=long,Julian=jd,dtm=m.dem)
      si.buffer<-raster(si,template=dem.buffer)
      # Calculate skyview for diffuse radiation
      sv<-skyview(m.dem)
      sv.buffer<-raster(sv,template=dem.buffer)
      # Calculate single value solar index for flat slope
      slope.flat<-0;aspect.flat<-0
      si.flat<-solarindex(slope=slope.flat,aspect=aspect.flat,localtime=hr,
                          Lat=lat,Long=long,Julian=jd,shadow=F)[1,1]
      
      # Calculate reference R values for flat surface (i.e. varies only with cloud cover)
      # OR JUST USE SIS!
      dir.flat.buffer<-si.flat*dnr.buffer # Rdir on flat surface = SID
      dif.flat.buffer<-sis.buffer-dir.flat.buffer #Rdif on flat surface
      #total.flat.buffer<-dif.flat.buffer+dir.flat.buffer
      
      # Calculate R values adjusted to DEM
      dir.buffer<-si.buffer*dnr.buffer # Rdir adjusted for slope and aspect and hillshading
      dif.buffer<-sv.buffer*dif.flat.buffer
      total.buffer<-dir.buffer+dif.buffer
      
      # Calculate correction ratio 
      # Calculate mean for whole block - in effect this might be better as mean of wider area?
      #total.flat.block<-mask(dir.flat.buffer,dem.block)
      sis.block<-mask(crop(sis.buffer,dem.block),dem.block)
      Rref<-cellStats(sis.block,stat="mean") # single value for whole block
      R100<-mask(crop(total.buffer,dem.block),dem.block)
      
      #par(mfrow=c(2,2))
      #plot(sis.block,main="SIS - input")
      #plot(mask(crop(dir.buffer,dem.block),dem.block),main="Direct Rad (si corrected)")
      plot(mask(crop(total.buffer,dem.block),dem.block),main="Total Rad (si,sv corrected)")
      plot(Rcorr.block,main=paste("Rcorr value at ",hr,":00",sep=""))
  #}
    } else  { Rcorr.block<-raster( array(1,c(dim(dem.block)[1],dim(dem.block)[2])),template=dem.block ) } # If NIGHTIME returns Rcor of 1.0
    Rcorr.block<-c(Rref,R100)  
  return(Rcorr.block)
  # Calculate modified temperature at 100m
  #t.block<-tref.block*Rcorr.block
  #plot(t.block,main=paste("Temperature at ",hr,":00",sep=""))
} # end function