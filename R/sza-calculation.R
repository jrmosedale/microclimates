
#Create raster for a given date and time describing sza etc

#### FUNCTION for calculating sza using sunvector and sunpos in insol package

sza.angle<- function(jd,lt,ln,tz=0){
  sv<-sunvector(jd,lt,ln,tz)
  sza<-sunpos(sv)[,2]
  return(sza)
}

# Using insol 
day<-30;month<-6;year<-2013; hr<-13
lt<-49.5;ln<-0
lt2<-51;ln2<-0

sza.angle(jd,lt,ln)

#Create test arrays of long, latt, jd
test.lt<-seq(49,52,by=0.2)
test.ln<-seq(-1,0,length=16)
jd.a<-rep(jd,16)

sza.angle(jd.a,test.lt,test.ln)





#### Convert N/E to lat long
# Convert from Eastings and Northings to Latitude and Longitude
 <-  spTransform(data, CRS (latlong))

# we also need to rename the columns
colnames(GP_SP_LL@coords)[colname(GP_SP_LL@coords) =="Easting"] <-"Longitude"
colnames(GP_SP_LL@coords)[colnames(GP_SP_LL@coords)=="Northing"] <-"Latitude"

# Could create sza matrix for every day/hr of year 