

for (jd in start.jd:(start.jd+30)){
  daydate<-paste(DMYjd(jd)$day,"/",DMYjd(jd)$month,"/",DMYjd(jd)$year,sep="")
  #load 5km temp file and find min/max
  max.infile<-paste(dir_temp,"MaxTemp_", DMYjd(jd)$year, "-",sprintf("%02d",DMYjd(jd)$month,sep=""),"-", sprintf("%02d",DMYjd(jd)$day,sep=""),"_ACTUAL.txt", sep="")
  min.infile<-paste(dir_temp,"MinTemp_", DMYjd(jd)$year, "-",sprintf("%02d",DMYjd(jd)$month,sep=""),"-", sprintf("%02d",DMYjd(jd)$day,sep=""),"_ACTUAL.txt", sep="")
  day.tmax<-crop(raster(max.infile, layer=1,xmn=-200000, ymn=-200000,nrows=180,ncols=290, res=c(5000,5000), crs="+init=epsg:27700"),dem.block)
  day.tmin<-crop(raster(min.infile, layer=2,xmn=-200000, ymn=-200000,nrows=180,ncols=290, res=c(5000,5000), crs="+init=epsg:27700"),dem.block)
  
  #load DNI
  
  #load Wstr
  
  #load wind.dir
  print(paste(daydate," Tmax=",round(cellStats(day.tmax,mean))," Tmin=",round(cellStats(day.tmin,mean)),sep="")) 
  if(round(cellStats(day.tmax,mean))>25) {print("VERY HOT")}
  if(round(cellStats(day.tmin,mean))<5) {print("VERY COLD")}
}