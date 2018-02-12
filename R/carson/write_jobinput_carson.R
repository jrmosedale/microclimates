# Write job array variable list 
#1 FOr monthyl jobs 1/1983-12/2014
out.file<-"C:/Data2015/dateinput_months.txt"
dataline<-rep("",((2014-1983)*12))  # = number of job lines
n<-1

for (year in 1983:2014){
  for (month in 1:12) {
      start.date<-paste(1,month,year,sep=",")
      if (month<12) end.jd<-JDdmy(1,month+1,year)-1 else end.jd<-JDdmy(1,1,year+1)-1
      end.date<-paste(DMYjd(end.jd)$day,DMYjd(end.jd)$month,DMYjd(end.jd)$year,sep=",")
      dataline[n]<-paste(start.date,end.date,sep=",")
      n<-n+1
  }    
}

print(dataline)
write(dataline,file=out.file)
