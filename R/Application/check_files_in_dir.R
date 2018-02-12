# Check file present for each hour of every day between start and end dates using JD
# Takes about 15 sec per year
# Start the clock!
#ptm <- proc.time()

dir_content<-dir_cal
print(dir_cal)
start<-JDdmy(1,1,1992)
end<-JDdmy(31,12,1992)
filelist<-list.files(path=dir_content,include.dirs=FALSE)
print(paste("Directory ",dir_content," contains ",length(filelist), " files", sep=""))
print(paste("Testing presence of files between ",DMYjd(start)$day,"/",DMYjd(start)$month,"/",DMYjd(start)$year," and ",DMYjd(end)$day,"/",DMYjd(end)$month,"/",DMYjd(end)$year,sep=""))
missing<-0
for (i in start:end){
  for (hr in 0:23){
    day<-DMYjd(i)$day; month<-DMYjd(i)$month; year<-DMYjd(i)$year
    searchname<-paste("CALhm",year,sprintf("%02d",month,sep=""),sprintf("%02d",day,sep=""),sprintf("%02d",hr,sep=""),"00002UD1000101UD.nc",sep="")
    #e.g. CALhm199305082000002UD1000101UD.nc
    # To print feedback on each file checked use:
    #if (searchname %in% filelist) print(paste(searchname," is present",sep="")) else print(paste("Missing file: ",searchname,sep=""))
    # To print summary information ONLY and save names of missing files use:
    if (!(searchname %in% filelist)) {
      #print(paste("Missing file: ",searchname,sep=""))
      missing<-missing+1
    }
}
}
print(paste("Total number of files missing= ",missing," from a total of ",(1+end-start)*24,sep=""))
# Stop the clock
#proc.time() - ptm
