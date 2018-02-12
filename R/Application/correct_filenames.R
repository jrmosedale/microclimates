# Check file present for each hour of every day between start and end dates using JD
# Takes about 15 sec per year
# Start the clock!
ptm <- proc.time()

dir_content<-"~/Documents/Exeter/Data2015/test/"
#dir_content<-dir_temp
print(dir_content)

filelist<-list.files(path=dir_content,include.dirs=FALSE)
print(paste("Directory ",dir_content," contains ",length(filelist), " files", sep=""))
newlist<-gsub("ACTUAL","Actual",filelist)
print(length(which(substr(filelist,20,25)=="ACTUAL")))
print(length(which(substr(filelist,20,25)=="Actual")))

# Add path for file replacement
dirlist<-rep(dir_content,length(filelist))
filelist<-paste(dirlist,filelist,sep="")
newlist<-paste(dirlist,newlist,sep="")

for (i in 1:length(filelist)){
  if (filelist[i]!=newlist[i]) file.rename(filelist[i],newlist[i])
}

filelist<-list.files(path=dir_content,include.dirs=FALSE)
print(length(which(substr(filelist,20,25)=="ACTUAL")))
print(length(which(substr(filelist,20,25)=="Actual")))           

# Stop the clock
proc.time() - ptm
