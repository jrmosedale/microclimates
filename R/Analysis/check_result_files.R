# Check results files 1 to 935 for specified year exist
# Carson job variables
args <-commandArgs(trailingOnly = TRUE)
print(args)
year <- as.integer(args[1])

root<-"/home/ISAD/jm622/Data2015/"   # Source data and output data
in.root<-"/data/jm622/ModelData/" # model input data 
dir_finalt<-paste(root,"Temp100m/",sep="")

# Check if results files for all cells in year
filelist<-list.files(dir_finalt)
cells<-c(1:935)
searchnames<-paste("block-",sprintf("%03d",cells,sep=""),"-",year,".R",sep="")
results<-searchnames %in% filelist

# Save search
fileout<-paste(dir_finalt,"files_for_",year,".csv",sep="")
print(fileout)
write.csv(results,file=fileout)
