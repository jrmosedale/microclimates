# Test code for accessing web data
library(XML)
library(RCurl)
library(httr) # install httr package and then use list
install.packages(  " ̃/Downloads/RHTMLForms_0.6-0.tar.gz",
  repos=NULL,
  type="source")
  
# Install from http://www.omegahat.org/RHTMLForms/

library(RHTMLForms) 






# Vineyard application - download table from website
ew_vineyards_url<-url("http://www.englishwine.com/vineyards.htm")

doc1<-xmlTreeParse("http://www.englishwine.com/vineyards.htm")
doc2<-xmlTreeParse("http://www.englishwine.com/vineyards.htm",isHTML=TRUE)
doc3<-htmlTreeParse("http://www.englishwine.com/vineyards.htm")

# Extract vineyard data
vineyards<-readHTMLTable("http://www.englishwine.com/vineyards.htm",which=1,
              header=TRUE, 
              colClasses=c("character","character","character","numeric","character"), 
              trim=TRUE,
              stringsAsFactors=FALSE)
out.file<-"~/Documents/Exeter/Data2015/Vineyards/ewvineyards.csv"
write.table(vineyards,file = out.file,sep=",")
# Submit data in HTML Form
# Using RCurl
postForm(uri, ..., .params = list(), .opts = curlOptions(url = uri),
         curl = getCurlHandle(), style = 'HTTPPOST')
getForm("https://google.co.uk/", .params = c(q="Wine and Climate Change", btnG="Search"))
getURL("https://google.co.uk")
