dd2dms=function(x) {
  # a funciton to convert latlons from decimal degree to degree minute and seconds
  library(stringr)
  latlon<-apply(x,2,as.character)
  out<-apply(latlon,2,FUN=function(x){
  dot<-str_locate(x,pattern="\\.")[,1]
  dd<-as.numeric(substr(x,1,dot-1))
  m<-as.numeric(substr(x,dot,nchar(x)))*60
  msc<-as.character(m)
  dot.m<-str_locate(msc,pattern="\\.")[,1]
  mm<-trunc(m)
  mm<-ifelse(mm<10,paste("0",mm,sep=""),mm)
  ss<-as.numeric(substr(msc,dot.m,nchar(msc)))*60
  paste(dd,mm,ss,sep=" ")
})
  out
}
