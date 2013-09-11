dms2dd<-function(latlon){
#
# a function to convert DMS to Decimal Degrees (DD)
#
# input: a vector or matrix/data.frame of class character
# exampe of input data
  #       [,1]             [,2]              
#   [1,] "41  03  21.726" "173   01  27.556"
#   [2,] "41  02  57.518" "173   02  59.480"
#   [3,] "41  02  47.859" "173   04  08.851"
#   [4,] "41  02  38.831" "173   05  26.622"
#   [5,] "41  02  32.047" "173   06  38.337"
#
# Sep 6, 2013,Weimin Jiang
# 
  if(is.null(dim(latlon))){
    lstr<-strsplit(latlon,split=" ")
    mtx<-sapply(lstr,FUN=function(x){
      tmp<-as.numeric(x)
      return(tmp[!is.na(tmp)])
    })
  dd<-mtx[1,]+(mtx[2,]+mtx[3,]/60)/60
  } else {
    lstrs<-apply(latlon,2,FUN=function(x){strsplit(x,split=" ")})
    mtxs<-lapply(lstrs,FUN=function(x){
      sapply(x,FUN=function(x){
      tmp<-as.numeric(x)
      return(tmp[!is.na(tmp)])
    })})
    ddlist<-lapply(mtxs,FUN=function(x){
      tmp<-x[1,]+(x[2,]+x[3,]/60)/60
    })
    dd<-do.call(cbind,ddlist)
}
    return(dd)
}

