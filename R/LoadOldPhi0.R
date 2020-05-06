LoadOldPhi0<-function(filename="PhiC.csv"){
  ncol.print <- function(dat) matrix(as.matrix(dat),ncol=ncol(dat),dimnames=NULL)
  PhiC<-ncol.print(as.matrix(fread(filename,head=F)))
  return(PhiC)
}
