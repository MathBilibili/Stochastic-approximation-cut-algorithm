LoadNewPhi0<-function(numsel=100,filename="PhiC.csv",CutModel){
  px <- CutModel$px
  py <- CutModel$py
  prox <- CutModel$prox
  rprox <- CutModel$rprox
  denT <- CutModel$denT
  tranT <- CutModel$tranT
  Z <- CutModel$Z
  Y <- CutModel$Y
  d_x <- CutModel$d_x
  d_y <- CutModel$d_y
  ncol.print <- function(dat) matrix(as.matrix(dat),ncol=ncol(dat),dimnames=NULL)

  Phistar<-ncol.print(as.matrix(fread(filename,head=F)))

  if(d_y!=1){
    da<-apply(Phistar,max,MARGIN = 1)
    xi<-apply(Phistar,min,MARGIN = 1)
    p<-cbind(Phistar,da,xi)
    stan<-function(v){
      v<-as.vector(v)
      if(v[length(v)]!=v[length(v)-1]){
        out<-(v[1:(length(v)-2)]-v[length(v)])/(v[length(v)-1]-v[length(v)])
      }else{
        out<-rep(0,(length(v)-2))
      }
      return(out)
    }
    Phistan<-t(apply(p,FUN=stan,MARGIN = 1))
  }else{
    Phistan<-Phistar
  }

  #Max-Min process
  MMP<-function(Phi,num_sel=100){
    if(d_y!=1){
      num_sel<-min(num_sel,dim(Phi)[1])
      ind<-seq(1,dim(Phi)[1])
      A<-sample(ind,size=1)
      Ac<-ind[!ind%in%A]
      for(i in 2:num_sel){
        X<-Phi[A,]
        if(is.vector(X)){
          X<-t(as.matrix(X))
        }
        Y<-Phi[Ac,]
        d<-dist2(X,Y)
        mi<-apply(d,min,MARGIN = 2)
        A[i]<-Ac[which.max(mi)]
        Ac<-ind[!ind%in%A]
      }
      return(A)
    }else{
      num_sel<-min(num_sel,dim(Phi)[1])
      ind<-seq(1,dim(Phi)[1])
      A<-sample(ind,size=1)
      Ac<-ind[!ind%in%A]
      for(i in 2:num_sel){
        X<-Phi[A,]
        Y<-Phi[Ac,]
        d<-dist2(X,Y)
        mi<-apply(d,min,MARGIN = 2)
        A[i]<-Ac[which.max(mi)]
        Ac<-ind[!ind%in%A]
      }
      return(A)
    }
  }

  MMP<-cmpfun(MMP)

  #the number of auxiliary \Phi_0 is given by 'num_sel'
  PhiC<-Phistar[MMP(Phistan,num_sel=numsel),] %>% as.matrix()

  return(PhiC)

}
