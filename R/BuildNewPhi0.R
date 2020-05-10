BuildNewPhi0<-function(numrun=1000,burnin=500,numsel=100,CutModel){
  px <- CutModel$px
  py <- CutModel$py
  prox <- CutModel$prox
  rprox <- CutModel$rprox
  denT <- CutModel$proy
  tranT <- CutModel$rproy
  Z <- CutModel$Z
  Y <- CutModel$Y
  d_x <- CutModel$d_x
  d_y <- CutModel$d_y

  init<-list(phi=rep(0.5,d_y))
  MAux<-function(init,Z,num_run=1000,burn_in=500,thin=1){
    phi<-init$phi
    sto.phi<-as.matrix(rep(phi,((num_run-burn_in)/thin))) #note the dimension of phi
    dim(sto.phi)<-c(((num_run-burn_in)/thin),length(phi))
    for(i in 1:num_run){
      if((i<=burn_in)|(i%%thin!=0)){
        phi_n<-rprox(phi)
        rate<-px(phi_n,Z)+prox(phi,phi_n)-px(phi,Z)-prox(phi_n,phi)
        alfa<-min(1,exp(rate))
        rpan<-runif(1)
        phi<-phi_n*sign(rpan<=alfa)+phi*sign(rpan>alfa)
      }else{
        phi_n<-rprox(phi)
        rate<-px(phi_n,Z)+prox(phi,phi_n)-px(phi,Z)-prox(phi_n,phi)
        alfa<-min(1,exp(rate))
        rpan<-runif(1)
        phi<-phi_n*sign(rpan<=alfa)+phi*sign(rpan>alfa)
        sto.phi[((i-burn_in)/thin),]<-phi
      }
    }
    out<-sto.phi[duplicated.matrix(sto.phi,MARGIN = 1)!=T,]
    return(out)
  }
  MAux<-cmpfun(MAux)     #byte compile



  #standardise element of phi
  Phistar<-MAux(init,Z,num_run=numrun,burn_in=burnin,thin=1)
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

  if(is.vector(Phistan)){Phistar<-as.matrix(Phistar);Phistan<-as.matrix(Phistan)}

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
  PhiC<-Phistar[MMP(Phistan,num_sel=numsel),]

  #store PhiC for future use
  write.table(PhiC,"PhiC.csv",row.names = F,col.names = F,sep = ',')

}
