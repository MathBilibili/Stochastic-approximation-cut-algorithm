library(spatstat)
library(compiler)
library(truncnorm)
#data specification

task_id_string <- Sys.getenv("SLURM_ARRAY_TASK_ID")
task_id <- as.numeric(task_id_string)

data<-read.csv("simu_data.csv")
Y<-apply(data,1,mean)
N<-dim(data)[1]
n<-dim(data)[2]
Z<-rep(0,length(Y))
for(i in 1:N){
  for(j in 1:n){
    Z[i]<-Z[i]+(data[i,j]-Y[i])^2
  }
}



#density 1 Z|phi \times phi
px<-function(phi,Z){
  PbZ<-rbind(phi,Z)
  out<-sum(apply(PbZ,FUN=function(x){log(((x[1])^(-(n+1)/2))*exp(-x[2]/(2*x[1])))},MARGIN = 2))
  return(out)
}
#unnormalizing density2 Y|theta&phi \times theta
py<-function(Y,theta,phi){
  PbY<-rbind(phi,Y)
  out<-log(1/theta+mean(phi)/n)+sum(apply(PbY,FUN=function(x){log((1/sqrt(theta+x[1]/n))*exp(-(x[2]^2)/(2*(theta+x[1]/n))))},MARGIN = 2))
  return(out)
}

#proposal 1
prox<-function(phi_n,phi){
  out<-sum(log(dtruncnorm(phi_n, a=0, b=Inf, mean = phi, sd = 0.1))) 
  return(out) 
}




#proposal 1 samling function
rprox<-function(phi){
  out<-rtruncnorm(1, a=0, b=Inf, mean = phi, sd = 0.1)
  return(out)
}





#Choice of Auxiliary Parameter Set chain
init<-list(phi=rep(1,N))
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


# #standardise element of phi
# Phistar<-MAux(init,Z,Y,num_run=5000,burn_in=0,thin=1,r=0.01)
# da<-apply(Phistar,max,MARGIN = 1)
# xi<-apply(Phistar,min,MARGIN = 1)
# p<-cbind(Phistar,da,xi)
# stan<-function(v){
#   v<-as.vector(v)
#   if(v[length(v)]!=v[length(v)-1]){
#     out<-(v[1:(length(v)-2)]-v[length(v)])/(v[length(v)-1]-v[length(v)])
#   }else{
#     out<-rep(0,(length(v)-2))
#   }
#   return(out)
# }
# Phistan<-t(apply(p,FUN=stan,MARGIN = 1))
# 
# #Max-Min process
# MMP<-function(Phi,num_sel=100){
#   num_sel<-min(num_sel,dim(Phi)[1])
#   ind<-seq(1,dim(Phi)[1])
#   A<-sample(ind,size=1)
#   Ac<-ind[!ind%in%A]
#   for(i in 2:num_sel){
#     X<-Phi[A,1]
#     Y<-Phi[A,2]
#     x2<-Phi[Ac,1]
#     y2<-Phi[Ac,2]
#     d<-crossdist(X,Y,x2,y2)
#     mi<-apply(d,min,MARGIN = 1)
#     A[i]<-Ac[which.max(mi)]
#     Ac<-ind[!ind%in%A]
#   }
#   return(A)
# }
# 
# MMP<-cmpfun(MMP)
# 
# 
# PhiC<-Phistar[MMP(Phistan,num_sel=500),]
# 
# write.csv(PhiC,"AuxPhi.csv")

Phistar<-read.csv("AuxPhi.csv",head=F)
Phistar<-as.matrix(Phistar)
ncol.print <- function(dat) matrix(as.matrix(dat),ncol=ncol(dat),dimnames=NULL)
Phistar<-ncol.print(Phistar)

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

#Max-Min process
MMP<-function(Phi,num_sel=100){
  num_sel<-min(num_sel,dim(Phi)[1])
  ind<-seq(1,dim(Phi)[1])
  A<-sample(ind,size=1)
  Ac<-ind[!ind%in%A]
  for(i in 2:num_sel){
    X<-Phi[A,1]
    Y<-Phi[A,2]
    x2<-Phi[Ac,1]
    y2<-Phi[Ac,2]
    d<-crossdist(X,Y,x2,y2)
    mi<-apply(d,min,MARGIN = 1)
    A[i]<-Ac[which.max(mi)]
    Ac<-ind[!ind%in%A]
  }
  return(A)
}

MMP<-cmpfun(MMP)


PhiC<-Phistar[MMP(Phistan,num_sel=70),]


#internal invariant distribution
p_inv<-function(t,phi,phist,wst){
  p_inv_inp<-cbind(phist,wst)
  inp_in<-dim(p_inv_inp)[2]-1
  psi<-function(x){
    outtt<-sign(identical(as.vector(x[1:inp_in]),phi))
    return(outtt)
  }
  out<-which(apply(p_inv_inp,MARGIN = 1,FUN=psi)==1)
  outtt<-as.numeric(py(Y,t,p_inv_inp[out,][1:inp_in])-log(p_inv_inp[out,][inp_in+1]))
  return(outtt)
}

#proposal sampling for t and phi
pro_tp<-function(P,pro_tp_inp){
  t<-pro_tp_inp$t
  phi<-pro_tp_inp$phi
  ran<-runif(1,0,1)
  if(ran<P){
    t_n<-rtruncnorm(1,a=0, b=Inf,mean = t,sd=0.1)      #proposal for t may be changed
    out<-list(t=t_n,phi=phi)
  }else{
    phi_n<-PhiC[sample(seq(1,dim(PhiC)[1]),1),]  
    out<-list(t=t,phi=phi_n)
  }
  return(out)
}

#proposal density for t and phi
dpro_tp<-function(P,x_n,x){
  t<-x$t
  t_n<-x_n$t
  phi<-x$phi
  phi_n<-x_n$phi
  if(identical(t,t_n)){
    out<-1/(dim(PhiC)[1])*(1-P)
  }else{
    out<-dtruncnorm(t_n,a=0, b=Inf,mean = t,sd=0.1)*P        #proposal for t may be changed
  }
  return(out)
}


#internal chain
no<-1000
MA_in<-function(H,W,n,pai,n_trun){
  t<-H$t
  I<-H$I
  P<-0.25
  inp<-list(t=t,phi=PhiC[I,])
  ou<-pro_tp(P,inp)
  rate<-p_inv(ou$t,ou$phi,PhiC,W)+log(dpro_tp(P,inp,ou))-p_inv(t,PhiC[I,],PhiC,W)-log(dpro_tp(P,ou,inp))
  alfa<-min(1,exp(rate))
  rpan<-runif(1)
  t_n<-ou$t*sign(rpan<=alfa)+inp$t*sign(rpan>alfa)
  phi_n<-ou$phi*sign(rpan<=alfa)+inp$phi*sign(rpan>alfa)
  iden<-function(x){
    Ou<-identical(x,phi_n)
    return(Ou)
  }
  I_n<-which(apply(PhiC,FUN = iden,MARGIN = 1))
  en<-rep(0,dim(PhiC)[1])
  en[I_n]<-1
  W_n<-exp(log(W)+(no/max(no,n))*(en-pai))        #no is n_0
  if(any(W_n>(10^(100+n_trun)))){
    W_n<-rep(1,length(W))            #truncated here
    n_trun<-n_trun+1
  }
  out<-list(t=t_n,I=I_n,W=W_n,n_trun=n_trun)
  return(out)
}

MA_in<-cmpfun(MA_in)

Digset<-c(3,4,10)
#sig_dig<-10                   #significant digit used to round t for high efficiency.
sig_dig<-Digset[task_id]

GRset<-seq(1,10,0.5)          #initial value set for separate chain




#construct sufficient auxiliary set
#init<-list(theta=2,phi=rep(2,N),t=2,I=1) 
init<-list(theta=GRset[task_id],phi=rep(GRset[task_id],N),t=2,I=1)          #t should be inversed vector.
MA_aux<-function(init,Z,Y,PhiC,num_run=1000,burn_in=500){
  theta<-init$theta
  phi<-init$phi
  n_trun<-1
  H<-list(t=init$t,I=init$I)
  W<-rep(1,dim(PhiC)[1])
  ColH<-list(t=init$t,I=init$I,w=W[1])      #Comulative information
  pai<-rep(1/length(W),length(W))          # sampling frequency
  Count_Tt<-1
  for(i in 1:num_run){
    if(i<=burn_in){
      InR<-MA_in(H,W,n=i,pai,n_trun)
      W<-InR$W
      H$t<-InR$t
      H$I<-InR$I
      n_trun<-InR$n_trun
    }else{
      j<-i-burn_in
      InR<-MA_in(H,W,n=i,pai,n_trun)
      W<-InR$W
      H$t<-InR$t
      H$I<-InR$I
      n_trun<-InR$n_trun
      
      ColH$t<-round(InR$t,digits = sig_dig)
      ColH$I<-InR$I
      ColH$w<-W[InR$I]
      #calculate sampling frequency
      bas<-rbind(ColH$I,ColH$w,ColH$t)
      if(j==1){
        Tt<-as.matrix(ColH$t)
      }else{
        Tt<-as.matrix(unique(cbind(Tt,ColH$t),MARGIN=2))
      }
      
      phi_n<-rprox(phi)
      if(j==1){
        deno<-bas[2]*exp(py(Y,ColH$t,phi_n))/exp(py(Y,ColH$t,PhiC[bas[1],]))
        numr<-bas[2]*exp(py(Y,ColH$t,phi_n))/exp(py(Y,ColH$t,PhiC[bas[1],]))
        fenzi_o<-exp(py(Y,Tt,phi_n))
        Ptau<-numr/deno
      }else{
        fenzi<-function(t){   
          return(exp(py(Y,t,phi_n)))
        }
        
        if(Count_Tt==dim(Tt)[2]){
          iidentical<-function(x){
            return(identical(x,ColH$t))
          }
          if(dim(Tt)[2]>1){
            fenzi_n<-apply(Tt,FUN = fenzi,MARGIN = 2)  
          }else{
            fenzi_n<-exp(py(Y,Tt,phi_n))
          }
          
          nchan<-which(apply(Tt,FUN=iidentical,MARGIN = 2))
          numr<-numr/fenzi_o
          numr<-numr*fenzi_n
          nchanadd<-bas[2]*exp(py(Y,ColH$t,phi_n))/exp(py(Y,ColH$t,PhiC[bas[1],]))
          numr[nchan]<-numr[nchan]+nchanadd
          fenzi_o<-fenzi_n
          deno<-sum(numr)
          Ptau<-numr/deno
        }else{
          fenzi_n<-apply(Tt,FUN = fenzi,MARGIN = 2)  
          nchanadd<-length(fenzi_n)
          numr<-numr/fenzi_o
          numr<-numr*fenzi_n[-nchanadd]
          numr[nchanadd]<-bas[2]*exp(py(Y,ColH$t,phi_n))/exp(py(Y,ColH$t,PhiC[bas[1],]))
          fenzi_o<-fenzi_n
          deno<-sum(numr)
          Ptau<-numr/deno
        }
      }
      Count_Tt<-dim(Tt)[2]
      
      tau_n<-Tt[,sample(x=seq(1,dim(Tt)[2]),size = 1,prob = Ptau)]
      if(length(tau_n)==1){
        theta_n<-runif(1,min=tau_n-5*10^(-sig_dig-1),max=tau_n+5*10^(-sig_dig-1))  #uniformly drawing \theta
      }else{
        for(q in 1:length(tau_n)){
          theta_n[q]<-runif(1,min=tau_n[q]-5*10^(-sig_dig-1),max=tau_n[q]+5*10^(-sig_dig-1))
        }
      }
      rate<-px(phi_n,Z)+prox(phi,phi_n)-px(phi,Z)-prox(phi_n,phi)
      alfa<-min(1,exp(rate))
      rpan<-runif(1)
      phi<-phi_n*sign(rpan<=alfa)+phi*sign(rpan>alfa)
      theta<-theta_n*sign(rpan<=alfa)+theta*sign(rpan>alfa)
    }
    print(c(i,InR$t))
  }
  MA_aux_out<-list(Tt=Tt,Ptau=Ptau,fenzi_o=fenzi_o,numr=numr,t=H$t,I=H$I,n_trun=n_trun,aux_num_run=num_run,W=W,theta=theta,phi=phi)
  return(MA_aux_out)
}

MA_aux<-cmpfun(MA_aux)

MA_aux_out<-MA_aux(init,Z,Y,PhiC,num_run = 11000,burn_in=10000)










#external chain
init<-list(theta=MA_aux_out$theta,phi=MA_aux_out$phi,t=MA_aux_out$t,I=MA_aux_out$I)          #t should be inversed vector.
MA_ex<-function(Aux_Tt,init,Z,Y,PhiC,num_run=1000,burn_in=500,thin=1){
  theta<-init$theta
  phi<-init$phi
  n_trun<-MA_aux_out$n_trun
  H<-list(t=init$t,I=init$I)
  W<-Aux_Tt$W
  ColH<-list(t=init$t,I=init$I,w=W[1])      #Comulative information
  pai<-rep(1/length(W),length(W))          # sampling frequency
  sto.phi<-as.matrix(rep(phi,((num_run-burn_in)/thin))) #note the dimension of phi
  dim(sto.phi)<-c(((num_run-burn_in)/thin),length(phi))
  sto.theta<-as.matrix(rep(theta,((num_run-burn_in)/thin))) #note the dimension of theta
  dim(sto.theta)<-c(((num_run-burn_in)/thin),length(theta))
  sto.time<-rep(0,((num_run-burn_in)/thin))
  Tt<-Aux_Tt$Tt
  Ptau<-Aux_Tt$Ptau
  fenzi_o<-Aux_Tt$fenzi_o
  numr<-Aux_Tt$numr
  InRadd<-Aux_Tt$aux_num_run
  Count_Tt<-dim(MA_aux_out$Tt)[2]
  for(i in 1:num_run){
    if((i<=burn_in)|(i%%thin!=0)){
      InR<-MA_in(H,W,n=i+InRadd,pai,n_trun)
      W<-InR$W
      H$t<-InR$t
      H$I<-InR$I
      n_trun<-InR$n_trun
      
      ColH$t<-round(InR$t,digits = sig_dig)
      ColH$I<-InR$I
      ColH$w<-W[InR$I]
      #calculate sampling frequency
      bas<-rbind(ColH$I,ColH$w,ColH$t)
      Tt<-as.matrix(unique(cbind(Tt,ColH$t),MARGIN=2))
      phi_n<-rprox(phi)
      
      fenzi<-function(t){   
        return(exp(py(Y,t,phi_n)))
      }
      
      if(Count_Tt==dim(Tt)[2]){
        iidentical<-function(x){
          return(identical(x,ColH$t))
        }
        if(dim(Tt)[2]>1){
          fenzi_n<-apply(Tt,FUN = fenzi,MARGIN = 2)  
        }else{
          fenzi_n<-exp(py(Y,Tt,phi_n))
        }
        
        nchan<-which(apply(Tt,FUN=iidentical,MARGIN = 2))
        numr<-numr/fenzi_o
        numr<-numr*fenzi_n
        nchanadd<-bas[2]*exp(py(Y,ColH$t,phi_n))/exp(py(Y,ColH$t,PhiC[bas[1],]))
        numr[nchan]<-numr[nchan]+nchanadd
        fenzi_o<-fenzi_n
        deno<-sum(numr)
        Ptau<-numr/deno
      }else{
        fenzi_n<-apply(Tt,FUN = fenzi,MARGIN = 2)  
        nchanadd<-length(fenzi_n)
        numr<-numr/fenzi_o
        numr<-numr*fenzi_n[-nchanadd]
        numr[nchanadd]<-bas[2]*exp(py(Y,ColH$t,phi_n))/exp(py(Y,ColH$t,PhiC[bas[1],]))
        fenzi_o<-fenzi_n
        deno<-sum(numr)
        Ptau<-numr/deno
      }
      Count_Tt<-dim(Tt)[2]
      
      tau_n<-Tt[,sample(x=seq(1,dim(Tt)[2]),size = 1,prob = Ptau)]
      if(length(tau_n)==1){
        theta_n<-runif(1,min=tau_n-5*10^(-sig_dig-1),max=tau_n+5*10^(-sig_dig-1))  #uniformly drawing \theta
      }else{
        for(q in 1:length(tau_n)){
          theta_n[q]<-runif(1,min=tau_n[q]-5*10^(-sig_dig-1),max=tau_n[q]+5*10^(-sig_dig-1))
        }
      }
      rate<-px(phi_n,Z)+prox(phi,phi_n)-px(phi,Z)-prox(phi_n,phi)
      alfa<-min(1,exp(rate))
      rpan<-runif(1)
      phi<-phi_n*sign(rpan<=alfa)+phi*sign(rpan>alfa)
      theta<-theta_n*sign(rpan<=alfa)+theta*sign(rpan>alfa)
    }else{
      earlier_time<-Sys.time()
      InR<-MA_in(H,W,n=i+InRadd,pai,n_trun)
      W<-InR$W
      H$t<-InR$t
      H$I<-InR$I
      n_trun<-InR$n_trun
      
      ColH$t<-round(InR$t,digits = sig_dig)
      ColH$I<-InR$I
      ColH$w<-W[InR$I]
      #calculate sampling frequency
      bas<-rbind(ColH$I,ColH$w,ColH$t)
      Tt<-as.matrix(unique(cbind(Tt,ColH$t),MARGIN=2))
      phi_n<-rprox(phi)
      fenzi<-function(t){   
        return(exp(py(Y,t,phi_n)))
      }
      
      if(Count_Tt==dim(Tt)[2]){
        iidentical<-function(x){
          return(identical(x,ColH$t))
        }
        if(dim(Tt)[2]>1){
          fenzi_n<-apply(Tt,FUN = fenzi,MARGIN = 2)  
        }else{
          fenzi_n<-exp(py(Y,Tt,phi_n))
        }
        
        nchan<-which(apply(Tt,FUN=iidentical,MARGIN = 2))
        numr<-numr/fenzi_o
        numr<-numr*fenzi_n
        nchanadd<-bas[2]*exp(py(Y,ColH$t,phi_n))/exp(py(Y,ColH$t,PhiC[bas[1],]))
        numr[nchan]<-numr[nchan]+nchanadd
        fenzi_o<-fenzi_n
        deno<-sum(numr)
        Ptau<-numr/deno
      }else{
        fenzi_n<-apply(Tt,FUN = fenzi,MARGIN = 2)  
        nchanadd<-length(fenzi_n)
        numr<-numr/fenzi_o
        numr<-numr*fenzi_n[-nchanadd]
        numr[nchanadd]<-bas[2]*exp(py(Y,ColH$t,phi_n))/exp(py(Y,ColH$t,PhiC[bas[1],]))
        fenzi_o<-fenzi_n
        deno<-sum(numr)
        Ptau<-numr/deno
      }
      Count_Tt<-dim(Tt)[2]
      
      tau_n<-Tt[,sample(x=seq(1,dim(Tt)[2]),size = 1,prob = Ptau)]
      if(length(tau_n)==1){
        theta_n<-runif(1,min=tau_n-5*10^(-sig_dig-1),max=tau_n+5*10^(-sig_dig-1))  #uniformly drawing \theta
      }else{
        for(q in 1:length(tau_n)){
          theta_n[q]<-runif(1,min=tau_n[q]-5*10^(-sig_dig-1),max=tau_n[q]+5*10^(-sig_dig-1))
        }
      }
      rate<-px(phi_n,Z)+prox(phi,phi_n)-px(phi,Z)-prox(phi_n,phi)
      alfa<-min(1,exp(rate))
      rpan<-runif(1)
      phi<-phi_n*sign(rpan<=alfa)+phi*sign(rpan>alfa)
      theta<-theta_n*sign(rpan<=alfa)+theta*sign(rpan>alfa)
      sto.phi[((i-burn_in)/thin),]<-phi
      sto.theta[((i-burn_in)/thin),]<-theta
      recent_time<-Sys.time()
      diff_time<-difftime(recent_time,earlier_time,tz="GMT",units="secs")
      sto.time[((i-burn_in)/thin)]<-diff_time
    }
    if(i %in% seq(10000,100000,1000)){
      Tem_out<-list(phi=sto.phi[1:((i-burn_in)/thin),],theta=sto.theta[1:((i-burn_in)/thin),],time=sto.time[1:((i-burn_in)/thin)],Tt=Tt[1,],Ptau_fenmu=numr/fenzi_o)
      Tem_RR<-data.frame(phi_ten= Tem_out$phi[,10],phi_one= Tem_out$phi[,1],theta= Tem_out$theta,time= Tem_out$time)
      write.csv(Tem_RR,paste("Result_",task_id,".csv",sep = ""))
      Tem_Pr<-data.frame(Tt=Tem_out$Tt,Ptau=Tem_out$Ptau_fenmu)
      write.csv(Tem_Pr,paste("TauProbability_",task_id,".csv",sep = ""))
    }
    print(c(i,rate,alfa,theta))
  }
  OUT<-list(phi=sto.phi,theta=sto.theta,Tt=Tt[1,],Ptau=Ptau,time=sto.time,Ptau_fenmu=numr/fenzi_o)
  return(OUT)
}

MA_ex<-cmpfun(MA_ex)

Result<-MA_ex(MA_aux_out,init,Z,Y,PhiC,num_run=100000,burn_in=0,thin=10)


RR<-data.frame(phi_ten=Result$phi[,10],phi_one=Result$phi[,1],theta=Result$theta,time=Result$time)
TauPro<-data.frame(Tt=Result$Tt,Ptau=Result$Ptau_fenmu)
write.csv(RR,paste("Result_",task_id,".csv",sep = ""))
write.csv(TauPro,paste("TauProbability_",task_id,".csv",sep = ""))
