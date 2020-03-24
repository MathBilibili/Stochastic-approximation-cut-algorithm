library(modeltools)
library(flexclust)
library(compiler)
library(truncnorm)
library(mvtnorm)
library(foreach)
library(iterators)
library(doParallel)
library(dplyr)
library(ggplot2)
library(tidyr)
library(data.table)

#indicate if you are using a C code to write density function
cpp_yes<-FALSE
#if cpp_yes == TRUE, name the cpp package here:
cpp_package<-'Plummer2015Example'
#indicate your system
is_Unix<-TRUE

#slurm array task id specification
task_id_string <- Sys.getenv("SLURM_ARRAY_TASK_ID")
task_id <- as.numeric(task_id_string) 

#specification of parallel computing
core_num<-detectCores()-1

if(is_Unix){
  cl<-makeCluster(core_num,type = 'FORK')  
  registerDoParallel(cl)
}else{
  #cl<-makeCluster(detectCores( ))
  cl<-makeCluster(core_num)  
  registerDoParallel(cl)
  dmvnorm<-dmvnorm
}


#sd scaler for proposal distribution of \theta in the internal chain
sd_scaler<-c(1,10,20)[task_id]

#function to select 1st-nth minimal number with their index.
min.n<-function(x,n,value=TRUE){
  if(value==TRUE){x[order(x)][n]} else {order(x)[n]}
} 

Min.n<-function(x,n){
  ou<-rep(0,n)
  for(i in 1:n){
    ou[i]<-min.n(x,i,value = F)
  }
  return(ou)
}

#specification of data
Z <- c(7, 6, 10, 10, 1, 1, 10, 4, 35, 0, 10, 8, 4)
Npart <- c(111, 71, 162, 188, 145, 215, 166, 37, 173,143, 229, 696, 93)
Y <-  c(16, 215, 362, 97, 76, 62, 710, 56, 133,28, 62, 413, 194)
Npop <- c(26983, 250930, 829348, 157775, 150467, 352445, 553066, 26751, 75815, 150302, 354993, 3683043, 507218)



d_x<-2       # dimension of theta
d_y<-13       # dimension of phi

#loading the density function py and px
if(cpp_yes){
  library(cpp_package,character.only=TRUE)
}else{
  #density 1 Z|phi \times phi
  px<-function(phi,Z){
    PbZ<-rbind(phi,Z,Npart)
    out<-sum(apply(PbZ,FUN=function(y){dbeta(x=y[1],shape1 = (1+y[2]), shape2 = (1+y[3]-y[2]), log = T)},MARGIN = 2))
    return(out)
  }
  #unnormalizing density2 Y|theta&phi \times theta
  py<-function(Y,theta,phi){
    PbY<-rbind(phi,Y,Npop)
    out<- log(1/200)+sum(apply(PbY,FUN=function(y){dpois(x=y[2],lambda = (y[3]*0.001*exp(theta[1]+theta[2]*y[1])),log = T)},MARGIN = 2))
    return(out)
  }
}

#proposal 1
prox<-function(phi_n,phi){
  out<-sum(log(dtruncnorm(phi_n, a=0, b=1, mean = phi, sd = 0.005))) 
  return(out) 
}




#proposal 1 samling function
rprox<-function(phi){
  out<-rtruncnorm(1, a=0, b=1, mean = phi, sd = 0.005)
  return(out)
}


###################################################
#Do you need to build a new auxiliary \Phi_0 set?
is_newPhi0<-TRUE
if(is_newPhi0){
  #Choice of Auxiliary Parameter Set chain
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
  Phistar<-MAux(init,Z,num_run=15000,burn_in=10000,thin=1)
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
    Phistan<-as.matrix(Phistar)
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
  Phistar<-as.matrix(Phistar)
  PhiC<-Phistar[MMP(Phistan,num_sel=500),]
  
  #store PhiC for future use    
  write.csv(PhiC,"AuxPhi.csv",row.names = F)
}
###################################################


ncol.print <- function(dat) matrix(as.matrix(dat),ncol=ncol(dat),dimnames=NULL)
###################################################
#Load and select auxiliary \Phi_0 set from existing file?
#make sure your have removed all column name and row name in your file.csv!
load_newPhi0<-TRUE
if(load_newPhi0){
  Phistar<-ncol.print(as.matrix(fread("PhiC.csv",head=F)))
  
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
    Phistan<-as.matrix(Phistar)
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
  Phistar<-as.matrix(Phistar)
  PhiC<-Phistar[MMP(Phistan,num_sel=80),]
}
###################################################


###################################################
#directly load auxiliary \Phi_0 from existing file without further selection
#make sure your have removed all column name and row name in your file.csv!
load.NoSelect<-TRUE
if(load.NoSelect){
  PhiC<-ncol.print(as.matrix(fread("PhiC.csv",head=F)))
}
###################################################




###################################################
#Only use this when you want to constrain you \Phi_0 around its median.
constrain_Phi<-FALSE
if(constrain_Phi){
  #quantile used:
  Phi_qt<-0.25
  
  library(matrixStats)
  me_PhiC<-colMedians(PhiC)
  di_PhiC<-dist2(PhiC,me_PhiC)
  chs_PhiC<-rep(1,1000)
  for(i in 1:1000){
    chs_PhiC[i]<-sign(di_PhiC[i]>quantile(di_PhiC,Phi_qt) & di_PhiC[i]<quantile(di_PhiC,(1-Phi_qt)))
  }
  PhiC<-PhiC[chs_PhiC==1,]
}
###################################################







#test mod of density2 Y|theta&phi \times theta
#library(plotly)
# a<-rep(0,1000)
# b<-rep(0,1000)
# c<-rep(0,1000)
# for(k in 1:dim(PhiC)[1]){
#   xxx<-seq(-2.5,-1,0.05)
#   yyy<-seq(8,20,0.5)
#   PPP<-PhiC[k,]
#   PP<-rep(0,length(xxx)*length(yyy))
#   dim(PP)<-c(length(xxx),length(yyy))
#   for(i in 1:length(xxx)){
#     for(j in 1:length(yyy)){
#       PP[i,j]<-py(Y,c(xxx[i],yyy[j]),PPP)
#     }
#   }
#   a[k]<-xxx[which(PP == max(PP), arr.ind = TRUE)[1]]
#   b[k]<-yyy[which(PP == max(PP), arr.ind = TRUE)[2]]
#   c[k]<-max(PP)
#   print(k)
# }

# PP[which(PP == max(PP), arr.ind = TRUE)[1],which(PP == max(PP), arr.ind = TRUE)[2]]<-PP[which(PP == max(PP), arr.ind = TRUE)[1],which(PP == max(PP), arr.ind = TRUE)[2]]+500
# plot_ly(
#   x = yyy, y = xxx,
#   z = PP, type = "heatmap"
# )

PhiC<-as.matrix(PhiC)

#renew proposal distribution for \phi
#proposal 1
prox<-function(phi_n,phi){
  out<-sum(log(dtruncnorm(phi_n, a=0, b=1, mean = phi, sd = apply(PhiC,MARGIN = 2,sd))/10)) 
  return(out) 
}




#proposal 1 samling function
rprox<-function(phi){
  out<-rtruncnorm(1, a=0, b=1, mean = phi, sd = apply(PhiC,MARGIN = 2,sd)/10)
  return(out)
}


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


#proposal sampling for \t (which is \theta in the internal chain) and \phi
#write your own proposal distribution for \t
tranT<-function(t){
  re<-rmvnorm(1,mean = t,sigma = matrix(c(0.1648181,-0.3979341,-0.3979341,2.737874),ncol=2,nrow=2)/sd_scaler) %>% t()
  return(re)
}

pro_tp<-function(P,pro_tp_inp){
  I<-pro_tp_inp$I
  I_n<-rep(0,length(I))
  t<-pro_tp_inp$t
  t_n<-rep(1,d_x)
  phi<-pro_tp_inp$phi
  ran<-runif(1,0,1)
  if(ran<P){
    coin<-c(1,0)
    t_n<-tranT(t)
    out<-list(t=t_n,phi=phi,I=I,coin=coin)
  }else{
    coin<-c(0,1)
    I_n<-sample(seq(1,(dim(PhiC)[1]-1)),1)
    if(I_n<I){
      I_n<-I_n
    }else{
      I_n<-I_n+1
    }
    phi_n<-PhiC[I_n,]
    out<-list(t=t,phi=phi_n,I=I_n,coin=coin)
  }
  return(out)
}

#proposal density for t and phi
#write your own proposal density for \t
denT<-function(t_n,t,P){
  re<-dmvnorm(as.numeric(t_n),mean = as.numeric(t),sigma = matrix(c(0.1648181,-0.3979341,-0.3979341,2.737874),ncol=2,nrow=2)/sd_scaler)*P 
  return(re)
}

dpro_tp<-function(P,x_n,x){
  t<-x$t
  t_n<-x_n$t
  phi<-x$phi
  phi_n<-x_n$phi
  I<-x$I
  I_n<-x_n$I
  if(I_n>I){
    I_n<-I_n-1
  }else{
    I_n<-I_n
  }
  if(identical(t,t_n)){
    out<-1/(dim(PhiC)[1]-1)*(1-P)
  }else{
    #out<-dtruncnorm(t_n[1],a=-100, b=100,mean = t[1],sd=0.005)*dtruncnorm(t_n[2],a=-1000, b=1000,mean = t[2],sd=0.1)*P        #proposal for t may be changed
    out<-denT(t_n,t,P)
  }
  return(out)
}


#internal chain
no<-20000     #no is n_0
#indicate the accelerate parameter for speedy convergence of the normalizing constant W
acce_pa<-10

MA_in<-function(H,W,n,pai,n_trun){
  t<-H$t
  I<-H$I
  P<-0.75
  inp<-list(t=t,phi=PhiC[I,],I=I)
  ou<-pro_tp(P,inp)
  rate<-p_inv(ou$t,ou$phi,PhiC,W)+log(dpro_tp(P,inp,ou))-p_inv(t,PhiC[I,],PhiC,W)-log(dpro_tp(P,ou,inp))
  alfa<-min(1,exp(rate))
  rpan<-runif(1)
  t_n<-ou$t*sign(rpan<=alfa)+inp$t*sign(rpan>alfa)
  phi_n<-ou$phi*sign(rpan<=alfa)+inp$phi*sign(rpan>alfa)
  coin<-ou$coin*sign(rpan<=alfa)+c(0,0)
  iden<-function(x){
    Ou<-identical(x,phi_n)
    return(Ou)
  }
  I_n<-which(apply(PhiC,FUN = iden,MARGIN = 1))
  en<-rep(0,dim(PhiC)[1])
  en[I_n]<-1
  W_n<-exp(log(W)+acce_pa*(no/max(no,n))*(en-pai))        #no is n_0
  if(any(W_n>(10^(250+n_trun)))){
    W_n<-rep(1,length(W))            #truncated here
    n_trun<-n_trun+1
  }
  out<-list(t=t_n,I=I_n,W=W_n,n_trun=n_trun,coin=coin)
  return(out)
}

MA_in<-cmpfun(MA_in)

#set precision parameter to round \theta for high efficiency
#can be set to large number if you are confident in the performance of parallel computing
sig_dig<-c(3,2)   
#Digset<-c(3,4,10)
#if(d_x==1){
#  sig_dig<-Digset[task_id]
#}else{
#  sig_dig<-rep(0,d_x)
#  for(s in 1:d_x){
#    sig_dig[s]<-Digset[task_id,s]       #multiple setting of presicion parameter.
#  }
#}

#GRset<-seq(1,5.5,0.5)          #initial value set for separate chain

#initial value for normalizing constant W
www<-rep(1,dim(PhiC)[1])

#trimmed mean parameter, we have given up this approach, so always keep it as 1
uuu<-1

#construct sufficient auxiliary set
init<-list(theta=c(-2,13),phi=apply(PhiC,MARGIN = 2,median),t=as.matrix(c(-2,13)),I=1)          #t should be inversed vector.
MA_aux<-function(init,Z,Y,PhiC,num_run=1000,burn_in=500){
  theta<-init$theta
  theta_n<-rep(0,d_x)
  phi<-init$phi
  n_trun<-1
  H<-list(t=init$t,I=init$I)
  W<-www
  ColH<-list(t=init$t,I=init$I,w=W[1])      #Comulative information
  pai<-rep(1/length(W),length(W))          # sampling frequency
  Count_Tt<-1
  st.W<-rep(0,length(W))
  st.I<-rep(0,num_run)
  sto.autheta.i<-as.matrix(rep(theta,(num_run-burn_in))) #note the dimension of theta
  dim(sto.autheta.i)<-c((num_run-burn_in),length(theta))
  coin<-c(0,0)
  foreach(i = icount(num_run)) %do% {
    if(i<=burn_in){
      for(k in 1:1){
        InR<-MA_in(H,W,n=((i-1)*1+k),pai,n_trun)
        W<-InR$W
        H$t<-InR$t
        H$I<-InR$I
        n_trun<-InR$n_trun 
      }
      coin<-coin+as.numeric(InR$coin)
    }else{
      j<-i-burn_in
      for(k in 1:1){
        InR<-MA_in(H,W,n=((i-1)*1+k),pai,n_trun)
        W<-InR$W
        H$t<-InR$t
        H$I<-InR$I
        n_trun<-InR$n_trun 
      }
      coin<-coin+as.numeric(InR$coin)
      
      ColH$t<-round(InR$t,digits = sig_dig)
      ColH$I<-InR$I
      ColH$w<-W[InR$I]
      
      sto.autheta.i[(i-burn_in),]<-ColH$t
      
      #calculate sampling frequency
      bas<-rbind(ColH$I,ColH$w,ColH$t)
      if(j==1){
        Tt<-as.matrix(ColH$t)
      }else{
        Tt<-as.matrix(unique(cbind(Tt,ColH$t),MARGIN=2))
      }
      
      phi_n<-rprox(phi)
      if(j==1){
        log.deno<-log(bas[2])+py(Y,ColH$t,phi_n)-py(Y,ColH$t,PhiC[bas[1],])
        log.numr<-log(bas[2])+py(Y,ColH$t,phi_n)-py(Y,ColH$t,PhiC[bas[1],])
        log.fenzi_o<-py(Y,Tt,phi_n)
        Ptau<-1
        rpt<-1
      }else{
        log.fenzi<-function(t){  
          return(max(py(Y,t,phi_n),-4000))
        }
        
        if(Count_Tt==dim(Tt)[2]){
          iidentical<-function(x){
            return(identical(x,as.numeric(ColH$t)))
          }
          if(dim(Tt)[2]>1){
            log.fenzi_n<-apply(Tt,FUN = log.fenzi,MARGIN = 2)  
          }else{
            log.fenzi_n<-py(Y,Tt,phi_n)
          }
          
          nchan<-which(apply(Tt,FUN=iidentical,MARGIN = 2))
          rpt[nchan]<-rpt[nchan]+1
          log.numr<-log.numr-log.fenzi_o+log.fenzi_n
          log.nchanadd<-log(bas[2])+py(Y,ColH$t,phi_n)-py(Y,ColH$t,PhiC[bas[1],])
          if(exp(log.numr[nchan])+exp(log.nchanadd)==0){
            log.numr[nchan]<-max(log.numr[nchan],log.nchanadd)
          }else{
            log.numr[nchan]<-log(exp(log.numr[nchan])+exp(log.nchanadd))
          }
          log.fenzi_o<-log.fenzi_n
          
          ###########
          if(TRUE){
            log.numr.ii<-log.numr
          }else{
            log.numr.i<-log.numr-log(rpt)
            aaa<-quantile(log.numr.i,uuu)
            bbb<-min(log.numr.i)
            ccut<-function(x){
              if(x>aaa){
                x<-bbb
              }
              return(x)
            }
            log.numr.ii<-sapply(log.numr.i, ccut)+log(rpt)
          }
          ##########
          
          if(max(log.numr.ii)<300){                                       #scale the density
            log.numr.scale<-log.numr.ii+(300-max(log.numr.ii))
          }else{
            log.numr.scale<-log.numr.ii
          }
          #print(c(max(log.numr.scale),min(log.numr.scale),py(Y,ColH$t,phi_n),min(log(W))))
          deno<-sum(exp(log.numr.scale))
          Ptau<-exp(log.numr.scale)/deno
        }else{
          rpt<-append(rpt,1)
          log.fenzi_n<-apply(Tt,FUN = log.fenzi,MARGIN = 2)  
          nchanadd<-length(log.fenzi_n)
          log.numr<-log.numr-log.fenzi_o+log.fenzi_n[-nchanadd]
          log.numr[nchanadd]<-log(bas[2])+py(Y,ColH$t,phi_n)-py(Y,ColH$t,PhiC[bas[1],])
          log.fenzi_o<-log.fenzi_n
          ###########
          if(TRUE){
            log.numr.ii<-log.numr
          }else{
            log.numr.i<-log.numr-log(rpt)
            aaa<-quantile(log.numr.i,uuu)
            bbb<-min(log.numr.i)
            ccut<-function(x){
              if(x>aaa){
                x<-bbb
              }
              return(x)
            }
            log.numr.ii<-sapply(log.numr.i, ccut)+log(rpt)
          }
          ##########
          
          if(max(log.numr.ii)<300){                                       #scale the density
            log.numr.scale<-log.numr.ii+(300-max(log.numr.ii))
          }else{
            log.numr.scale<-log.numr.ii
          }
          #print(c(max(log.numr.scale),min(log.numr.scale),py(Y,ColH$t,phi_n),min(log(W))))
          deno<-sum(exp(log.numr.scale))
          Ptau<-exp(log.numr.scale)/deno
        }
      }
      Count_Tt<-dim(Tt)[2]
      
      tau_n<-Tt[,sample(x=seq(1,dim(Tt)[2]),size = 1,prob = Ptau)]
      if(length(tau_n)==1){
        theta_n<-runif(1,min=tau_n-5*10^(-sig_dig-1),max=tau_n+5*10^(-sig_dig-1))  #uniformly drawing \theta
      }else{
        for(q in 1:length(tau_n)){
          theta_n[q]<-runif(1, min=tau_n[q]-5*10^(-sig_dig[q]-1),max=tau_n[q]+5*10^(-sig_dig[q]-1))
        }
      }
      rate<-px(phi_n,Z)+prox(phi,phi_n)-px(phi,Z)-prox(phi_n,phi)
      alfa<-min(1,exp(rate))
      rpan<-runif(1)
      phi<-phi_n*sign(rpan<=alfa)+phi*sign(rpan>alfa)
      theta<-theta_n*sign(rpan<=alfa)+theta*sign(rpan>alfa)
    }
    st.W<-log(W)
    st.I[i]<-InR$I
    ac_pro<-coin/i
    if(i %in% seq(1,num_run,100)){
      print(c(i,InR$t,InR$I,ac_pro))
    }
  }
  MA_aux_out<-list(Tt=Tt,autheta=sto.autheta.i,Ptau=Ptau,rpt=rpt,log.fenzi_o=log.fenzi_o,log.numr=log.numr,t=H$t,I=H$I,n_trun=n_trun,aux_num_run=num_run,W=W,theta=theta,phi=phi,coin=coin)
  return(MA_aux_out)
}

MA_aux<-cmpfun(MA_aux)

MA_aux_out<-MA_aux(init,Z,Y,PhiC,num_run = 1501000,burn_in=1500000)

print(paste('internal chain finished at',Sys.time()))

#store the normalizing constant W
write.csv(log(MA_aux_out$W),paste("logW",task_id,".csv",sep=""))

#store the record of auxiliary \phi and \theta
RR.in<-data.table(auphi=MA_aux_out$auphi,autheta=MA_aux_out$autheta)
write.csv(RR.in,paste("Result_aux",task_id,".csv",sep = ""))










#external chain
init<-list(theta=MA_aux_out$theta,phi=MA_aux_out$phi,t=MA_aux_out$t,I=MA_aux_out$I)          #t should be inversed vector.
MA_ex<-function(Aux_Tt,init,Z,Y,PhiC,num_run=1000,burn_in=500,thin=1){
  rpt<-MA_aux_out$rpt
  theta<-init$theta
  theta_n<-rep(0,d_x)
  phi<-init$phi
  n_trun<-MA_aux_out$n_trun
  H<-list(t=init$t,I=init$I)
  W<-Aux_Tt$W
  coin<-Aux_Tt$coin
  ColH<-list(t=init$t,I=init$I,w=W[1])      #Comulative information
  ColH<-ColH
  pai<-rep(1/length(W),length(W))          # sampling frequency
  sto.phi<-as.matrix(rep(phi,((num_run-burn_in)/thin))) #note the dimension of phi
  dim(sto.phi)<-c(((num_run-burn_in)/thin),length(phi))
  sto.theta<-as.matrix(rep(theta,((num_run-burn_in)/thin))) #note the dimension of theta
  dim(sto.theta)<-c(((num_run-burn_in)/thin),length(theta))
  sto.auphi<-rep(0,((num_run-burn_in)/thin)) #note the dimension of phi
  sto.autheta<-as.matrix(rep(theta,((num_run-burn_in)/thin))) #note the dimension of theta
  dim(sto.autheta)<-c(((num_run-burn_in)/thin),length(theta))
  sto.time<-rep(0,((num_run-burn_in)/thin))
  Tt<-Aux_Tt$Tt
  Ptau<-Aux_Tt$Ptau
  log.fenzi_o<-Aux_Tt$log.fenzi_o
  log.numr<-Aux_Tt$log.numr
  InRadd<-Aux_Tt$aux_num_run
  Count_Tt<-dim(MA_aux_out$Tt)[2]
  if(!is_Unix){
    clusterExport(cl, list('py','Y','Npop','dmvnorm'), envir=environment())
  }
  foreach(i = icount(num_run)) %do%{
    if((i<=burn_in)|(i%%thin!=0)){
      for(k in 1:1){
        InR<-MA_in(H,W,n=(((i-1)*1+k)+InRadd),pai,n_trun)
        W<-InR$W
        H$t<-InR$t
        H$I<-InR$I
        n_trun<-InR$n_trun 
      }
      coin<-coin+as.numeric(InR$coin)
      
      ColH$t<-round(InR$t,digits = sig_dig)
      ColH$I<-InR$I
      ColH$w<-W[InR$I]
      #calculate sampling frequency
      bas<-rbind(ColH$I,ColH$w,ColH$t)
      Tt<-as.matrix(unique(cbind(Tt,ColH$t),MARGIN=2))
      phi_n<-rprox(phi)
      
      log.fenzi<-function(t){   
        return(max(py(Y,t,phi_n),-4000))
      }
      
      if(Count_Tt==dim(Tt)[2]){
        iidentical<-function(x){
          return(identical(x,as.numeric(ColH$t)))
        }
        if(dim(Tt)[2]>1){
          if(is_Unix){
            log.fenzi_n<-as.numeric(mclapply(lapply(seq_len(ncol(Tt)), function(i) Tt[, i]),FUN = log.fenzi))
          }else{
            clusterExport(cl, list('phi_n'), envir=environment())
            log.fenzi_n<-parApply(cl,Tt,FUN = log.fenzi,MARGIN = 2)
          }
        }else{
          log.fenzi_n<-py(Y,Tt,phi_n)
        }
        
        if(is_Unix){
          nchan<-which(unlist(mclapply(lapply(seq_len(ncol(Tt)), function(i) Tt[, i]),FUN = iidentical)))
        }else{
          clusterExport(cl, list('identical','ColH'),envir=environment())
          nchan<-which(parApply(cl,Tt,FUN=iidentical,MARGIN = 2))
        }
        rpt[nchan]<-rpt[nchan]+1
        log.numr<-log.numr-log.fenzi_o+log.fenzi_n
        log.nchanadd<-log(bas[2])+py(Y,ColH$t,phi_n)-py(Y,ColH$t,PhiC[bas[1],])
        nchanadd<-bas[2]*exp(py(Y,ColH$t,phi_n)-py(Y,ColH$t,PhiC[bas[1],]))
        if(exp(log.numr[nchan])+exp(log.nchanadd)==0){
          log.numr[nchan]<-max(log.numr[nchan],log.nchanadd)
        }else{
          log.numr[nchan]<-log(exp(log.numr[nchan])+exp(log.nchanadd))
        }
        log.fenzi_o<-log.fenzi_n
        ###########
        if(TRUE){
          log.numr.ii<-log.numr
        }else{
          log.numr.i<-log.numr-log(rpt)
          aaa<-quantile(log.numr.i,uuu)
          bbb<-min(log.numr.i)
          ccut<-function(x){
            if(x>aaa){
              x<-bbb
            }
            return(x)
          }
          log.numr.ii<-sapply(log.numr.i, ccut)+log(rpt)
        }
        ##########
        
        if(max(log.numr.ii)<300){                                       #scale the density
          log.numr.scale<-log.numr.ii+(300-max(log.numr.ii))
        }else{
          log.numr.scale<-log.numr.ii
        }
        deno<-sum(exp(log.numr.scale))
        Ptau<-exp(log.numr.scale)/deno
      }else{
        rpt<-append(rpt,1)
        if(is_Unix){
          log.fenzi_n<-as.numeric(mclapply(lapply(seq_len(ncol(Tt)), function(i) Tt[, i]),FUN = log.fenzi))
        }else{
          clusterExport(cl, list('phi_n'), envir=environment())
          log.fenzi_n<-parApply(cl,Tt,FUN = log.fenzi,MARGIN = 2)
        }
        nchanadd<-length(log.fenzi_n)
        log.numr<-log.numr-log.fenzi_o+log.fenzi_n[-nchanadd]
        log.numr[nchanadd]<-log(bas[2])+py(Y,ColH$t,phi_n)-py(Y,ColH$t,PhiC[bas[1],])
        log.fenzi_o<-log.fenzi_n
        ###########
        if(TRUE){
          log.numr.ii<-log.numr
        }else{
          log.numr.i<-log.numr-log(rpt)
          aaa<-quantile(log.numr.i,uuu)
          bbb<-min(log.numr.i)
          ccut<-function(x){
            if(x>aaa){
              x<-bbb
            }
            return(x)
          }
          log.numr.ii<-sapply(log.numr.i, ccut)+log(rpt)
        }
        ##########
        
        if(max(log.numr.ii)<300){                                       #scale the density
          log.numr.scale<-log.numr.ii+(300-max(log.numr.ii))
        }else{
          log.numr.scale<-log.numr.ii
        }
        deno<-sum(exp(log.numr.scale))
        Ptau<-exp(log.numr.scale)/deno
      }
      Count_Tt<-dim(Tt)[2]
      
      tau_n<-Tt[,sample(x=seq(1,dim(Tt)[2]),size = 1,prob = Ptau)]
      if(length(tau_n)==1){
        theta_n<-runif(1,min=tau_n-5*10^(-sig_dig-1),max=tau_n+5*10^(-sig_dig-1))  #uniformly drawing \theta
      }else{
        for(q in 1:length(tau_n)){
          theta_n[q]<-runif(1,min=tau_n[q]-5*10^(-sig_dig[q]-1),max=tau_n[q]+5*10^(-sig_dig[q]-1))
        }
      }
      rate<-px(phi_n,Z)+prox(phi,phi_n)-px(phi,Z)-prox(phi_n,phi)
      alfa<-min(1,exp(rate))
      rpan<-runif(1)
      phi<-phi_n*sign(rpan<=alfa)+phi*sign(rpan>alfa)
      theta<-theta_n*sign(rpan<=alfa)+theta*sign(rpan>alfa)
    }else{
      earlier_time<-Sys.time()
      for(k in 1:1){
        InR<-MA_in(H,W,n=(((i-1)*1+k)+InRadd),pai,n_trun)
        W<-InR$W
        H$t<-InR$t
        H$I<-InR$I
        n_trun<-InR$n_trun 
      }
      coin<-coin+as.numeric(InR$coin)
      
      ColH$t<-round(InR$t,digits = sig_dig)
      ColH$I<-InR$I
      ColH$w<-W[InR$I]
      
      sto.auphi[((i-burn_in)/thin)]<-InR$I
      sto.autheta[((i-burn_in)/thin),]<-ColH$t
      
      #calculate sampling frequency
      bas<-rbind(ColH$I,ColH$w,ColH$t)
      Tt<-as.matrix(unique(cbind(Tt,ColH$t),MARGIN=2))
      phi_n<-rprox(phi)
      log.fenzi<-function(t){   
        return(max(py(Y,t,phi_n),-4000))
      }
      
      if(Count_Tt==dim(Tt)[2]){
        iidentical<-function(x){
          return(identical(x,as.numeric(ColH$t)))
        }
        if(dim(Tt)[2]>1){
          if(is_Unix){
            log.fenzi_n<-as.numeric(mclapply(lapply(seq_len(ncol(Tt)), function(i) Tt[, i]),FUN = log.fenzi))
          }else{
            clusterExport(cl, list('phi_n'), envir=environment())
            log.fenzi_n<-parApply(cl,Tt,FUN = log.fenzi,MARGIN = 2)
          }
        }else{
          log.fenzi_n<-py(Y,Tt,phi_n)
        }
        
        if(is_Unix){
          nchan<-which(unlist(mclapply(lapply(seq_len(ncol(Tt)), function(i) Tt[, i]),FUN = iidentical)))
        }else{
          clusterExport(cl, list('identical','ColH'),envir=environment())
          nchan<-which(parApply(cl,Tt,FUN=iidentical,MARGIN = 2))
        }
        rpt[nchan]<-rpt[nchan]+1
        log.numr<-log.numr-log.fenzi_o+log.fenzi_n
        log.nchanadd<-log(bas[2])+py(Y,ColH$t,phi_n)-py(Y,ColH$t,PhiC[bas[1],])
        nchanadd<-bas[2]*exp(py(Y,ColH$t,phi_n)-py(Y,ColH$t,PhiC[bas[1],]))
        if(exp(log.numr[nchan])+exp(log.nchanadd)==0){
          log.numr[nchan]<-max(log.numr[nchan],log.nchanadd)
        }else{
          log.numr[nchan]<-log(exp(log.numr[nchan])+exp(log.nchanadd))
        }
        log.fenzi_o<-log.fenzi_n
        ###########
        if(TRUE){
          log.numr.ii<-log.numr
        }else{
          log.numr.i<-log.numr-log(rpt)
          aaa<-quantile(log.numr.i,uuu)
          bbb<-min(log.numr.i)
          ccut<-function(x){
            if(x>aaa){
              x<-bbb
            }
            return(x)
          }
          log.numr.ii<-sapply(log.numr.i, ccut)+log(rpt)
        }
        ##########
        
        if(max(log.numr.ii)<300){                                       #scale the density
          log.numr.scale<-log.numr.ii+(300-max(log.numr.ii))
        }else{
          log.numr.scale<-log.numr.ii
        }
        
        deno<-sum(exp(log.numr.scale))
        Ptau<-exp(log.numr.scale)/deno
      }else{
        rpt<-append(rpt,1)
        if(is_Unix){
          log.fenzi_n<-as.numeric(mclapply(lapply(seq_len(ncol(Tt)), function(i) Tt[, i]),FUN = log.fenzi))
        }else{
          clusterExport(cl, list('phi_n'), envir=environment())
          log.fenzi_n<-parApply(cl,Tt,FUN = log.fenzi,MARGIN = 2)
        }
        nchanadd<-length(log.fenzi_n)
        log.numr<-log.numr-log.fenzi_o+log.fenzi_n[-nchanadd]
        log.numr[nchanadd]<-log(bas[2])+py(Y,ColH$t,phi_n)-py(Y,ColH$t,PhiC[bas[1],])
        log.fenzi_o<-log.fenzi_n
        ###########
        if(TRUE){
          log.numr.ii<-log.numr
        }else{
          log.numr.i<-log.numr-log(rpt)
          aaa<-quantile(log.numr.i,uuu)
          bbb<-min(log.numr.i)
          ccut<-function(x){
            if(x>aaa){
              x<-bbb
            }
            return(x)
          }
          log.numr.ii<-sapply(log.numr.i, ccut)+log(rpt)
        }
        ##########
        
        if(max(log.numr.ii)<300){                                       #scale the density
          log.numr.scale<-log.numr.ii+(300-max(log.numr.ii))
        }else{
          log.numr.scale<-log.numr.ii
        }
        deno<-sum(exp(log.numr.scale))
        Ptau<-exp(log.numr.scale)/deno
      }
      Count_Tt<-dim(Tt)[2]
      
      tau_n<-Tt[,sample(x=seq(1,dim(Tt)[2]),size = 1,prob = Ptau)]
      if(length(tau_n)==1){
        theta_n<-runif(1,min=tau_n-5*10^(-sig_dig-1),max=tau_n+5*10^(-sig_dig-1))  #uniformly drawing \theta
      }else{
        for(q in 1:length(tau_n)){
          theta_n[q]<-runif(1,min=tau_n[q]-5*10^(-sig_dig[q]-1),max=tau_n[q]+5*10^(-sig_dig[q]-1))
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
    
    debug <- F   #Debug
    if(debug=T){
    if(i %in% seq(1000,500000,1000)){
      Tem_out<-list(phi=sto.phi[1:((i-burn_in)/thin),],auxphi=sto.auphi[1:((i-burn_in)/thin)],theta=sto.theta[1:((i-burn_in)/thin),],auxtheta=sto.autheta[1:((i-burn_in)/thin),],time=sto.time[1:((i-burn_in)/thin)],Tt=Tt[2,],log.Ptau_fenmu=log.numr-log.fenzi_o)
      Tem_RR<-data.table(phi= Tem_out$phi,auxphi= Tem_out$auxphi,theta= Tem_out$theta,auxtheta= Tem_out$auxtheta,time= Tem_out$time)
      write.csv(Tem_RR,paste("Result",task_id,".csv",sep = ""))
      Tem_Pr<-data.table(Tt=Tem_out$Tt,Ptau=Tem_out$log.Ptau_fenmu)
      write.csv(Tem_Pr,paste("TauProbability",task_id,".csv",sep = ""))
      pdf(paste("theta_plot",task_id,".pdf",sep=""))
      fig<-ggplot(Tem_RR,aes(x=theta.V2,y=theta.V1))+geom_point(col='lightblue')
      print(fig)
      dev.off()
      print('Record result')
    }
    }
    ac_pro<-coin/(i+InRadd)
    
    print(c(i,theta,n_trun,ac_pro))
    #if(sign(rpan<=alfa)==1){
    #  xx<-seq(-2.5,-1,0.005)
    #  yy<-seq(8,23,0.05)
    #  PP<-rep(0,length(xx)*length(yy))
    #  dim(PP)<-c(length(xx),length(yy))
    #  for(s in 1:dim(Tt)[2]){
    #    PP[which.min(abs(xx-Tt[1,][s])),which.min(abs(yy-Tt[2,][s]))]<-PP[which.min(abs(xx-Tt[1,][s])),which.min(abs(yy-Tt[2,][s]))]+Ptau[s]
    #  }
    #  heatmap(PP,Rowv = NA,Colv = NA,scale = 'none')
    #}
  }
  OUT<-list(phi=sto.phi,theta=sto.theta,auphi=sto.auphi,autheta=sto.autheta,Tt=Tt[1,],Ptau=Ptau,time=sto.time,log.Ptau_fenmu=log.numr-log.fenzi_o)
  return(OUT)
}

MA_ex<-cmpfun(MA_ex)

Result<-MA_ex(MA_aux_out,init,Z,Y,PhiC,num_run=500000,burn_in=0,thin=1)

stopCluster(cl)

RR<-data.table(phi=Result$phi,auphi=Result$auphi,theta=Result$theta,autheta=Result$autheta,time=Result$time)
TauPro<-data.table(Tt=Result$Tt,Ptau=Result$log.Ptau_fenmu)
write.csv(RR,paste("Result",task_id,".csv",sep = ""))
write.csv(TauPro,paste("TauProbability",task_id,".csv",sep = ""))
