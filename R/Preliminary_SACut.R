Preliminary_SACut<-function(init=list(theta=c(-2,13),phi=rep(1,13),t=as.matrix(c(-2,13)),I=1),
                            PhiC,numrun=1000,auxrun=500,no=1000,acce_pa=1, sig_dig, CutModel){
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

  PhiC<-as.matrix(PhiC)

  #initial value for normalizing constant W
  www<-rep(1,dim(PhiC)[1])

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

  MA_aux_out<-MA_aux(init,Z,Y,PhiC,num_run = numrun,burn_in=auxrun)

  return(MA_aux_out)

}
