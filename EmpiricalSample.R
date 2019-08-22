library(truncnorm)
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

d_x<-1       # dimension of theta
d_y<-N       # dimension of phi

#density 1 Z|phi \times phi
px<-function(phi,Z){
  PbZ<-rbind(phi,Z)
  out<-sum(apply(PbZ,FUN=function(x){log(((x[1])^(-(n+1)/2))*exp(-x[2]/(2*x[1])))},MARGIN = 2))
  return(out)
}
#unnormalizing density2 Y|theta&phi \times theta
py<-function(Y,theta,phi){
  PbY<-rbind(phi,Y)
  out<-log(1/(theta+mean(phi)/n))+sum(apply(PbY,FUN=function(x){log((1/sqrt(theta+x[1]/n))*exp(-(x[2]^2)/(2*(theta+x[1]/n))))},MARGIN = 2))
  return(out)
}

#proposal 1
prox<-function(phi_n,phi){
  out<-sum(log(dtruncnorm(phi_n, a=0, b=100, mean = phi, sd = 0.1))) 
  return(out) 
}




#proposal 1 samling function
rprox<-function(phi){
  out<-rtruncnorm(1, a=0, b=100, mean = phi, sd = 0.1)
  return(out)
}



x1<-read.csv(paste0('TauProbability_',task_id,".csv"))

init<-list(phi=rep(1,d_y))
CanPhi<-function(init,Z,num_run=1000,burn_in=500,thin=1){
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
    print(i)
  }
  out<-sto.phi
  return(out)
}
Phi1<-CanPhi(init,Z,num_run=150000,burn_in=50000,thin=10)

sam1<-rep(0,10000)
Tt<-as.matrix(x1$Tt)
for(i in 1:10000){
  fenzi<-function(t){   
    return(exp(py(Y,t,Phi1[i,])))
  }
  fenzi_n<-apply(Tt,FUN = fenzi,MARGIN = 1)  
  Ptau<-x1$Ptau*fenzi_n
  sam1[i]<-x1$Tt[sample(x=seq(1,length(x1$Tt)),size = 1,prob = Ptau)]+runif(1,-0.0005,0.0005)
  print(i)
}


write.csv(sam1,paste("EpiricalSample_",task_id,".csv",sep = ""))
