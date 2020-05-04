library(compiler)
library(truncnorm)
library(data.table)
#data specification
task_id_string <- Sys.getenv("SLURM_ARRAY_TASK_ID")
task_id <- as.numeric(task_id_string)



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

Z <- as.numeric(fread('Z.txt')[[1]])
ind_var<-seq(1,50)
#Y <-  rnorm(50,mean = (2+ind_var),sd=3)
Y <- as.numeric(fread('Y.txt')[[1]])



d_x<-1       # dimension of theta
d_y<-1       # dimension of phi

#density 1 Z|phi \times phi
px<-function(phi,Z){
  out<-sum(sapply(Z,FUN=function(y){dnorm(y,mean = phi,sd=1,log = T)}))+dnorm(phi,mean = 0,sd=100,log = T)
  return(out)
}
#unnormalizing density2 Y|theta&phi \times theta
py<-function(Y,theta,phi){
  PbY<-rbind(Y,ind_var)
  out<-sum(apply(PbY,FUN=function(y){dnorm(y[1],mean = theta+phi*y[2],sd = 3,log = T)},MARGIN = 2))+dnorm(theta,mean = 0,sd=100,log = T)
  return(out)
}

#proposal 1
prox<-function(phi_n,phi){
  out<-dnorm(phi_n, mean=phi, sd=0.25,log = T) 
  return(out) 
}




#proposal 1 samling function
rprox<-function(phi){
  out<-rnorm(1, mean=phi, sd=0.25) 
  return(out)
}




#proposal 2
proy<-function(theta_n,theta){
  out<-dnorm(theta_n,mean = theta,sd=0.00025,log = T)    #small sd for slow convergence.
  return(out) 
}

#proposal 1 samling function
rproy<-function(theta){
  out<-rnorm(1, mean = theta, sd = 0.00025)
  return(out)
}


MC1<-function(phi,Z){
  phi_n<-rprox(phi)
  rate<-px(phi_n,Z)+prox(phi,phi_n)-px(phi,Z)-prox(phi_n,phi)
  alfa<-min(1,exp(rate))
  rpan<-runif(1)
  out<-phi_n*sign(rpan<=alfa)+phi*sign(rpan>alfa)
  return(out)
}



MC2<-function(theta,phi,Y,N){              #N is number of internal iterations
  for(k in 1:N){
    theta_n<-rproy(theta)
    rate<-py(Y,theta_n,phi)+proy(theta,theta_n)-py(Y,theta,phi)-proy(theta_n,theta)
    alfa<-min(1,exp(rate))
    rpan<-runif(1)
    theta<-theta_n*sign(rpan<=alfa)+theta*sign(rpan>alfa)
  }
  out<-theta
  return(out)
}




init<-list(theta=1,phi=0.1)  
MC<-function(init,Y,Z,inter=100,num_run=1000,burn_in=500,thin=1){
  theta<-init$theta
  phi<-init$phi
  sto.phi<-as.matrix(rep(phi,((num_run-burn_in)/thin))) #note the dimension of phi
  dim(sto.phi)<-c(((num_run-burn_in)/thin),length(phi))
  sto.theta<-as.matrix(rep(theta,((num_run-burn_in)/thin))) #note the dimension of theta
  dim(sto.theta)<-c(((num_run-burn_in)/thin),length(theta))
  for(i in 1:num_run){
    if((i<=burn_in)|(i%%thin!=0)){
      phi_n<-MC1(phi,Z)
      theta_n<-MC2(theta,phi_n,Y,inter)
      phi<-phi_n
      theta<-theta_n
    }else{
      phi_n<-MC1(phi,Z)
      theta_n<-MC2(theta,phi_n,Y,inter)
      phi<-phi_n
      theta<-theta_n
      sto.phi[((i-burn_in)/thin),]<-phi
      sto.theta[((i-burn_in)/thin),]<-theta
    }
    print(c(i,theta))
  }
  out<-list(theta=sto.theta,phi=sto.phi)
  return(out)
}


inter_set<-c(rep(1,9),rep(100,9),rep(1000,9),rep(1500,9),rep(3000,9),rep(5000,9))
index2<-rep(c(1,2,3,4,5,6,7,8,9),6)
index<-c(rep(1,9),rep(2,9),rep(3,9),rep(4,9),rep(5,9),rep(6,9))

q1<-Sys.time()
#Result<-MC(init,Y,Z,inter=inter_set[task_id],num_run=50000,burn_in=20000,thin=1)
inter_num<-c(10,500,1000,1500,2000,3000)
Result<-MC(init,Y,Z,inter=inter_num[task_id],num_run=50000,burn_in=0,thin=1)
q2<-Sys.time()
print(paste("Martin_",inter_num[task_id],"_",index2[task_id],".csv",sep = ""))
print(q2-q1)
write.csv(Result,paste("Martin_",inter_num[task_id],"_",index2[task_id],".csv",sep = ""))


############## analyse
analys<-F
if(analys==T){
index2<-c(0,1,2,3,4,5,6,7,8,9)

d<-list()
m1<-rep(0,10)
m2<-rep(0,10)
for(i in index2){
  d[[i+1]]<-read.csv(paste("Martin_",2,"_",i,".csv",sep = ""))
  m1[i+1]<-mean(d[[i+1]]$theta.1)
  m2[i+1]<-mean(d[[i+1]]$theta.2)
}

mean(m1)
mean(m2)
}
