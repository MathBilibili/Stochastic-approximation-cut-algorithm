library(truncnorm)

task_id_string <- Sys.getenv("SLURM_ARRAY_TASK_ID")
task_id <- as.numeric(task_id_string) 

data<-read.csv("simu_data.csv")
y.m<-apply(data,1,mean)
N<-dim(data)[1]
n<-dim(data)[2]
s<-rep(0,length(y.m))
for(i in 1:N){
  for(j in 1:n){
    s[i]<-s[i]+(data[i,j]-y.m[i])^2
  }
}

post<-function(theta,phi,y.m,s){
  a<-log(1/(theta+mean(phi)/n))
  b<-0
  for(k in 1:N){
    b<-b+log(((phi[k])^(-(n+1)/2))*exp(-s[k]/(2*phi[k]))*(1/sqrt(theta+phi[k]/n))*exp(-(y.m[k]^2)/(2*(theta+phi[k]/n))))
  }
  out<-a+b
  return(out)
}

MC<-function(init,y.m,s,num_run=1000,burn_in=500,thin=1){
  theta<-init[1]
  phi<-init[2:length(init)]
  theta.n<-init[1]
  phi.n<-init[2:length(init)]
  sto.theta<-as.matrix(rep(0,((num_run-burn_in)/thin)))
  dim(sto.theta)<-c(((num_run-burn_in)/thin),length(theta))
  sto.phi<-as.matrix(rep(phi,((num_run-burn_in)/thin))) #note the dimension of phi
  dim(sto.phi)<-c(((num_run-burn_in)/thin),length(phi))
  for(i in 1:num_run){
    if((i<=burn_in)|(i%%thin!=0)){
      theta.n<-rtruncnorm(1, a=0, b=100, mean = theta, sd = 0.1)
      phi.n<-rtruncnorm(1, a=0, b=100, mean = phi, sd = 0.05)
      rate<-post(theta.n,phi.n,y.m,s)+log(prod(dtruncnorm(phi,a=0,b=100,mean=phi.n,sd=0.05)))+log(dtruncnorm(theta,a=0,b=100,mean=theta.n,sd=0.1))-post(theta,phi,y.m,s)-log(prod(dtruncnorm(phi.n,a=0,b=100,mean=phi,sd=0.05)))-log(dtruncnorm(theta.n,a=0,b=100,mean=theta,sd=0.1))
      alfa<-min(1,exp(rate))
      rpan<-runif(1)
      theta<-theta.n*sign(rpan<=alfa)+theta*sign(rpan>alfa)
      phi<-phi.n*sign(rpan<=alfa)+phi*sign(rpan>alfa)
    }else{
      theta.n<-rtruncnorm(1, a=0, b=100, mean = theta, sd = 0.1)
      phi.n<-rtruncnorm(1, a=0, b=100, mean = phi, sd = 0.05)
      rate<-post(theta.n,phi.n,y.m,s)+log(prod(dtruncnorm(phi,a=0,b=100,mean=phi.n,sd=0.05)))+log(dtruncnorm(theta,a=0,b=100,mean=theta.n,sd=0.1))-post(theta,phi,y.m,s)-log(prod(dtruncnorm(phi.n,a=0,b=100,mean=phi,sd=0.05)))-log(dtruncnorm(theta.n,a=0,b=100,mean=theta,sd=0.1))
      alfa<-min(1,exp(rate))
      rpan<-runif(1)
      theta<-theta.n*sign(rpan<=alfa)+theta*sign(rpan>alfa)
      phi<-phi.n*sign(rpan<=alfa)+phi*sign(rpan>alfa)
      sto.phi[((i-burn_in)/thin),]<-phi
      sto.theta[((i-burn_in)/thin),]<-theta
    }
    print(c(i,theta,phi[10],rate))
  }
  out<-list(phi=sto.phi[,10],theta=sto.theta)
  return(out)
}

init<-rep(5.5,N+1)
Result<-MC(init,y.m,s,num_run=110000,burn_in=10000,thin=10)
write.csv(Result,paste("Liu_mcmc",task_id,".csv",sep = ""))
