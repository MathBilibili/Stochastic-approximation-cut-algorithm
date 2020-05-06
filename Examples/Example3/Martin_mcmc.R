library(flexclust)
library(compiler)
library(truncnorm)
#data specification

#indicate if you are using a C code to write density function
cpp_yes<-TRUE
#if cpp_yes == TRUE, name the cpp package here:
cpp_package<-'Plummer2015Example'

#slurm array task id specification
task_id_string <- Sys.getenv("SLURM_ARRAY_TASK_ID")
task_id <- as.numeric(task_id_string) 

Z <- c(7, 6, 10, 10, 1, 1, 10, 4, 35, 0, 10, 8, 4)
Npart <- c(111, 71, 162, 188, 145, 215, 166, 37, 173,143, 229, 696, 93)
Y <-  c(16, 215, 362, 97, 76, 62, 710, 56, 133,28, 62, 413, 194)
Npop <- c(26983, 250930, 829348, 157775, 150467, 352445, 553066, 26751, 75815, 150302, 354993, 3683043, 507218)



d_x<-2       # dimension of theta
d_y<-length(Npart)       # dimension of phi

ncol.print <- function(dat) matrix(as.matrix(dat),ncol=ncol(dat),dimnames=NULL)
PhiC<-ncol.print(as.matrix(read.csv("PhiC.csv",head=F)))

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

post<-function(theta,phi,Y,Z){
  out<-px(phi,Z)+py(Y,theta,phi)
  return(out)
}

MC<-function(init.t,init.p,Y,Z,num_run=1000,burn_in=500){
  theta<-init.t
  phi<-init.p
  sto.theta<-rep(0,(num_run-burn_in)*d_x)
  dim(sto.theta)<-c((num_run-burn_in),d_x)
  sto.phi<-rep(0,(num_run-burn_in)*d_y)
  dim(sto.phi)<-c((num_run-burn_in),d_y)
  sd.p<-apply(PhiC,MARGIN = 2,sd)/10
  for(i in 1:num_run){
    if(i<=burn_in){
      theta.n<-rnorm(2, mean = theta, sd = c(0.01,0.25))
      phi.n<-rtruncnorm(1, a=0, b=1, mean = phi, sd = sd.p)
      rate<-post(theta.n,phi.n,Y,Z)+log(prod(dtruncnorm(phi.n,a=0,b=1,mean=phi,sd=sd.p)))+log(prod(dnorm(theta.n,mean=theta,sd=c(0.01,0.25))))-post(theta,phi,Y,Z)-log(prod(dtruncnorm(phi,a=0,b=1,mean=phi.n,sd=sd.p)))-log(prod(dnorm(theta,mean=theta.n,sd=c(0.01,0.25))))
      alfa<-min(1,exp(rate))
      rpan<-runif(1)
      theta<-theta.n*sign(rpan<=alfa)+theta*sign(rpan>alfa)
      phi<-phi.n*sign(rpan<=alfa)+phi*sign(rpan>alfa)
    }else{
      theta.n<-rnorm(2, mean = theta, sd = c(0.01,0.25))
      phi.n<-rtruncnorm(1, a=0, b=1, mean = phi, sd = sd.p)
      rate<-post(theta.n,phi.n,Y,Z)+log(prod(dtruncnorm(phi.n,a=0,b=1,mean=phi,sd=sd.p)))+log(prod(dnorm(theta.n,mean=theta,sd=c(0.01,0.25))))-post(theta,phi,Y,Z)-log(prod(dtruncnorm(phi,a=0,b=1,mean=phi.n,sd=sd.p)))-log(prod(dnorm(theta,mean=theta.n,sd=c(0.01,0.25))))
      alfa<-min(1,exp(rate))
      rpan<-runif(1)
      theta<-theta.n*sign(rpan<=alfa)+theta*sign(rpan>alfa)
      phi<-phi.n*sign(rpan<=alfa)+phi*sign(rpan>alfa)
      sto.theta[(i-burn_in),]<-theta
      sto.phi[(i-burn_in),]<-phi
    }
    print(c(i,rate,theta))
  }
  out<-cbind(sto.theta,sto.phi)
  return(out)
}

init.t<-c(-2,13)
init.p<-apply(PhiC,MARGIN = 2,median)
Result<-MC(init.t,init.p,Y,Z,num_run=140000,burn_in=0)
write.csv(Result,paste("Martin_mcmc_",task_id,".csv",sep = ""))
