library(unbiasedmcmc)
library(doParallel)
library(data.table)
library(foreach)
library(dplyr)
library(tidyr)

is_Unix<-TRUE
core_num<-10

if(is_Unix){
  cl<-makeCluster(core_num,type = 'FORK')  
  registerDoParallel(cl)
}else{
  #cl<-makeCluster(detectCores( ))
  cl<-makeCluster(core_num)  
  registerDoParallel(cl)
}

Z <- as.numeric(fread('Z.txt')[[1]])
ind_var <-  as.numeric(fread('ind_var.txt')[[1]])
#theta_X <- matrix(0,nrow = 500, ncol = 20)
#for(i in 1:500){
#  theta_X[i,] <- rnorm(20, mean = 2, sd = 1)
#}
#write.csv(theta_X,file = 'theta_X.csv')

theta_X <- as.matrix(read.csv('theta_X.csv'))

#theta <- rep(0,20)
#for(i in 1:20){
#  theta[i]<-sin(i)
#}
#Y<-rep(0,500)
#for(i in 1:500){
#  Y[i] <- rnorm(1, mean = (sum(theta*theta_X[i,])+ind_var[i]),sd=3)
#}
#write.table(Y,file = 'Ycandidate.txt',col.names = F,row.names = F,quote = F)
Y <-  as.numeric(fread('Y.txt')[[1]])

d_x<-20    # dimension of theta
d_y<-1      # dimension of phi

num_phi<-50000

print(paste('Number of phi sampled is',num_phi))

#unnormalizing density2 Y|theta&phi \times theta
py<-function(Y,theta,phi){
  PbY<-rbind(Y,ind_var,t(theta_X))
  out<-sum(apply(PbY,FUN=function(y){dnorm(y[1],mean = sum(theta*y[3:(d_x+2)])+phi*y[2],sd = 3,log = T)},MARGIN = 2))+sum(dnorm(theta,mean = 0,sd=100,log = T))
  return(out)
}

#posterior sampler phi|Z (prior: normal(0,1))
post_phi<-function(n,Z){
  mean <- sum(Z)*(1/(1+length(Z)))
  var <- 1/(1+length(Z))
  re<-rnorm(n,mean = mean,sd=sqrt(var))
  return(re)
}

q1<-Sys.time()

samples_phi <- post_phi(num_phi,Z)

meetingtimes <- rep(0, num_phi)

if(is_Unix==FALSE)(clusterExport(cl, list('get_mh_kernels','sample_meetingtime','sample_coupled_chains','Y','ind_var','dnorm','theta_X','d_x'), envir=environment()))


meetingtimes <- foreach(i = icount(num_phi)) %dopar% {
  phi<-samples_phi[i]
  
  target <- function(theta){
    evals <- py(Y,theta,phi)
    return(evals)
  }
  
  #kernels <- get_mh_kernels(target, (1e-10)*diag(20))
  kernels <- get_mh_kernels(target, diag(20)*(10^(-5)))
  # Markov kernel of the chain
  single_kernel <- kernels$single_kernel
  # Markov kernel of the coupled chain
  coupled_kernel <- kernels$coupled_kernel
  # initial distribution, towards the right-most mode of the target
  rinit <- function(){
    chain_state <- rnorm(20, mean = 1, sd = 0.00001)
    current_pdf <- target(chain_state)
    return(list(chain_state = chain_state, current_pdf = current_pdf))
  }
  
  #meetingtimes[i] <- sample_meetingtime(single_kernel, coupled_kernel, rinit)$meetingtime
  sample_meetingtime(single_kernel, coupled_kernel, rinit)$meetingtime
  
} %>% as.numeric()

q2<-Sys.time()
print(q2-q1)

write.csv(meetingtimes,file = 'meetingtimes.csv')

pdf(file="meetingtimes.pdf")
hist(meetingtimes, breaks = 100, col="violet")
dev.off()

max_iter <- quantile(meetingtimes,0.95)[[1]]*10 %>% round()
print(paste('0.95 upper CI meetingtimes is',max_iter/10))
print(paste('0.99 upper CI meetingtimes is',quantile(meetingtimes,0.99)[[1]]))
coupledchains <- list()

if(is_Unix==FALSE)(clusterExport(cl, list('max_iter','sample_coupled_chains','H_bar'), envir=environment()))


coupledEst <- foreach(irep = icount(num_phi)) %dopar% {
  phi<-samples_phi[irep]
  
  target <- function(theta){
    evals <- py(Y,theta,phi)
    return(evals)
  }
  
  #kernels <- get_mh_kernels(target, (1e-10)*diag(20))
  kernels <- get_mh_kernels(target, diag(20)*(10^(-5)))
  # Markov kernel of the chain
  single_kernel <- kernels$single_kernel
  # Markov kernel of the coupled chain
  coupled_kernel <- kernels$coupled_kernel
  # initial distribution, towards the right-most mode of the target
  
  rinit <- function(){
    chain_state <- rnorm(20, mean = 1, sd = 0.00001)
    current_pdf <- target(chain_state)
    return(list(chain_state = chain_state, current_pdf = current_pdf))
  }
  
  coupled_chain <- sample_coupled_chains(single_kernel, coupled_kernel, rinit, m = max_iter, lag = 1)
  H_bar(coupled_chain, h = function(x) x, k = round(max_iter*0.5), m = max_iter)
}

q3<-Sys.time()

theta_est <- rep(0,d_x)

for(i in 1:length(coupledEst)){
  theta_est <- theta_est + coupledEst[[i]]
}

theta_est<-theta_est/length(coupledEst)

true_theta <- rep(0,d_x)
for(i in 1:d_x){
  true_theta[i]<-sin(i)
}

print(q3-q1)
print(theta_est)
print(paste('MSE is',sum((theta_est-true_theta)^2)))

stopCluster(cl)
