#specification of default data
Z <- c(7, 6, 10, 10, 1, 1, 10, 4, 35, 0, 10, 8, 4)
Npart <- c(111, 71, 162, 188, 145, 215, 166, 37, 173,143, 229, 696, 93)
Y <-  c(16, 215, 362, 97, 76, 62, 710, 56, 133,28, 62, 413, 194)
Npop <- c(26983, 250930, 829348, 157775, 150467, 352445, 553066, 26751, 75815, 150302, 354993, 3683043, 507218)

#setting default likelihood function \times prior
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

#setting the default proposal distribution for phi
#proposal density
prox<-function(phi_n,phi){
  out<-sum(log(dtruncnorm(phi_n, a=0, b=1, mean = phi, sd = 0.005)))
  return(out)
}
#samling function
rprox<-function(phi){
  out<-rtruncnorm(1, a=0, b=1, mean = phi, sd = 0.005)
  return(out)
}

#setting the default proposal distribution for theta in the internal chain
#proposal density
denT<-function(t_n,t,P){
  re<-dmvnorm(as.numeric(t_n),mean = as.numeric(t),sigma = matrix(c(0.1648181,-0.3979341,-0.3979341,2.737874),ncol=2,nrow=2)/10)*P
  return(re)
}

#samling function
tranT<-function(t){
  re<-rmvnorm(1,mean = t,sigma = matrix(c(0.1648181,-0.3979341,-0.3979341,2.737874),ncol=2,nrow=2)/10) %>% t()
  return(re)
}


#setting the default dimension of \theta and \phi
d_x<-2       # dimension of theta
d_y<-13       # dimension of phi

CutModel<-function(px = px, py = py, prox = prox, rprox = rprox, denT = denT, tranT = tranT, Z = Z, Y = Y, d_x = d_x, d_y = d_y, cpp_yes = FALSE, cpp_package = ''){
  if(cpp_yes){
    rm(px)
    rm(py)
    library(cpp_package,character.only=TRUE)
    return(list(px=px,py=py, prox = prox, rprox = rprox, denT = denT, tranT = tranT, Z = Z, Y = Y, d_x = d_x, d_y = d_y))
  }else{
    return(list(px=px,py=py, prox = prox, rprox = rprox, denT = denT, tranT = tranT, Z = Z, Y = Y, d_x = d_x, d_y = d_y))
  }
}
