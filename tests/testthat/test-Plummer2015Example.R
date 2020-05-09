context('Plummer2015Example')

library(modeltools)
library(flexclust)
library(compiler)
library(truncnorm)
library(mvtnorm)
library(foreach)
library(iterators)
library(doParallel)
library(dplyr)
library(tidyr)
library(data.table)

test_that('Plummer2015Example',{
  Z <- c(7, 6, 10, 10, 1, 1, 10, 4, 35, 0, 10, 8, 4)
  Npart <- c(111, 71, 162, 188, 145, 215, 166, 37, 173,143, 229, 696, 93)
  Y <-  c(16, 215, 362, 97, 76, 62, 710, 56, 133,28, 62, 413, 194)
  Npop <- c(26983, 250930, 829348, 157775, 150467, 352445, 553066, 26751, 75815, 150302, 354993, 3683043, 507218)

  d_x<-2       # dimension of theta
  d_y<-13       # dimension of phi

  px<-function(phi,Z){
    PbZ<-rbind(phi,Z,Npart)
    out<-sum(apply(PbZ,FUN=function(y){dbeta(x=y[1],shape1 = (1+y[2]), shape2 = (1+y[3]-y[2]), log = T)},MARGIN = 2))
    return(out)
  }

  py<-function(Y,theta,phi){
    PbY<-rbind(phi,Y,Npop)
    out<- log(1/200)+sum(apply(PbY,FUN=function(y){dpois(x=y[2],lambda = (y[3]*0.001*exp(theta[1]+theta[2]*y[1])),log = T)},MARGIN = 2))
    return(out)
  }

  prox<-function(phi_n,phi){
    out<-sum(log(dtruncnorm(phi_n, a=0, b=1, mean = phi, sd = 0.005)))
    return(out)
  }

  rprox<-function(phi){
    out<-rtruncnorm(1, a=0, b=1, mean = phi, sd = 0.005)
    return(out)
  }

  proy<-function(theta_n,theta){
    re<-dmvnorm(as.numeric(theta_n),mean = as.numeric(theta),sigma = matrix(c(0.1648181,-0.3979341,-0.3979341,2.737874),ncol=2,nrow=2)/10)
    return(re)
  }

  rproy<-function(theta){
    re<-rmvnorm(1,mean = theta,sigma = matrix(c(0.1648181,-0.3979341,-0.3979341,2.737874),ncol=2,nrow=2)/10) %>% t()
    return(re)
  }

  cutmodel <- SACut::CutModel(px = px, py = py, prox = prox, rprox = rprox, proy = proy, rproy = rproy,
                              Z = Z, Y = Y, d_x = d_x, d_y = d_y)

  expect_true(is.list(cutmodel))

  #num_core <- detectCores()

  comenvir <- SACut::ComEnvir(is_Unix = FALSE, core_num = 1, clusterExport = list('py','Y','Npop','dmvnorm'))

  expect_true(is.list(comenvir))

  SACut::BuildNewPhi0(numrun = 1000, burnin = 0,numsel = 10,CutModel = cutmodel)

  PhiC_test <- read.csv('PhiC.csv')

  expect_true(is.data.frame(PhiC_test))

  unlink('PhiC.csv')

  PhiC_test2 <- SACut::LoadNewPhi0(numsel=10,filename="PhiC_test.csv",CutModel=cutmodel)

  expect_equal(dim(PhiC_test2),c(10,13))

  PhiC <- SACut::LoadOldPhi0(filename="PhiC_test.csv")

  expect_true(is.matrix(PhiC))

  init<-list(theta=c(-2,13),phi=apply(PhiC,MARGIN = 2,median),t=as.matrix(c(-2,13)),I=1)

  PreRun <- SACut::Preliminary_SACut(init=init, PhiC,numrun=250,auxrun=200,no=10,acce_pa=1, sig_dig=c(3,2), CutModel=cutmodel)

  expect_true(is.list(PreRun))

  expect_silent(check_pre_conv(PreRun,PhiC))

  SACut::SACut(pre_values=PreRun, PhiC=PhiC,numrun=100,burnin=0,thin=1, no=10,acce_pa=1, sig_dig=c(3,2),
               filename='Result.csv', Comenvir=comenvir, CutModel=cutmodel)

  result <- read.csv('Result.csv')

  expect_true(is.data.frame(result))

  expect_equal(dim(result),c(100,20))

  unlink('Result.csv')

})
