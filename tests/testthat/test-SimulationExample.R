context('SimulationExample')

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
library(progress)

ind_var<<-seq(1,50)

test_that('SimulationExample',{
  Z <- as.numeric(fread('Z.txt')[[1]])
  ind_var<-seq(1,50)
  Y <-  as.numeric(fread('Y.txt')[[1]])

  d_x<-1    # dimension of theta
  d_y<-1      # dimension of phi

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

  proy<-function(t_n,t){
    re<-dnorm(as.numeric(t_n),mean = as.numeric(t),sd=2)
    return(re)
  }

  rproy<-function(t){
    re<-rnorm(1,mean = t, sd=2) %>% t()
    return(re)
  }

  cutmodel <- SACut::CutModel(px = px, py = py, prox = prox, rprox = rprox, proy = proy, rproy = rproy,
                              Z = Z, Y = Y, d_x = d_x, d_y = d_y)

  expect_true(is.list(cutmodel))

  comenvir <- SACut::ComEnvir(is_Unix = FALSE, core_num = 2, clusterExport = list('py','Y','ind_var','dnorm'))

  expect_true(is.list(comenvir))

  SACut::BuildNewPhi0(numrun = 1000, burnin = 0,numsel = 10,CutModel = cutmodel)

  PhiC_test <- read.csv('PhiC.csv')

  expect_true(is.data.frame(PhiC_test))

  unlink('PhiC.csv')

  PhiC_test2 <- SACut::LoadNewPhi0(numsel=10,filename="PhiC_test2.csv",CutModel=cutmodel)

  expect_equal(dim(PhiC_test2),c(10,1))

  PhiC <- SACut::LoadOldPhi0(filename="PhiC_test2.csv")

  init<-list(theta=1,phi=apply(PhiC,MARGIN = 2,median),t=as.matrix(1),I=1)

  PreRun <- SACut::Preliminary_SACut(init=init, PhiC,numrun=250,auxrun=200,no=10,acce_pa=1, sig_dig=4, CutModel=cutmodel)

  expect_true(is.list(PreRun))

  SACut::SACut(pre_values=PreRun, PhiC=PhiC,numrun=120,burnin=10,thin=1, no=10,acce_pa=1, sig_dig=4,
               filename='Result.csv', storage_step=10, Comenvir=comenvir, CutModel=cutmodel)

  result <- read.csv('Result.csv')

  expect_equal(dim(result),c(110,6))

  unlink('Result.csv')


})
