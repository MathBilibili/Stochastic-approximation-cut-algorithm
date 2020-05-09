# SACut: an R package for Stochastic Approximation Cut Algorithm (SACut)
[![](https://travis-ci.com/MathBilibili/Stochastic-approximation-cut-algorithm.svg?branch=master)](https://travis-ci.com/MathBilibili/Stochastic-approximation-cut-algorithm)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

This is an [`R`][R] package to conduct the Stochastic Approximation Cut Algorithm. It also contains code that replicates the results in the paper.
<img align="right" width="200" height="200" src="https://user-images.githubusercontent.com/24710640/81212775-495b2880-8fcd-11ea-9319-52ac4fd15f4f.png">

## Authors
Yang Liu and Robert Goudie

MRC Biostatistics Unit, University of Cambridge

## Introduction of Stochastic Approximation Cut Algorithm
The potential effect of partial misspecification of Bayesian modelling is a concern. Recent studies have proposed the idea of modularized models, and the cut model is proposed to prevent feedback of information from the suspect module. This leads to the cut posterior distribution which normally does not have a closed-form. Previous studies have proposed algorithms to sample from this distribution, but these algorithms have unclear theoretical convergence properties. To address this convergence problem, the novel Stochastic Approximation Cut algorithm (SACut) is proposed as an alternative. The algorithm is divided into two parallel chains. The main chain targets an approximation of the cut distribution; the auxiliary chain forms a proposal distribution used in the main chain. We prove convergence of the samples drawn by the proposed algorithm and present the exact limit. Although SACut is biased, since the main chain does not target the exact cut distribution, we prove this bias can be reduced geometrically by increasing a user-chosen tuning parameter. In addition, parallel processing can be easily adopted for SACut, unlike existing algorithms. This greatly reduces computation time.

## Citation

Please cite the following paper when using SACut:

Yang Liu, Robert Goudie. Stochastic Approximation Cut Algorithm for Inference in Modularised Bayesian Models. DOI: XXX

## Installation
Simply download and install the SACut package from Github.
```r
devtools::install_github('MathBilibili/Stochastic-approximation-cut-algorithm')
```

SACut depends on the packages `modeltools`, `flexclust`, `compiler`, `doParallel`, `dplyr`, `tidyr`, and `data.table` which can be installed via:
```r
install.packages(c("modeltools", "flexclust", "compiler", "doParallel", "dplyr", "tidyr", "data.table"))
```

Additionally you may need `truncnorm` and `mvtnorm` to write likelihood or proposal distribution of truncated normal distribution and multivariate normal distribution, and `coda` to assess the convergence of the Markov chain.
```r
install.packages(c("truncnorm", "mvtnorm", "coda"))
```

## Example of usage
Here we explore the usage of SACut on the HPV study discussed in [Plummer (2015)][Plummer2015]. First, we load the data and the dimension of parameter `theta` and `phi`
```r
Z <- c(7, 6, 10, 10, 1, 1, 10, 4, 35, 0, 10, 8, 4)
Npart <- c(111, 71, 162, 188, 145, 215, 166, 37, 173,143, 229, 696, 93)
Y <-  c(16, 215, 362, 97, 76, 62, 710, 56, 133,28, 62, 413, 194)
Npop <- c(26983, 250930, 829348, 157775, 150467, 352445, 553066, 26751, 75815, 150302, 354993, 3683043, 507218)

d_x<-2       # dimension of theta
d_y<-13       # dimension of phi
```
where `Z` is the number of people with HPV infection out of a sample of `NPart`, and `Y` is the number of cancer cases from `T=0.001*Npop` person-years at 13 cities. Two modules are difined as:

<img src="https://latex.codecogs.com/gif.latex?Z_i\sim&space;\mathbf{B}(N_i,\varphi_i)" title="Z_i\sim \mathbf{B}(N_i,\varphi_i)" />

<img src="https://latex.codecogs.com/gif.latex?Y_i\sim&space;\mathbf{Poisson}(T_i(\exp(\theta_1&plus;\theta_2\varphi_i)))" title="Y_i\sim \mathbf{Poisson}(T_i(\exp(\theta_1+\theta_2\varphi_i)))" />

SACut is applied to prevent the feedback from second module to the estimation of parameter `phi`. The log-density of the joint distribution (likelihood and prior) for first module is loaded:
```r
px<-function(phi,Z){
    PbZ<-rbind(phi,Z,Npart)
    out<-sum(apply(PbZ,FUN=function(y){dbeta(x=y[1],shape1 = (1+y[2]), shape2 = (1+y[3]-y[2]), log = T)},MARGIN = 2))
    return(out)
  }
```
and log-density of the joint distribution (likelihood and prior) for the second module is loaded:
```r
py<-function(Y,theta,phi){
    PbY<-rbind(phi,Y,Npop)
    out<- log(1/200)+sum(apply(PbY,FUN=function(y){dpois(x=y[2],lambda = (y[3]*0.001*exp(theta[1]+theta[2]*y[1])),log = T)},MARGIN = 2))
    return(out)
  }
```
We use the truncated normal distribution as the proposal distribution for parameter `phi` and the log-density and random generation are:
```r
prox<-function(phi_n,phi){
  out<-sum(log(dtruncnorm(phi_n, a=0, b=1, mean = phi, sd = 0.005))) 
  return(out) 
}

rprox<-function(phi){
  out<-rtruncnorm(1, a=0, b=1, mean = phi, sd = 0.005)
  return(out)
}
```
In the auxiliary chain, we use the multivariate normal distribution as the proposal distribution for auxiliary parameter `theta` (note that, this is not the parameter `theta` in the main chain), the density (not log) and random generation are:
```r
proy<-function(theta_n,theta){
  re<-dmvnorm(as.numeric(theta_n),mean = as.numeric(theta),sigma = matrix(c(0.1648181,-0.3979341,-0.3979341,2.737874),ncol=2,nrow=2)/10)
  return(re)
}

rproy<-function(theta){
  re<-rmvnorm(1,mean = theta,sigma = matrix(c(0.1648181,-0.3979341,-0.3979341,2.737874),ncol=2,nrow=2)/10) %>% t()
  return(re)
}
```
Now we have defined every components of the cut distribution, then call function `CutModel` to build the cut model:
```r
cutmodel <- SACut::CutModel(px = px, py = py, prox = prox, rprox = rprox, proy = proy, rproy = rproy,
Z = Z, Y = Y, d_x = d_x, d_y = d_y)
```
Here we do not import `px` and `py` from other package. In the case that they are written and built within other R package for high computational speed:
```r
devtools::install_github('MathBilibili/Stochastic-approximation-cut-algorithm/Examples/Example3/Plummer2015Example/')
library(Plummer2015Example)

cutmodel <- SACut::CutModel(px = px, py = py, prox = prox, rprox = rprox, proy = proy, rproy = rproy,
Z = Z, Y = Y, d_x = d_x, d_y = d_y, cpp_yes = TRUE, cpp_package = 'Plummer2015Example')
```

After setting the cut model, we set the parallel environment for the computation by function `ComEnvir`. For a Windows device with 4 core, a cluster of `PSOCK` is created. We also need claim the list of variables that are exported to clusters (no need for Linux).
```r
comenvir <- SACut::ComEnvir(is_Unix = FALSE, core_num = 4, clusterExport = list('py','Y','Npop','dmvnorm'))
```

We then load the  existing auxiliary parameter set `phi0` from a `.CSV` file by function `LoadOldPhi0`, the file (`PhiC.csv`) can be download from <https://github.com/MathBilibili/Stochastic-approximation-cut-algorithm/tree/master/Examples/Example3>. Note that, this parameter set can be created by calling function `BuildNewPhi0`.
```r
PhiC <- LoadOldPhi0(filename="PhiC.csv")
```

As suggested by [Liang et al. (2016)][Liang2016], the auxiliary chain is conducted solely with sufficient iterations so that it is converged when the main Markov chain (external chain) start to run. Hence, we conduct a preliminary run with 1501000 iterations by calling function `Preliminary_SACut`. The shrink magnitude `no` is set to be large enough to ensure that the auxiliary chain can completely go through every `phi0` before it converges, `acce_pa` is set to be 10 to speed up the convergence. The precision parameter is set to be 3 for `theta_1` and 2 for `theta_2`.
```r
init<-list(theta=c(-2,13),phi=apply(PhiC,MARGIN = 2,median),t=as.matrix(c(-2,13)),I=1)

PreRun <- SACut::Preliminary_SACut(init=init, PhiC,numrun=1501000,auxrun=1500000,no=20000,acce_pa=10, sig_dig=c(3,2), CutModel=cutmodel)
```

Finally, we are able to run the auxiliary chain and the main chain in parallel by calling the function `SACut`. The total number of iterations is 140000 and we retain only every 100 sample after discarding the first 40000 samples. The result is stored in file `Result.csv`.
```r
SACut::SACut(pre_values=PreRun, PhiC=PhiC,numrun=140000,burnin=40000,thin=100, no=20000,acce_pa=10, sig_dig=c(3,2),
filename='Result.csv', Comenvir=comenvir, CutModel=cutmodel)
```



[R]: http://www.r-project.org "The R Project for Statistical Computing"
[Plummer2015]:https://link.springer.com/article/10.1007/s11222-014-9503-z "Cuts in Bayesian graphical models"
[Liang2016]:https://www.tandfonline.com/doi/full/10.1080/01621459.2015.1009072 "An Adaptive Exchange Algorithm for Sampling From Distributions With Intractable Normalizing Constants"
