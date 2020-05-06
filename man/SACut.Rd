\name{SACut}
\alias{SACut}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to run the Stochastic Approximation Cut (SACut) algorithm.
}
\description{
The Stochastic Approximation Cut (SACut) algorithm is proposed to draw samples from the cut distribution discussed in Liu et al. (2009) and Plummer (2015). The algorithm is divided into two parallel chains. The main chain targets an approximation of the cut distribution; the auxiliary chain forms a proposal distribution used in the main chain.
}
\usage{
SACut(pre_values, PhiC,numrun=1000,burnin=500,thin=1, no=1000,acce_pa=1, sig_dig, filename, Comenvir, CutModel)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pre_values}{
    The initial condition for the SACut algorithm. It is obtained by running a preliminary chain. Hence, it should be a output from function Preliminary_SACut().
  }
  \item{PhiC}{
    The auxiliary parameter set. It is obtained by calling function LoadNewPhi0() or LoadOldPhi0(). Note that, this should be the same as the PhiC in function Preliminary_SACut().
  }
  \item{numrun}{
    The number of total iterations of the SACut algorithm.
  }
  \item{burnin}{
    The number of iterations that are discarded for burn-in.
  }
  \item{thin}{
    This is the thinning parameter for the Markov chain.
  }
  \item{no}{
    The shrink magnitude no determines the sequence of steps used in the stochastic approximation method. It should be set to a large number for a complex problem as large no forces the auxiliary chain to visit each distribution (condition on different phi_0) quickly. It is the t_0 defined in Liang et al. (2007) and Liang et al. (2016). See the reference for details. Note that, this should be the same as the no in function Preliminary_SACut().
  }
  \item{acce_pa}{
    The accelerate parameter acce_pa is set for speedy convergence of the auxiliary chain by increasing the size of the sequence of steps used in the stochastic approximation method. Normally it is set to be 1, but you can try higher value if the auxiliary chain converges slowly. Note that, this should be the same as the acce_pa in function Preliminary_SACut().
  }
  \item{sig_dig}{
    The precision parameter \eqn{\kappa}. It has same dimension as \eqn{\theta}. Note that, this should be the same as the sig_dig in function Preliminary_SACut().
  }
  \item{filename}{
    This is the name of the output file. Note that, it should be a .CSV file (i.e. 'XX.CSV').
  }
  \item{Comenvir}{
    This is the parallel computational environment. It should be an output from function ComEnvir().
  }
  \item{CutModel}{
    The cut model used, it is an output of function CutModel().
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A .CSV file which contains samples of theta and phi is produced at the current path.
}
\references{
Liu, F., M. Bayarri, J. Berger, et al. (2009). Modularization in Bayesian Analysis, with Emphasis on Analysis of Computer Models. Bayesian Analysis 4 (1), 119-150.

Plummer, M. (2015). Cuts in Bayesian Graphical Models. Statistics and Computing 25 (1), 37-43.

Liang, F., C. Liu, and R. J. Carroll (2007). Stochastic Approximation in Monte Carlo Computation. Journal of the American Statistical Association 102 (477), 305-320.

Liang, F., I. H. Jin, Q. Song, and J. S. Liu (2016). An Adaptive Exchange Algorithm for Sampling from Distributions with Intractable Normalizing Constants. Journal of the American Statistical Association 111 (513), 377-393.

}
\author{
Yang Liu
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the

