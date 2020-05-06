\name{Preliminary_SACut}
\alias{Preliminary_SACut}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to run a preliminary chain for the SACut algorithm.
}
\description{
As suggested by Liang et al. (2016), the auxiliary chain is conducted solely with sufficient iterations so that it is converged when the main Markov chain (external chain) start to run. Normally, this is followed by some further iterations when the auxiliary chain and main chain are run in parallel so a set of auxiliary variables is stored as the initial condition for the SACut algorithm.
}
\usage{
Preliminary_SACut(init=list(theta=c(-2,13),phi=rep(1,13),t=as.matrix(c(-2,13)),I=1),
                            PhiC,numrun=1000,auxrun=500,no=1000,acce_pa=1, sig_dig, CutModel)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{init}{
    The initial condition for the preliminary chain. It should be a list defined as list(theta=c(...),phi=c(...),t=as.matrix(c(...)),I=1), where theta is the initial value of parameter \eqn{\theta}, phi is the initial value of parameter phi, t is the initial value of auxiliary theta (with same dimension as \eqn{\theta}), and I is the index of initial value of auxiliary phi_0.
  }
  \item{PhiC}{
    The auxiliary parameter set. It is obtained by calling function LoadNewPhi0() or LoadOldPhi0().
  }
  \item{numrun}{
    The number of total iterations of the preliminary chain.
  }
  \item{auxrun}{
    The number of iterations of the auxiliary chain when it is conducted solely. It should be less than numrun. You can set auxrun=numrun-100.
  }
  \item{no}{
    The shrink magnitude no determines the sequence of steps used in the stochastic approximation method. It should be set to a large number for a complex problem as large no forces the auxiliary chain to visit each distribution (condition on different phi_0) quickly. It is the t_0 defined in Liang et al. (2007) and Liang et al. (2016). See the reference for details.
  }
  \item{acce_pa}{
    The accelerate parameter acce_pa is set for speedy convergence of the auxiliary chain by increasing the size of the sequence of steps used in the stochastic approximation method. Normally it is set to be 1, but you can try higher value if the auxiliary chain converges slowly.
  }
  \item{sig_dig}{
    The precision parameter \eqn{\kappa}. It has same dimension as \eqn{\theta}.
  }
  \item{CutModel}{
    The cut model used, it is an output of function CutModel().
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A list of parameters used as the initial conditions for the SACut algorithm.
}
\references{
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
% R documentation directory.

