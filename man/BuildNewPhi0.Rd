\name{BuildNewPhi0}
\alias{BuildNewPhi0}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to build a new auxiliary parameter \Phi_0 set
}
\description{
This function builds an auxiliary parameter set for parameter phi, following the process discussed in Liang et al. (2016). The idea is to run a MCMC chain for phi only, then select a set of values from the samples of phi which represents the support of distribution p(phi|Z).
}
\usage{
BuildNewPhi0(numrun=1000,burnin=500,numsel=100,CutModel)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{numrun}{
    The number of iterations of the MCMC chain for phi
  }
  \item{burnin}{
    The number of samples that are removed
  }
  \item{numsel}{
    The number of auxiliary parameter phi_0 selected
  }
  \item{CutModel}{
    The cut model used, it is an output of function CutModel().
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A .CSV file which contains auxiliary parameter phi_0 is produced at current path.
}
\references{
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
