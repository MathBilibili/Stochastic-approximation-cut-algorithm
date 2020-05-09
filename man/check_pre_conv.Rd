\name{check_pre_conv}
\alias{check_pre_conv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to visually check the convergence of auxiliary chain.
}
\description{
As suggested by Liang et al. (2016), when the auxiliary chain has converged, the number of visited auxiliary phi should be equal in theorem. Hence, this function provides a histogram to help visually check if every auxiliary phi has been equally visited in the auxiliary chain.
}
\usage{
check_pre_conv(pre_values,PhiC)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pre_values}{
    The initial condition for the SACut algorithm. It is obtained by running a preliminary chain. Hence, it should be an output from function Preliminary_SACut().
  }
  \item{PhiC}{
    The auxiliary parameter set. It is obtained by calling function LoadNewPhi0() or LoadOldPhi0(). Note that, this should be the same as the PhiC in function Preliminary_SACut().
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
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
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
