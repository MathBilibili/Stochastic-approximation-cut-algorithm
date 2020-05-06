\name{LoadNewPhi0}
\alias{LoadNewPhi0}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to load an existing auxiliary parameter file and select a smaller auxiliary parameter set.
}
\description{
This function loads an existing auxiliary parameter set from a .CSV file at current path. It then selects a smaller auxiliary parameter set from it based on the Max-Min procedure introduced in Liang et al. (2016). Before using this function, please make sure your have removed all column name and row name in your .CSV file!
}
\usage{
LoadNewPhi0(numsel=100,filename="PhiC.csv",CutModel)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{numsel}{
    The number of auxiliary parameter phi_0 selected
  }
  \item{filename}{
    The name of the .CSV file which contains the existing auxiliary parameter set. Please make sure your have removed all column name and row name in your .CSV file!
  }
  \item{CutModel}{
    The cut model used, it is an output of function CutModel().
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
An auxiliary parameter set phi_0.
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

