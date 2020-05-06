\name{LoadOldPhi0}
\alias{LoadOldPhi0}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to load an existing auxiliary parameter file without further selection.
}
\description{
This function loads an existing auxiliary parameter set from a .CSV file at current path. Before using this function, please make sure your have removed all column name and row name in your .CSV file!
}
\usage{
LoadOldPhi0(filename="PhiC.csv")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filename}{
    The name of the .CSV file which contains the existing auxiliary parameter set. Please make sure your have removed all column name and row name in your .CSV file!
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
An auxiliary parameter set phi_0.
}
\references{
%% ~put references to the literature/web site here ~
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
## The function is currently defined as
function(filename="PhiC.csv"){
  ncol.print <- function(dat) matrix(as.matrix(dat),ncol=ncol(dat),dimnames=NULL)
  PhiC<-ncol.print(as.matrix(fread(filename,head=F)))
  return(PhiC)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

