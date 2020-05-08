\name{CutModel}
\alias{CutModel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Function to set the cut model
}
\description{
 This function combines all probability density and observations involved in the cut distribution. Including two density functions for phi and theta and the observations Z and Y.
}
\usage{
CutModel(px = px, py = py, prox = prox, rprox = rprox, proy = proy, rproy = rproy, Z = Z, Y = Y, d_x = d_x, d_y = d_y, cpp_yes = FALSE, cpp_package = '')
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{px}{
    The density1: Z|phi * phi
}
  \item{py}{
    The unnormalized density2: Y|theta&phi * theta
  }
  \item{prox}{
    The proposal density for phi
  }
  \item{rprox}{
    The proposal samling function for phi
  }
  \item{proy}{
    The proposal density for theta in the auxiliary chain
  }
  \item{rproy}{
    The proposal samling function for theta in the auxiliary chain
  }
  \item{Z}{
    The observations of Z
  }
  \item{Y}{
    The observations of Y
  }
  \item{d_x}{
    The dimension of theta
  }
  \item{d_y}{
    The dimension of phi
  }
  \item{cpp_yes}{
    Logic value, indicate if you are using a C code to write density function px and py and build them within a package XX.
  }
  \item{cpp_package}{
    If cpp_yes == TRUE, name the package XX here.
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  It produce a list of density function and observations. This simplifies the process when these functions and observations are loaded in other functions.
\item{px}{
    The density1: Z|phi * phi
}
  \item{py}{
    The unnormalized density2: Y|theta&phi * theta
  }
  \item{prox}{
    The proposal density for phi
  }
  \item{rprox}{
    The proposal samling function for phi
  }
  \item{proy}{
    The proposal density for theta in the auxiliary chain
  }
  \item{rproy}{
    The proposal samling function for theta in the auxiliary chain
  }
  \item{Z}{
    The observations of Z
  }
  \item{Y}{
    The observations of Y
  }
  \item{d_x}{
    The dimension of theta
  }
  \item{d_y}{
    The dimension of phi
  }
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
function(px = px, py = py, Z = Z, Y = Y, d_x = d_x, d_y = d_y, cpp_yes = FALSE, cpp_package = ''){
  if(cpp_yes){
    rm(px)
    rm(py)
    library(cpp_package,character.only=TRUE)
    return(list(px=px,py=py, Z = Z, Y = Y, d_x = d_x, d_y = d_y))
  }else{
    return(list(px=px,py=py, Z = Z, Y = Y, d_x = d_x, d_y = d_y))
  }
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Cut }% use one of  RShowDoc("KEYWORDS")
