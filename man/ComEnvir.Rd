\name{ComEnvir}
\alias{ComEnvir}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to set the parallel computational environment
}
\description{
 This function is used to set the parallel environment. For Linux system, a cluster of 'FORK' is created. For Windows system, a cluster of 'PSOCK' is created.
}
\usage{
ComEnvir(is_Unix = TRUE, core_num = 1, clusterExport = list('py','Y','Npop','dmvnorm'))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{is_Unix}{
    Indicate if your system is Linux or Windows
  }
  \item{core_num}{
    The number of cores
  }
  \item{clusterExport}{
    If Windows is used, please state the list of values involved on the master process of the variables to variables of the same names in the global environment of each node. This normalily is a list that includes 'py' (the unnormalized density2: Y|theta&phi \times theta), observation 'Y' and other variables or functions.
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
If is_Unix == TRUE, a list of is_Unix is returned. If is_Unix == FALSE, the list of is_Unix and clusterExport is returned.
%% ...
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
function(is_Unix = TRUE, core_num = 1, clusterExport = list('py','Y','Npop','dmvnorm')){
  if(is_Unix){
    cl<-makeCluster(core_num,type = 'FORK')
    registerDoParallel(cl)
  }else{
    #cl<-makeCluster(detectCores( ))
    cl<-makeCluster(core_num)
    registerDoParallel(cl)
    return(clusterExport)
  }
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the

