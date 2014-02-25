\name{clickstream-package}
\alias{clickstream-package}
\alias{clickstream}
\docType{package}
\title{Analyzes click streams based on Markov chains
}
\description{
This package allows modeling click streams with Markov chains. It supports to model click streams as zero-order, first-order or higher-order Markov chains. 
}
\details{
\tabular{ll}{
Package: \tab clickstream\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2014-02-25\cr
License: \tab GPL-2\cr
Depends: \tab R (>= 3.0), methods\cr
}

}
\author{
Michael Scholz

Maintainer: Michael Scholz <michael.scholz@uni-passau.de>

}
\references{
For modeling click streams with Markov chains, see Ching, W.-K. et al.: \emph{Markov Chains -- Models, Algorithms and Applications}, 2nd edition, Springer, 2013.
}

\keyword{ click stream }
\keyword{ Markov chain }
\examples{
# fitting a simple Markov chain and predicting the next click
clickstreams<-c("User1,h,c,c,p,c,h,c,p,p,c,p,p,o",
               "User2,i,c,i,c,c,c,d",
               "User3,h,i,c,i,c,p,c,c,p,c,c,i,d",
               "User4,c,c,p,c,d",
               "User5,h,c,c,p,p,c,p,p,p,i,p,o",
               "User6,i,h,c,c,p,p,c,p,c,d")
csf<-tempfile()
writeLines(clickstreams, csf)
cls<-readClickstreams(csf, header=TRUE)
mc<-fitMarkovChain(cls)
startPattern<-new("Pattern", sequence=c("h", "c"))
predict(mc, startPattern)
plot(mc)
}
