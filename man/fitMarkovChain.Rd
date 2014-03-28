\name{fitMarkovChain}
\alias{fitMarkovChain}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fits a list of clickstreams to a Markov chain
}
\description{
This function fits a list of clickstreams to a Markov chain. Zero-order, first-order as well as higher-order Markov chains are supported. For estimating higher-order Markov chains this function solves the following quadratic programming problem:\cr
\deqn{\min ||\sum_{i=1}^k X-\lambda_i Q_iX||}{min ||\sum X-\lambda_i Q_iX||}
\deqn{\mathrm{s.t.}}{s.t.}
\deqn{\sum_{i=1}^k \lambda_i = 1}{sum \lambda_i = 1}
\deqn{\lambda_i \ge 0}{\lambda_i \ge 0}
The distribution of states is given as \eqn{X}. \eqn{\lambda_i} is the lag parameter for lag \eqn{i} and \eqn{Q_i} the transition matrix.
}
\usage{
fitMarkovChain(clickstreamList, order = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{clickstreamList}{
A list of clickstreams for which a Markov chain is fitted.
}
  \item{order}{
(Optional) The order of the Markov chain that is fitted from the clickstreams. Per default, Markov chains with \code{order=1} are fitted. It is also possible to fit zero-order Markov chains (\code{order=0}) and higher-order Markov chains. 
}
}
\details{
For solving the quadratic programming problem of higher-order Markov chains, an augmented Lagrange multiplier method from the package \code{\link{Rsolnp}} is used.
}
\value{
Returns a \code{MarkovChain} object.
}
\references{
This method implements the parameter estimation method presented in Ching, W.-K. et al.: \emph{Markov Chains -- Models, Algorithms and Applications}, 2nd edition, Springer, 2013.
}
\author{
Michael Scholz \email{michael.scholz@uni-passau.de}
}
\note{
At least half of the clickstreams need to consist of as many clicks as the order of the Markov chain that should be fitted.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link[=MarkovChain-class]{MarkovChain}}, \code{\link[=Rsolnp]{Rsolnp}}
}
\examples{
# fitting a simple Markov chain
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
print(mc)
}

