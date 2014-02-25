\name{fitMarkovChain}
\alias{fitMarkovChain}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fits a list of click streams to a Markov chain
}
\description{
This function fits a list of click streams to a Markov chain. Zero-order, first-order as well as higher-order Markov chains are supported. For estimating higher-order Markov chains this function solves the following linear programming problem:\cr
\deqn{\min w}{min w}
\deqn{\mathrm{s.t.}}{s.t.}
\deqn{w \ge X-[QX]\lambda}{w \ge X-[QX]\lambda}
\deqn{w \ge -X+[QX]\lambda}{w \ge -X+[QX]\lambda}
\deqn{w \ge 0}{w \ge 0}
\deqn{\sum_{i=1}^k \lambda_i = 1}{sum \lambda_i = 1}
\deqn{\lambda_i \ge 0}{\lambda_i \ge 0}
}
\usage{
fitMarkovChain(clickstreamList, order = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{clickstreamList}{
A list of click streams for which a Markov chain is fitted.
}
  \item{order}{
(Optional) The order of the Markov chain that is fitted from the click streams. Per default, Markov chains with \code{order=1} are fitted. It is also possible to fit zero-order Markov chains (\code{order=0}) and higher-order Markov chains. 
}
}
\details{
For solving the linear programming problem of higher-order Markov chains, an augmented Lagrange multiplier method from the package \code{\link{Rsolnp}} is used.
}
\value{
Returns a \code{MarkovChain} object.
}
\references{
This method implements the parameter estimation method presented in Ching, W.-K. et al.: \emph{Markov Chains -- Models, Algorithms and Applications}, 2nd edition, Springer, 2013.
}
\author{
Michael Scholz
}
\note{
At least half of the click streams need to consist of as many clicks as the order of the Markov chain that should be fitted.
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

