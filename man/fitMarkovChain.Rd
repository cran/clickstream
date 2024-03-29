% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Fitting.r
\name{fitMarkovChain}
\alias{fitMarkovChain}
\title{Fits a List of Clickstreams to a Markov Chain}
\usage{
fitMarkovChain(clickstreamList, order = 1, verbose = TRUE, control = list())
}
\arguments{
\item{clickstreamList}{A list of clickstreams for which a Markov chain is
fitted.}

\item{order}{(Optional) The order of the Markov chain that is fitted from
the clickstreams. Per default, Markov chains with \code{order=1} are fitted.
It is also possible to fit zero-order Markov chains (\code{order=0}) and
higher-order Markov chains.}

\item{verbose}{(Optional) An optimal logical variable to indicate whether warnings
and infos should be printed.}

\item{control}{(Optional) The control list of optimization parameters. Parameter
\code{optimizer} specifies the type of solver used to solve the given
optimization problem. Possible values are "linear" (default) and "quadratic".
Parameter \code{use.lpSolve} determines whether lpSolve or linprog is used as
linear solver.}
}
\value{
Returns a \code{MarkovChain} object.
}
\description{
This function fits a list of clickstreams to a Markov chain. Zero-order,
first-order as well as higher-order Markov chains are supported. For
estimating higher-order Markov chains this function solves the following
linear or quadratic programming problem:\cr \deqn{\min ||\sum_{i=1}^k X-\lambda_i
Q_iX||}{min ||\sum X-\lambda_i Q_iX||} \deqn{\mathrm{s.t.}}{s.t.}
\deqn{\sum_{i=1}^k \lambda_i = 1}{sum \lambda_i = 1} \deqn{\lambda_i \ge
0}{\lambda_i \ge 0} The distribution of states is given as \eqn{X}.
\eqn{\lambda_i} is the lag parameter for lag \eqn{i} and \eqn{Q_i} the
transition matrix.
}
\details{
For solving the quadratic programming problem of higher-order Markov chains,
an augmented Lagrange multiplier method from the package
\code{\link{Rsolnp}} is used.
}
\note{
At least half of the clickstreams need to consist of as many clicks as
the order of the Markov chain that should be fitted.
}
\examples{

# fitting a simple Markov chain
clickstreams <- c("User1,h,c,c,p,c,h,c,p,p,c,p,p,o",
               "User2,i,c,i,c,c,c,d",
               "User3,h,i,c,i,c,p,c,c,p,c,c,i,d",
               "User4,c,c,p,c,d",
               "User5,h,c,c,p,p,c,p,p,p,i,p,o",
               "User6,i,h,c,c,p,p,c,p,c,d")
cls <- as.clickstreams(clickstreams, header = TRUE)
mc <- fitMarkovChain(cls)
show(mc)

}
\references{
This method implements the parameter estimation method presented
in Ching, W.-K. et al.: \emph{Markov Chains -- Models, Algorithms and
Applications}, 2nd edition, Springer, 2013.
}
\seealso{
\code{\link[=MarkovChain-class]{MarkovChain}},
\code{\link[=Rsolnp]{Rsolnp}}
}
\author{
Michael Scholz \email{michael.scholz@th-deg.de}
}
