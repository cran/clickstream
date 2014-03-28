\name{summary-methods}
\docType{methods}
\alias{summary-methods}
\alias{summary,MarkovChain-method}
\title{Prints a summary of a MarkovChain object}
\description{
 Generates a summary for a given \code{MarkovChain} object
}
\section{Methods}{
\describe{

\item{\code{signature(object = "MarkovChain")}}{
}
}}

\value{
Returns a \code{MarkovChainSummary} object.

\item{\code{desc}}{A short description of the \code{MarkovChain} object.}
\item{\code{observations}}{The number of observations from which the \code{MarkovChain} object has been fitted.}
\item{\code{k}}{The number of estimation parameters.}
\item{\code{logLikelihood}}{The LogLikelihood of the \code{MarkovChain} object.}
\item{\code{aic}}{Akaike's Information Criterion for the \code{MarkovChain} object}
\item{\code{bic}}{Bayesian Information Criterion for the \code{MarkovChain} object}
}
\author{
Michael Scholz \email{michael.scholz@uni-passau.de}
}
\keyword{methods}
