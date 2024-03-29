% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Pattern.r
\docType{methods}
\name{initialize,Pattern-method}
\alias{initialize,Pattern-method}
\title{Creates a new \code{Pattern} object}
\usage{
\S4method{initialize}{Pattern}(.Object, sequence, probability, absorbingProbabilities, ...)
}
\arguments{
\item{.Object}{Pattern (name of the class)}

\item{sequence}{Click sequence}

\item{probability}{Probability for the click sequence}

\item{absorbingProbabilities}{Probabilities that the sequence will finally end in one of the absorbing states}

\item{...}{Further arguments for the \code{CallNextMethod} function}
}
\description{
Creates a new \code{Pattern} object
}
\section{Methods}{
 \describe{
\item{list("signature(sequence = \"character\", probability = \"numeric\", absorbingProbabilities = \"numeric"))}{Creates a new \code{Pattern} object.}
}
}

\author{
Michael Scholz \email{michael.scholz@th-deg.de}
}
\keyword{methods}
