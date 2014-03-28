\name{MarkovChain-class}
\Rdversion{1.1}
\docType{class}
\alias{MarkovChain-class}
\alias{plot,MarkovChain-method}
\alias{print,MarkovChain-method}
\alias{show,MarkovChain-method}

\title{Class \code{"MarkovChain"}}
\description{
This S4 class describes \code{MarkovChain} objects.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("MarkovChain", ...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{states}:}{Names of the states}
    \item{\code{order}:}{Order of the Markov chain}
    \item{\code{transitions}:}{List of transition probability matrices}
    \item{\code{lambda}:}{Lag variables}    
    \item{\code{logLikelihood}:}{LogLikelihood of the fitted Markov chain}
    \item{\code{observations}:}{Number of observations from which the Markov chain has been fitted}
    \item{\code{start}:}{Start probabilities}
    \item{\code{end}:}{End probabilities}
    \item{\code{transientStates}:}{Names of the transient states}
    \item{\code{absorbingStates}:}{Names of the absorbing states}
    \item{\code{absorbingProbabilities}:}{Probabilities for the states to be absorbed by a specific absorbing state}
  }
}
\section{Methods}{
  \describe{
    \item{absorbingStates}{\code{signature(object = "MarkovChain")}: ... }
    \item{plot}{\code{signature(x = "MarkovChain")}: ... }
    \item{predict}{\code{signature(object = "MarkovChain")}: ... }
    \item{print}{\code{signature(x = "MarkovChain")}: ... }
    \item{randomClicks}{\code{signature(object = "MarkovChain")}: ... }
    \item{show}{\code{signature(object = "MarkovChain")}: ... }
    \item{states}{\code{signature(object = "MarkovChain")}: ... }
    \item{summary}{\code{signature(object = "MarkovChain")}: ... }
    \item{transientStates}{\code{signature(object = "MarkovChain")}: ... }
	 }
}

\author{
Michael Scholz \email{michael.scholz@uni-passau.de}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{fitMarkovChain}}
}
\examples{
# show MarkovChain definition
showClass("MarkovChain")

# fit a simple Markov chain from a list of click streams
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
\keyword{classes}
