\name{Pattern-class}
\Rdversion{1.1}
\docType{class}
\alias{Pattern-class}
\alias{[,Pattern-method}
\alias{+,Pattern,Pattern-method}
\alias{initialize,Pattern-method}
\alias{show,Pattern-method}

\title{Class \code{"Pattern"}}
\description{
This S4 class describes a click pattern consisting of a sequence of clicks and a probability of occurrence.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("Pattern", sequence, probability, ...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{sequence}:}{Sequence of clicks}
    \item{\code{probability}:}{(Optional) Probability of occurrence}
    \item{\code{absorbingProbabilities}:}{(Optional) Probability of absorbing in any of the absorbing states}
  }
}
\section{Methods}{
  \describe{
    \item{[}{\code{signature(x = "Pattern")}: ... }
    \item{+}{\code{signature(e1 = "Pattern", e2 = "Pattern")}: ... }
    \item{initialize}{\code{signature(.Object = "Pattern")}: ... }
    \item{show}{\code{signature(object = "Pattern")}: ... }
	 }
}
\author{
Michael Scholz \email{michael.scholz@uni-passau.de}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link[=MarkovChain-class]{randomClicks}}
}
\examples{
# show Pattern definition
showClass("Pattern")

# create simple Pattern objects
pattern1<-new("Pattern", sequence=c("h", "c", "p"))
pattern2<-new("Pattern", sequence=c("c", "p", "p"), probability=0.2)
pattern3<-new("Pattern", sequence=c("h", "p", "p"), probability=0.35, 
        absorbingProbabilities=data.frame(d=0.6, o=0.4))
}
\keyword{classes}
