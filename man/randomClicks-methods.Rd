\name{randomClicks-methods}
\docType{methods}
\alias{randomClicks}
\alias{randomClicks,MarkovChain-method}
\title{Generates a sequence of clicks}
\description{
Generates a sequence of clicks by randomly walking through the transition graph of a given \code{MarkovChain} object.
}
\section{Methods}{
\describe{

\item{\code{signature(object = "MarkovChain")}}{
%%  ~~describe this method here~~
}
}}

\arguments{
  \item{object}{
A \code{MarkovChain} object used for generating the next click(s)
}
  \item{startPattern}{
The first clicks of a user as \code{Pattern} object. A \code{Pattern} object with an empty sequence is also possible.
}
    \item{dist}{
(Optional) The number of clicks that should be generated (default is 1).
}
}
\author{
Michael Scholz \email{michael.scholz@uni-passau.de}
}
\seealso{
\code{\link{fitMarkovChain}}
}
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
}
\keyword{methods}
