\name{randomClickstreams}
\alias{randomClickstreams}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generates a list of click streams
}
\description{
Generates a list of click streams by randomly walking through a given transition matrix.
}
\usage{
randomClickstreams(states, startProbabilities, transitionMatrix, meanLength, n = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{states}{
Names of all possible states.
}
  \item{startProbabilities}{
Start probabilities for all states.
}
  \item{transitionMatrix}{
Matrix of transition probabilities.
}
  \item{meanLength}{
Average length of the click streams.
}
  \item{n}{
Number of click streams to be generated.
}
}
\value{
Returns a list of click streams.
}
\author{
Michael Scholz
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{fitMarkovChain}}, \code{\link{readClickstreams}}, \code{\link{print.Clickstreams}}
}
\examples{
# generate a simple list of click streams
states<-c("a", "b", "c")
startProbabilities<-c(0.2, 0.5, 0.3)
transitionMatrix<-matrix(c(0,0.4,0.6,0.3,0.1,0.6,0.2,0.8,0), nrow=3)
cls<-randomClickstreams(states, startProbabilities, transitionMatrix, meanLength=5, n=10)
print(cls)
}
