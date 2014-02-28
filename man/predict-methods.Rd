\name{predict-methods}
\docType{methods}
\alias{predict,MarkovChain-method}
\title{Predicts the next click(s) of a user}
\description{
This method predicts the next click(s) of a user. The first clicks of a user are given as \code{Pattern} object. The next click(s) are predicted based on the transition probabilities in the \code{MarkovChain} object. The probability distribution of the next click (n) is estimated as follows:\cr
\deqn{X^{(n)}=\sum_{i=1}^k \lambda_iQ_iX^{(n-i)}}{X^n=sum \lambda_i Q_iX^{n-i}}
}
\section{Methods}{
\describe{

\item{\code{signature(object = "MarkovChain")}}{
%%  ~~describe this method here~~
}
}}

\arguments{
  \item{object}{
A \code{MarkovChain} object used for predicting the next click(s)
}
  \item{startPattern}{
The first clicks of a user as \code{Pattern} object. A \code{Pattern} object with an empty sequence is also possible.
}
    \item{dist}{
(Optional) The number of clicks that should be predicted (default is 1).
}
    \item{ties}{
(Optional) The strategy for handling ties in predicting the next click. Possible strategies are \code{random} (default) and \code{first}.
}
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
