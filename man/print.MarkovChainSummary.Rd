\name{print.MarkovChainSummary}
\alias{print.MarkovChainSummary}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Prints the summary of a \code{MarkovChain} object
}
\description{
Prints the summary of a \code{MarkovChain} object. 
}
\usage{
\method{print}{MarkovChainSummary}(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{x}{
        A \code{MarkovChainSummary} object generated with \code{\link[=MarkovChain-class]{summary}}
    }
    
    \item{...}{
        Ignored parameters.
    }
}

\author{
Michael Scholz \email{michael.scholz@uni-passau.de}
}

\seealso{
\code{\link[=MarkovChain-class]{summary}}
}
\examples{
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
print(summary(mc))
}
