\name{print.Clickstreams}
\alias{print.Clickstreams}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Prints a list of click streams
}
\description{
Prints a list of click streams. 
}
\usage{
\method{print}{Clickstreams}(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{x}{
        A list of click streams.
    }
    
    \item{...}{
        Ignored parameters.
    }
}

\author{
Michael Scholz
}


\seealso{
\code{\link{readClickstreams}}, \code{\link{randomClickstreams}}
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
print(cls)
}
