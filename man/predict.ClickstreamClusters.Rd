\name{predict.ClickstreamClusters}
\alias{predict.ClickstreamClusters}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Predicts the cluster for a given Pattern object
}
\description{
Predicts the cluster for a given \code{Pattern} object. Potential clusters need to be identified with the method \code{clusterClickstreams} before predicting the cluster.
}
\usage{
\method{predict}{ClickstreamClusters}(object, pattern, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{object}{
        A \code{ClickstreamClusters} object containing the potential clusters. A \code{ClickstreamClusters} object represents the result of a cluster analysis on a list of clickstreams (see \code{\link{clusterClickstreams}}).
    }
    \item{pattern}{
        The first clicks of a user as \code{Pattern} object.
    }    
    \item{...}{
        Ignored parameters.
    }
}

\value{
    Returns the index of the clusters to which the given \code{Pattern} object most probably belongs to.
}

\author{
Michael Scholz \email{michael.scholz@uni-passau.de}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{clusterClickstreams}}, \code{\link{print.ClickstreamClusters}}
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
clusters<-clusterClickstreams(cls, order=0, centers=2)
pattern<-new("Pattern", sequence=c("h", "c"))
predict(clusters, pattern)
}
