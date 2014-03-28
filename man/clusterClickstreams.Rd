\name{clusterClickstreams}
\alias{clusterClickstreams}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Perform k-means clustering on a list of clickstreams
}
\description{
Perform k-means clustering on a list of clickstreams. For each clickstream a transition matrix of a given order is computed. These transition matrices are used as input for performing k-means clustering.
}
\usage{
clusterClickstreams(clickstreamList, order = 0, centers, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{clickstreamList}{
A list of clickstreams for which the cluster analysis is performed.
}
  \item{order}{
The order of the transition matrices used as input for clustering (default is 0). 
}
  \item{centers}{
The number of clusters.
}
  \item{...}{
Additional parameters for k-means clustering (see \code{\link{kmeans}}).
}
}

\value{
This method returns a \code{ClickstreamClusters} object (S3-class). It is a list with the following components:
    \item{clusters}{
        A list of \code{Clickstream} objects representing the resulting clusters.    
    }
    \item{centers}{
        A matrix of cluster centres.
    }
    \item{totss}{
        The total sum of squares.
    }
    \item{withinss}{
        Vector of within-cluster sum of squares, one component per cluster.
    }
    \item{tot.withinss}{
        Total within-cluster sum of squares, i.e., \code{sum(withinss)}.
    }
    \item{betweenss}{
        The between-cluster sum of squares, i.e., \code{totss - tot.withinss}.
    }
}
\author{
Michael Scholz \email{michael.scholz@uni-passau.de}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{print.ClickstreamClusters}}, \code{\link{summary.ClickstreamClusters}}
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
print(clusters)
}
