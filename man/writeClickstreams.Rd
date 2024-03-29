% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Clickstream.r
\name{writeClickstreams}
\alias{writeClickstreams}
\title{Writes a List of Clickstreams to File}
\usage{
writeClickstreams(
  clickstreamList,
  file,
  header = TRUE,
  sep = ",",
  quote = TRUE
)
}
\arguments{
\item{clickstreamList}{The list of clickstreams to be written.}

\item{file}{The name of the file which the clickstreams are written to.}

\item{header}{A logical flag indicating whether the name of each clickstream
element should be used as first element.}

\item{sep}{The character used to separate clicks (default is \dQuote{,}).}

\item{quote}{A logical flag indicating whether each element of a clickstream
will be surrounded by double quotes (default is \code{TRUE}.}
}
\description{
Writes a list of clickstream to a csv-file.
}
\examples{

clickstreams <- c("User1,h,c,c,p,c,h,c,p,p,c,p,p,o",
               "User2,i,c,i,c,c,c,d",
               "User3,h,i,c,i,c,p,c,c,p,c,c,i,d",
               "User4,c,c,p,c,d",
               "User5,h,c,c,p,p,c,p,p,p,i,p,o",
               "User6,i,h,c,c,p,p,c,p,c,d")
cls <- as.clickstreams(clickstreams, header = TRUE)
clusters <- clusterClickstreams(cls, order = 0, centers = 2)
writeClickstreams(cls, file = "clickstreams.csv", header = TRUE, sep = ",")

# Remove the clickstream file
unlink("clickstreams.csv")

}
\seealso{
\code{\link{readClickstreams}}, \code{\link{clusterClickstreams}}
}
\author{
Michael Scholz \email{michael.scholz@th-deg.de}
}
