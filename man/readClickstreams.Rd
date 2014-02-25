\name{readClickstreams}
\alias{readClickstreams}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Reads a list of click streams from a file
}
\description{
Reads a list of click stream from a csv-file.
}
\usage{
readClickstreams(file, sep = ",", header = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file}{
The name of the file which the click streams are to be read from. Each line of the file appears as one click stream. If it does not contain an absolute path, the file name is relative to the current working directory, \code{\link{getwd}}.
}
  \item{sep}{
The character used to separate clicks (default is \dQuote{,}).
}
  \item{header}{
A logical flag indicating whether the first entry of each line in the file is the name of the click stream user.
}
}

\value{
A list of click streams
}
\author{
Michael Scholz
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{print.Clickstreams}}, \code{\link{randomClickstreams}}
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
