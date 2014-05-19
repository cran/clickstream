#' Reads a list of clickstreams from a file
#' 
#' Reads a list of clickstream from a csv-file.
#' 
#' 
#' @param file The name of the file which the clickstreams are to be read from.
#' Each line of the file appears as one click stream. If it does not contain an
#' absolute path, the file name is relative to the current working directory,
#' \code{\link{getwd}}.
#' @param sep The character used to separate clicks (default is \dQuote{,}).
#' @param header A logical flag indicating whether the first entry of each line
#' in the file is the name of the clickstream user.
#' @return A list of clickstreams. Each element is a vector of characters
#' representing the clicks. The name of each list element is either the header
#' of a clickstream file or a unique number.
#' @author Michael Scholz \email{michael.scholz@@uni-passau.de}
#' @seealso \code{\link{print.Clickstreams}}, \code{\link{randomClickstreams}}
#' @examples
#' 
#' clickstreams<-c("User1,h,c,c,p,c,h,c,p,p,c,p,p,o",
#'                "User2,i,c,i,c,c,c,d",
#'                "User3,h,i,c,i,c,p,c,c,p,c,c,i,d",
#'                "User4,c,c,p,c,d",
#'                "User5,h,c,c,p,p,c,p,p,p,i,p,o",
#'                "User6,i,h,c,c,p,p,c,p,c,d")
#' csf<-tempfile()
#' writeLines(clickstreams, csf)
#' cls<-readClickstreams(csf, header=TRUE)
#' print(cls)
#' 
#' @export readClickstreams
readClickstreams=function(file, sep=",", header=FALSE) {   
    count=max(count.fields(file, sep=sep))
    dat=read.table(file, sep=sep, header=F, fill=T, stringsAsFactors=F, col.names=c(1:count))
    
    if (header) {
        nams=dat[, 1]
        dat=dat[, -1]
    }
    
    ldat=split(dat, seq_len(nrow(dat)))
    ldat=llply(.data=ldat, .fun=function(x) x[x != ""])
    if (header)
        names(ldat)=nams
    else
        names(ldat)=seq(1, length(dat[,1]), 1)
    class(ldat)="Clickstreams"
    return(ldat)
}

.writeClickstream=function(name, clickstreamList, file, header, sep, quote) {  
    clickstream=clickstreamList[[name]]
    if (header)
        clickstream=c(name, clickstream)
    write.table(t(clickstream), file=file, sep=sep, append=T, row.names=F, col.names=F, quote=quote)
}



#' Writes a list of clickstreams to a file
#' 
#' Writes a list of clickstream to a csv-file.
#' 
#' 
#' @param clickstreamList The list of clickstreams to be written.
#' @param file The name of the file which the clickstreams are written to.
#' @param sep The character used to separate clicks (default is \dQuote{,}).
#' @param header A logical flag indicating whether the name of each clickstream
#' element should be used as first element.
#' @param quote A logical flag indicating whether each element of a clickstream
#' will be surrounded by double quotes (default is \code{TRUE}.
#' @author Michael Scholz \email{michael.scholz@@uni-passau.de}
#' @seealso \code{\link{readClickstreams}}, \code{\link{clusterClickstreams}}
#' @examples
#' 
#' clickstreams<-c("User1,h,c,c,p,c,h,c,p,p,c,p,p,o",
#'                "User2,i,c,i,c,c,c,d",
#'                "User3,h,i,c,i,c,p,c,c,p,c,c,i,d",
#'                "User4,c,c,p,c,d",
#'                "User5,h,c,c,p,p,c,p,p,p,i,p,o",
#'                "User6,i,h,c,c,p,p,c,p,c,d")
#' csf<-tempfile()
#' writeLines(clickstreams, csf)
#' cls<-readClickstreams(csf, header=TRUE)
#' clusters<-clusterClickstreams(cls, order=0, centers=2)
#' writeClickstreams(cls, file="clickstreams.csv", header=TRUE, sep=",")
#' 
#' @export writeClickstreams
writeClickstreams=function(clickstreamList, file, header=TRUE, sep=",", quote=TRUE) {
    l_ply(.data=names(clickstreamList), .fun=.writeClickstream, clickstreamList, file, header, sep, quote)
}

.randomClickstream=function(i, states, startProbabilities, transitionMatrix, meanLength) {   
    if (is.matrix(transitionMatrix)) {
        transitionMatrix=as.data.frame(transitionMatrix)
    }
    names(transitionMatrix)=states
    row.names(transitionMatrix)=states
    len=0
    while (len==0) {
        len=rpois(1, meanLength) 
    }
    start=runif(1, 0, 1)
    cs=cumsum(startProbabilities)
    previous=states[length(cs[cs<start])+1]
    chain=previous
    if (len>1) {
        for (i in seq(2, len, 1)) {
            nextState=runif(1, 0, 1)
            if (sum(transitionMatrix[,previous]) > 0) {
                cs=cumsum(transitionMatrix[,previous])     
                previous=states[length(cs[cs<nextState])+1]
                chain=c(chain, previous)
            } else {
                break
            }
        }
    }
    return(chain)
}





#' Generates a list of clickstreams
#' 
#' Generates a list of clickstreams by randomly walking through a given
#' transition matrix.
#' 
#' 
#' @param states Names of all possible states.
#' @param startProbabilities Start probabilities for all states.
#' @param transitionMatrix Matrix of transition probabilities.
#' @param meanLength Average length of the click streams.
#' @param n Number of click streams to be generated.
#' @return Returns a list of clickstreams.
#' @author Michael Scholz \email{michael.scholz@@uni-passau.de}
#' @seealso \code{\link{fitMarkovChain}}, \code{\link{readClickstreams}},
#' \code{\link{print.Clickstreams}}
#' @examples
#' 
#' # generate a simple list of click streams
#' states<-c("a", "b", "c")
#' startProbabilities<-c(0.2, 0.5, 0.3)
#' transitionMatrix<-matrix(c(0,0.4,0.6,0.3,0.1,0.6,0.2,0.8,0), nrow=3)
#' cls<-randomClickstreams(states, startProbabilities, transitionMatrix, meanLength=5, n=10)
#' print(cls)
#' 
#' @export randomClickstreams
randomClickstreams=function(states, startProbabilities, transitionMatrix, meanLength, n=100) {  
    s1=sum(aaply(.data=transitionMatrix, .margins=2, .fun=sum)==1)
    s2=sum(aaply(.data=transitionMatrix, .margins=2, .fun=sum)==0)
    if (!is.character(states)) 
        stop("states has to be a character vector")
    if (!is.numeric(startProbabilities))
        stop("startProbabilities has to be numeric")
    if (!is.matrix(transitionMatrix) && !is.data.frame(transitionMatrix))
        stop("transitionMatrix has to be a matrix or data frame")
    if (!is.numeric(meanLength))
        stop("meanLength has to be numeric")
    if (!is.numeric(n))
        stop("n has to be numeric")
    if (sum(startProbabilities)!=1)
        stop("The sum of startProbabilities has to be equal to 1")
    if (s1+s2<length(states)) {
        stop("The colums in transitionMatrix must sum up to 0 (absorbing states) or 1 (all other)")
    }
    clickstreamList=alply(.data=seq(1,n,1), .margins=1, .fun=.randomClickstream, states, startProbabilities, transitionMatrix, meanLength)
    class(clickstreamList)="Clickstreams"
    return(clickstreamList)
}

.printClick=function(click) {
    cat(click)
}

.printClickstream=function(pos, clickstreams) {
    clickstream=clickstreams[pos]
    cat(names(clickstreams)[pos])
    cat(": ")
    a_ply(.data=clickstream, .margins=1, .fun=.printClick)
    cat("\n")
}





#' Prints a list of clickstreams
#' 
#' Prints a list of clickstreams.
#' 
#' 
#' @param x A list of clickstreams.
#' @param ...  Ignored parameters.
#' @method print Clickstreams
#' @author Michael Scholz \email{michael.scholz@@uni-passau.de}
#' @seealso \code{\link{readClickstreams}}, \code{\link{randomClickstreams}}
#' @examples
#' 
#' clickstreams<-c("User1,h,c,c,p,c,h,c,p,p,c,p,p,o",
#'                "User2,i,c,i,c,c,c,d",
#'                "User3,h,i,c,i,c,p,c,c,p,c,c,i,d",
#'                "User4,c,c,p,c,d",
#'                "User5,h,c,c,p,p,c,p,p,p,i,p,o",
#'                "User6,i,h,c,c,p,p,c,p,c,d")
#' csf<-tempfile()
#' writeLines(clickstreams, csf)
#' cls<-readClickstreams(csf, header=TRUE)
#' print(cls)
#' 
#' @export print.Clickstreams
print.Clickstreams=function(x, ...) {
    cat("Clickstreams\n\n")
    a_ply(.data=seq(1, length(x), 1), .margins=1, .fun=.printClickstream, x)    
}



#' Prints a summary of a list of clickstreams
#' 
#' Returns a summary of a list of clickstreams.
#' 
#' 
#' @param object A list of clickstreams.
#' @param ...  Ignored parameters.
#' @method summary Clickstreams
#' @author Michael Scholz \email{michael.scholz@@uni-passau.de}
#' @seealso \code{\link{readClickstreams}}, \code{\link{randomClickstreams}}
#' @examples
#' 
#' clickstreams<-c("User1,h,c,c,p,c,h,c,p,p,c,p,p,o",
#'                "User2,i,c,i,c,c,c,d",
#'                "User3,h,i,c,i,c,p,c,c,p,c,c,i,d",
#'                "User4,c,c,p,c,d",
#'                "User5,h,c,c,p,p,c,p,p,p,i,p,o",
#'                "User6,i,h,c,c,p,p,c,p,c,d")
#' csf<-tempfile()
#' writeLines(clickstreams, csf)
#' cls<-readClickstreams(csf, header=TRUE)
#' summary(cls)
#' 
#' @export summary.Clickstreams
summary.Clickstreams=function(object, ...) {
    clicks=table(unlist(object))
    cat("Observations: ")
    cat(length(object))
    cat("\n\n")
    cat("Click Frequencies:")
    print(clicks)
}





#' Prints a ClickstreamClusters object
#' 
#' Prints a \code{ClickstreamClusters} object. A \code{ClickstreamClusters}
#' object represents the result of a cluster analysis on a list of clickstreams
#' (see \code{\link{clusterClickstreams}}).
#' 
#' 
#' @param x A \code{ClickstreamClusters} object returned by
#' \code{\link{clusterClickstreams}}.
#' @param ...  Ignored parameters.
#' @method print ClickstreamClusters
#' @author Michael Scholz \email{michael.scholz@@uni-passau.de}
#' @seealso \code{\link{clusterClickstreams}},
#' \code{\link{summary.ClickstreamClusters}}
#' @examples
#' 
#' clickstreams<-c("User1,h,c,c,p,c,h,c,p,p,c,p,p,o",
#'                "User2,i,c,i,c,c,c,d",
#'                "User3,h,i,c,i,c,p,c,c,p,c,c,i,d",
#'                "User4,c,c,p,c,d",
#'                "User5,h,c,c,p,p,c,p,p,p,i,p,o",
#'                "User6,i,h,c,c,p,p,c,p,c,d")
#' csf<-tempfile()
#' writeLines(clickstreams, csf)
#' cls<-readClickstreams(csf, header=TRUE)
#' clusters<-clusterClickstreams(cls, order=0, centers=2)
#' print(clusters)
#' 
#' @export print.ClickstreamClusters
print.ClickstreamClusters=function(x, ...) {
    print(x$clusters)
}



#' Prints a summary of a ClickstreamCluster object
#' 
#' Prints a summary of a \code{ClickstreamCluster} object. A
#' \code{ClickstreamClusters} object represents the result of a cluster
#' analysis on a list of clickstreams (see \code{\link{clusterClickstreams}}).
#' 
#' 
#' @param object A \code{ClickstreamClusters} object returned by
#' \code{\link{clusterClickstreams}}.
#' @param ...  Ignored parameters.
#' @method summary ClickstreamClusters
#' @author Michael Scholz \email{michael.scholz@@uni-passau.de}
#' @seealso \code{\link{clusterClickstreams}},
#' \code{\link{print.ClickstreamClusters}}
#' @examples
#' 
#' clickstreams<-c("User1,h,c,c,p,c,h,c,p,p,c,p,p,o",
#'                "User2,i,c,i,c,c,c,d",
#'                "User3,h,i,c,i,c,p,c,c,p,c,c,i,d",
#'                "User4,c,c,p,c,d",
#'                "User5,h,c,c,p,p,c,p,p,p,i,p,o",
#'                "User6,i,h,c,c,p,p,c,p,c,d")
#' csf<-tempfile()
#' writeLines(clickstreams, csf)
#' cls<-readClickstreams(csf, header=TRUE)
#' clusters<-clusterClickstreams(cls, order=0, centers=2)
#' summary(clusters)
#' 
#' @export summary.ClickstreamClusters
summary.ClickstreamClusters=function(object, ...) {
    cat("Centers:\n")
    for (i in 1:dim(object$centers)[1]) {
        cat("Cluster ", i, ":\n", sep="")
        transition=data.frame(matrix(object$centers[i,], ncol=length(object$states)))
        if (sqrt(dim(object$centers)[2]) == length(object$states)) {
            row.names(transition)=object$states
        }
        names(transition)=object$states
        print(transition)
        cat("\n\n")
    }
    cat("\n")
    cat("Total SS:", object$totss, "\n")
    cat("Within SS:", object$withinss, "\n")
    cat("Total Within SS:", object$tot.withinss, "\n")
    cat("Between SS:", object$betweenss, "\n")    
}





#' Predicts the cluster for a given Pattern object
#' 
#' Predicts the cluster for a given \code{Pattern} object. Potential clusters
#' need to be identified with the method \code{clusterClickstreams} before
#' predicting the cluster.
#' 
#' 
#' @param object A \code{ClickstreamClusters} object containing the potential
#' clusters. A \code{ClickstreamClusters} object represents the result of a
#' cluster analysis on a list of clickstreams (see
#' \code{\link{clusterClickstreams}}).
#' @param pattern The first clicks of a user as \code{Pattern} object.
#' @param ...  Ignored parameters.
#' @method predict ClickstreamClusters
#' @return Returns the index of the clusters to which the given \code{Pattern}
#' object most probably belongs to.
#' @author Michael Scholz \email{michael.scholz@@uni-passau.de}
#' @seealso \code{\link{clusterClickstreams}},
#' \code{\link{print.ClickstreamClusters}}
#' @examples
#' 
#' clickstreams<-c("User1,h,c,c,p,c,h,c,p,p,c,p,p,o",
#'                "User2,i,c,i,c,c,c,d",
#'                "User3,h,i,c,i,c,p,c,c,p,c,c,i,d",
#'                "User4,c,c,p,c,d",
#'                "User5,h,c,c,p,p,c,p,p,p,i,p,o",
#'                "User6,i,h,c,c,p,p,c,p,c,d")
#' csf<-tempfile()
#' writeLines(clickstreams, csf)
#' cls<-readClickstreams(csf, header=TRUE)
#' clusters<-clusterClickstreams(cls, order=0, centers=2)
#' pattern<-new("Pattern", sequence=c("h", "c"))
#' predict(clusters, pattern)
#' 
#' @export predict.ClickstreamClusters
predict.ClickstreamClusters=function(object, pattern, ...) {
    pos=aaply(.data=pattern@sequence, .margins=1, .fun=function(x) which(dimnames(object$centers)[[2]]==x))
    likelihood=aaply(.data=object$centers[,pos], .margins=1, .fun=prod)
    return(as.numeric(which(likelihood==max(likelihood))))
}
