#' Analyzes Clickstreams Based on Markov Chains
#' 
#' This package allows modeling clickstreams with Markov chains. It supports to
#' model clickstreams as zero-order, first-order or higher-order Markov chains.
#' 
#' \tabular{ll}{ Package: \tab clickstream\cr Type: \tab Package\cr Version:
#' \tab 1.3.3\cr Date: \tab 2023-09-27\cr License: \tab GPL-2\cr Depends: \tab
#' R (>= 3.0), methods\cr }
#' 
#' @name clickstream-package
#' @aliases clickstream-package clickstream
#' @docType package
#' @author Michael Scholz \email{michael.scholz@@th-deg.de}
#' @author Theo van Kraay \email{theo.vankraay@@hotmail.com}
#' @references Scholz, M. (2016) R Package clickstream: Analyzing Clickstream Data with Markov Chains, \emph{Journal of Statistical Software}, \bold{74}, 4, pages 1--17 .
#' 
#' Ching, W.-K.and Huang, X. and Ng, M.K. and Siu, T.-K. (2013) \emph{Markov Chains -- Models, Algorithms and Applications}, 2nd edition, New York: Springer-Verlag.
#' @import methods Rsolnp arules plyr linprog ggplot2 ClickClust parallel
#' @importFrom stats runif rpois kmeans cor pchisq
#' @importFrom utils read.table write.table count.fields
#' @importFrom igraph E graph.adjacency
#' @importFrom reshape2 melt
#' @importFrom data.table data.table dcast.data.table as.data.table
#' @importFrom MASS ginv
#' @keywords click stream Markov chain
#' @examples
#' 
#' # fitting a simple Markov chain and predicting the next click
#' clickstreams <- c("User1,h,c,c,p,c,h,c,p,p,c,p,p,o",
#'                "User2,i,c,i,c,c,c,d",
#'                "User3,h,i,c,i,c,p,c,c,p,c,c,i,d",
#'                "User4,c,c,p,c,d",
#'                "User5,h,c,c,p,p,c,p,p,p,i,p,o",
#'                "User6,i,h,c,c,p,p,c,p,c,d")
#' cls <- as.clickstreams(clickstreams, header = TRUE)
#' mc <- fitMarkovChain(cls)
#' startPattern <- new("Pattern", sequence = c("h", "c"))
#' predict(mc, startPattern)
#' plot(mc)
#' 
NULL
