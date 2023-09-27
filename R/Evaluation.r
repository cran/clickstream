#' Class \code{EvaluationResult}
#'
#' @name EvaluationResult-class
#' @aliases EvaluationResult-class
#' @docType class
#' @section Objects from the Class: Objects can be created by calls of the form
#' \code{new("EvaluationResult", ...)}. This S4 class describes \code{EvaluationResult} objects.
#' @author Michael Scholz \email{michael.scholz@@th-deg.de}
#' @seealso \code{\link{mcEvaluate}}
#' @keywords classes
#' @examples
#'
#' # show EvaluationResult definition
#' showClass("EvaluationResult")
#'
#' @export
setClass(
    "EvaluationResult",
    representation(
        totalclicks = "numeric",
        observed = "numeric",
        expected = "numeric",
        residual = "numeric",
        residualSquared = "numeric",
        component = "numeric",
        predictedNextClick = "character",
        patternSequence = "character",
        probability = "numeric"
    )
)

#' Shows an \code{EvaluationResult} object
#'
#' @export
#' @docType methods
#' @rdname EvaluationResult-method 
#' @param object An instance of the \code{EvaluationResult}-class
#' @section Methods: \describe{
#' \item{list("signature(object = \"EvaluationResult\")")}{ Shows an \code{EvaluationResult} object. } }
#' @author Michael Scholz \email{michael.scholz@@th-deg.de}
#' @keywords methods
setMethod("show", "EvaluationResult",
          function(object) {
              cat("\n\n")
            }
)


.patternMatch <- function(p,s) {
    gg <- gregexpr(paste0("(?=", p, ")"), s, perl=TRUE)[[1]]
    if (length(gg)==1 && gg==-1) 
        return(0) 
    else 
        return(length(gg))
}


#' Evaluates the number of occurrences of predicted next clicks
#' 
#' @description Evaluates the number of occurrences of predicted next clicks vs. total number of starting pattern occurrences
#' in a given clickstream. The predicted next click can be a markov chain of any order.
#' @param mc a markovchain object (this should have been built from a set of training data) 
#' @param startPattern the starting pattern we want to predict next click on, and evaluate observed occurrences in test data. 
#' @param testCLS clickstream object with test data 
#' @author Theo van Kraay \email{theo.vankraay@@hotmail.com}
#' @examples
#' training <- c("User1,h,c,c,p,c,h,c,p,p,c,p,p,o",
#'               "User2,i,c,i,c,c,c,d",
#'               "User3,h,i,c,i,c,p,c,c,p,c,c,i,d",
#'               "User4,c,c,p,c,d")
#' 
#' test <- c("User1,h,h,h,h,c,c,p,p,h,c,p,p,c,p,p,o",
#'           "User2,i,c,i,c,c,c,d",
#'           "User4,c,c,c,c,d,c,c,c,c")
#' 
#' csf <- tempfile()
#' writeLines(training, csf)
#' trainingCLS <- readClickstreams(csf, header = TRUE)
#' unlink(csf)
#' 
#' csf <- tempfile()
#' writeLines(test, csf)
#' testCLS <- readClickstreams(csf, header = TRUE)
#' unlink(csf)
#' 
#' mc <- fitMarkovChain(trainingCLS, order = 1)
#' startPattern <- new("Pattern", sequence = c("c","c")) 
#' res <- mcEvaluate(mc, startPattern, testCLS)
#' res
#' @export mcEvaluate
mcEvaluate = function(mc, startPattern, testCLS) {
    pred <- predict(mc, startPattern)
    expec <- pred@probability
    predicted <- as.character(pred@sequence)
    vec_totalPages <- vector()
    vec_observed <- vector()
    for (i in testCLS){
        
        clicks <- paste(i, collapse = ",")
        pattern <- paste(startPattern@sequence, collapse = ",")
        pattern <- as.character(pattern)
        patternchar <- paste(pattern, ",", sep = "")
        totalPages <- .patternMatch(patternchar, clicks)
        vec_totalPages <- c(vec_totalPages, totalPages)
        expectedSeq <- paste(pattern, predicted, sep = ",")
        obs <- .patternMatch(expectedSeq, clicks)
        vec_observed <- c(vec_observed, obs)
    }
    TotalExpectded <- sum(vec_totalPages) * expec
    observed = sum(vec_observed)
    expected = TotalExpectded
    residual = observed-expected
    residualSquared = residual^2
    component = residualSquared/expected
    result = new(
        "EvaluationResult", expected = expected, observed = observed, residual = residual,
        residualSquared = residualSquared, component = component,
        predictedNextClick = predicted, patternSequence = pattern, totalclicks = sum(vec_totalPages),
        probability = expec
    )
    return(result)
}

.getNOrderPatterns = function(trainingCLS, order) {
    secondOrder = function(patternlist = NULL) {
        if (is.null(patternlist)){
            vec<-unlist(trainingCLS)
            dedupe <- vec[which(!duplicated(vec))]
        }
        else{
            dedupe <- patternlist
        }
        listPatterns <- list()
        OuterlistPatterns <- list()
        for (a in dedupe){
            value <- a[[1]]
            vecPatternsList <- list()
            for(b in dedupe){
                vecPattern <- c(value, b)
                vecPatternsList <- c(vecPatternsList, list(vecPattern))
            }
            OuterlistPatterns <- c(OuterlistPatterns, list(vecPatternsList))
        }
        OuterlistPatterns <- unlist(OuterlistPatterns, recursive = FALSE)
    }
    for(i in 2:order){
        if(i>2){
            patternlist <- secondOrder(patternlist = patternlist)
        }
        else{
            patternlist <- secondOrder()
        }
    }
    return(patternlist)
}

#' Evaluates all next page clicks in a clickstream training data set against a test data
#' 
#' @description Evaluates all next page clicks in a clickstream training data set against a test data. Handles higher order by cycling through every possible pattern permutation. Produces a report of observed and expected values in a matrix.
#' @param mc A markovchain object that corresponds to a list of clusters.
#' @param trainingCLS Clickstream object with training data (this should be the data used to build the markov chain object).
#' @param testCLS Clickstream object with test data.
#' @param includeChiSquare Should the result include the chi-square value?
#' @param returnChiSquareOnly Should the result only consist of the chi-square value?
#' @seealso \code{\link[=mcEvaluate]{mcEvaluate}}
#' @author Theo van Kraay \email{theo.vankraay@@hotmail.com}
#' @examples
#' training <- c("User1,h,c,c,p,c,h,c,p,p,c,p,p",
#'               "User2,i,c,i,c,c,c,d")
#' 
#' test <- c("User1,h,c,c,p,c,h,c,d,p,c,d,p",
#'              "User2,i,c,i,p,c,c,d")
#' 
#' csf <- tempfile()
#' writeLines(training, csf)
#' trainingCLS <- readClickstreams(csf, header = TRUE)
#' unlink(csf)
#' 
#' csf <- tempfile()
#' writeLines(test, csf)
#' testCLS <- readClickstreams(csf, header = TRUE)
#' unlink(csf)
#' 
#' mc <- fitMarkovChain(trainingCLS, order = 2)
#' mcEvaluateAll(mc, trainingCLS, testCLS)
#' @export mcEvaluateAll
mcEvaluateAll = function(mc, trainingCLS, testCLS, includeChiSquare = TRUE, returnChiSquareOnly = FALSE){
    results <- data.frame( "totalclicks" = character(), "observed" = character(), "expected" = character(), "residual" = character(), "residualSquared" = character(), "component" = numeric(), "predictedNextClick" = character(), "patternSequence" = character(), "probability" = character(),  stringsAsFactors=FALSE)
    if(mc@order==1){
        vec<-unlist(trainingCLS)
        dedupe <- vec[which(!duplicated(vec))]
    }
    else{
        dedupe <- .getNOrderPatterns(trainingCLS, order = mc@order)
    }
    for (d in dedupe){
        if(mc@order==1){
            value <- d[[1]]
        }
        else{
            value <- d
        }
        startPattern <- new("Pattern", sequence = c(value)) 
        res <- mcEvaluate(mc, startPattern, testCLS)
        if (res@totalclicks != 0 && res@expected > 0){
            vec_results <- c(res@totalclicks, res@observed, res@expected, res@residual, res@residualSquared, res@component, res@patternSequence, res@predictedNextClick,res@probability)
            results[nrow(results) + 1, ] <- c(res@totalclicks, res@observed, res@expected, res@residual, res@residualSquared, res@component, res@predictedNextClick, res@patternSequence, res@probability)
        }
    }   
    ChiSquare <- sum(as.numeric(results$component))
    if (includeChiSquare == TRUE){
        results[nrow(results) + 1, ] <- c(0, 0, 0, 0, "variance:", ChiSquare, 0, 0, 0)
    }
    if (returnChiSquareOnly == TRUE){
        results <- ChiSquare
    }
    return(results)
}


#' Evaluates all next page clicks in a clickstream training data set against a test data
#' 
#' @description Evaluates all next page clicks in a clickstream training data set against a test data on the basis of a set of pre-computed Markov chains and corresponding clusters. Handles higher order by cycling through every possible pattern permutation. Produces and produces a report of observed and expected values in a matrix
#' @param markovchains A list of MarkovChain-objects.
#' @param clusters The list of clusters.
#' @param trainingCLS Clickstream object with training data (this should be the data used to build the markov chain object).
#' @param testCLS Clickstream object with test data.
#' @param includeChiSquare Should the result include the chi-square value?
#' @param returnChiSquareOnly Should the result only consist of the chi-square value?
#' @seealso \code{\link[=mcEvaluateAll]{mcEvaluateAll}}
#' @author Theo van Kraay \email{theo.vankraay@@hotmail.com}
#' @examples
#' training <- c("User1,h,c,c,p,c,h,c,h,o,p,p,c,p,p,o",
#'               "User2,i,c,i,c,c,c,o,o,o,i,d",
#'               "User3,h,i,c,i,c,o,i,p,c,c,p,c,c,i,d",
#'               "User4,c,c,p,c,d,o,i,h,o,o")
#' 
#' test <- c("User1,h,c,c,p,p,h,o,i,c,p,p,c,p,p,o",
#'           "User2,i,c,i,c,c,c,d",
#'           "User4,c,c,c,c,d")
#' 
#' csf <- tempfile()
#' writeLines(training, csf)
#' trainingCLS <- readClickstreams(csf, header = TRUE)
#' unlink(csf)
#' 
#' csf <- tempfile()
#' writeLines(test, csf)
#' testCLS <- readClickstreams(csf, header = TRUE)
#' unlink(csf)
#' 
#' clusters <- clusterClickstreams(trainingCLS, centers = 2, order = 1)
#' markovchains <- fitMarkovChains(clusters, order = 2)
#' mcEvaluateAllClusters(markovchains, clusters, testCLS, trainingCLS)
#' @export mcEvaluateAllClusters
mcEvaluateAllClusters = function(markovchains, clusters, testCLS, trainingCLS, includeChiSquare = TRUE, returnChiSquareOnly = FALSE){
    results <- data.frame( "totalclicks" = character(), "observed" = character(), "expected" = character(), "residual" = character(), "residualSquared" = character(), "component" = numeric(), "predictedNextClick" = character(), "patternSequence" = character(), "probability" = character(),  stringsAsFactors=FALSE)
    order = markovchains[[1]]@order
    if (order==1){
        vec<-unlist(trainingCLS)
        dedupe <- vec[which(!duplicated(vec))]
    } else {
        dedupe <- .getNOrderPatterns(trainingCLS, order = order)
    }   
    for (d in dedupe) {
        if(order==1){
            value <- d[[1]]
        } else{
            value <- d
        }
        startPattern <- new("Pattern", sequence = c(value)) 
        mc <- getOptimalMarkovChain(startPattern, markovchains, clusters)
        res <- mcEvaluate(mc, startPattern, testCLS)
        if (res@totalclicks != 0){
            vec_results <- c(res@totalclicks, res@observed, res@expected, res@residual, res@residualSquared, res@component, res@patternSequence, res@predictedNextClick,res@probability)
            results[nrow(results) + 1, ] <- c(res@totalclicks, res@observed, res@expected, res@residual, res@residualSquared, res@component, res@predictedNextClick, res@patternSequence, res@probability)
        }
    }
    ChiSquare <- sum(as.numeric(results$component))
    if (includeChiSquare == TRUE){
        results[nrow(results) + 1, ] <- c(0, 0, 0, 0, "variance:", ChiSquare, 0, 0, 0)
    }
    if (returnChiSquareOnly == TRUE){
        results <- ChiSquare
    }
    return(results)
}


.getSequenceTotal = function(startPattern, nextSequence, testCLS){
    nextseq <- nextSequence
    vec_totalPages <- NULL
    vec_observed <- NULL
    for (i in testCLS){
        clicks <- paste(i, collapse = ",")
        clicks <- paste(clicks, ",", sep = "")
        pattern <- paste(startPattern@sequence, collapse = ",")
        pattern <- as.character(pattern)
        patternchar <- paste(pattern, ",", sep = "")
        totalPages <- .patternMatch(patternchar, clicks)    
        vec_totalPages <- c(vec_totalPages, totalPages)
        expectedSeq <- paste(pattern, nextseq, sep = ",")
        obs <- .patternMatch(expectedSeq, clicks)
        vec_observed <- c(vec_observed, obs)
    }
    result = list(observed = sum(vec_observed), nextSequence = nextseq, 
                  startPattern = pattern, totalclicks = sum(vec_totalPages))
    return(result)
}

.getSequenceMatrix = function(mc, testCLS){
    empiricalMatrix <- mc
    cols <- colnames(empiricalMatrix)
    rows <- rownames(empiricalMatrix)
    for (c in rows){
        for(r in cols){
            value <- empiricalMatrix[[c,r]]
            startPattern <- new("Pattern", sequence = c(c))
            nextSequence <- r
            value <- .getSequenceTotal(startPattern, nextSequence, testCLS)
            empiricalMatrix[[c,r]] <- value$observed
        }
    }
    return(empiricalMatrix)
}


#' Calculates the chi-square statistic
#' 
#' @description Calculates the chi-Square statistic, p-value, and degrees of freedom, for the first-order transition matrix of a \code{MarkovChain} object compared with observed state changes. 
#' @param cls The clickstream object.
#' @param mc The Markov chain against which to compare the clickstream data. Please note that the first-order transition matrix is used for performing the chi-square test.
#' @author Theo van Kraay \email{theo.vankraay@@hotmail.com}
#' @examples
#' clickstreams <- c("User1,h,c,c,p,c,h,c,p,p,c,p,p,o",
#'                  "User2,i,c,i,c,c,c,d",
#'                  "User3,h,i,c,i,c,p,c,c,p,c,c,i,d",
#'                  "User4,c,c,p,c,d")
#'
#' csf <- tempfile()
#' writeLines(clickstreams, csf)
#' cls <- readClickstreams(csf, header = TRUE)
#' unlink(csf)
#' 
#' mc <- fitMarkovChain(cls)
#' chiSquareTest(cls, mc)
#' @export chiSquareTest
chiSquareTest <- function(cls, mc) {
    object <- as.data.frame(t(mc@transitions[[1]]))
    data <- cls
    data <- .getSequenceMatrix(object, data)
    data <- data[match(rownames(data),names(object)),]
    data <- data[,match(colnames(data),names(object))]
    cols <- colSums(data)
    statistic <- 0
    for (i in 1:length(object)) {
        for (j in 1:length(object)) {
            if (data[i, j]>0 && object[i, j]>0) statistic <- statistic + data[i, j]*log(data[i, j]/(cols[i]*object[i, j]))
        }
    }
    statistic <- statistic * 2
    nullElements <- sum(object == 0)
    dof <- length(object) * (length(object) - 1) - nullElements
    p.value <- 1 - pchisq(q = statistic, df = dof)
    out <- list(statistic = statistic, dof = dof, pvalue = p.value)
    return(out)
}
