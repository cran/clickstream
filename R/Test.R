#' evaluates the number of occurances of predicted next click vs total number of starting pattern occurances
#' in a given clickstream. The predicted next click can be from a markov chain of any order.
#' @export
#' @param mc a markovchain object (this should have been built from a set of training data) 
#' @param startPattern the starting pattern we want to predict next click on, and evaluate observed occurances in test data. 
#' @param testCLS clickstream object with test data 
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
#' 
#' csf <- tempfile()
#' writeLines(test, csf)
#' testCLS <- readClickstreams(csf, header = TRUE)
#' 
#' mc <- fitMarkovChain(trainingCLS, order = 1)
#' startPattern <- new("Pattern", sequence = c("c","c")) 
#' res <- mcEvaluate(mc, startPattern, testCLS)
#' res

mcEvaluate = function(mc, startPattern, testCLS){
  setClass(
    "Result",
    representation(
      totalclicks = "numeric",
      observed = "numeric",
      expected = "numeric",
      predictedNextClick = "character",
      patternSequence = "character",
      probability = "numeric"
    )
  )
  patternMatch <- function(p,s) {
    gg <- gregexpr(paste0("(?=",p,")"),s,perl=TRUE)[[1]]
    if (length(gg)==1 && gg==-1) 0 else length(gg)
  }
  pred <- predict(mc, startPattern)
  expec <- pred@probability
  predicted <- as.character(pred@sequence)
  vec_totalPages <- NULL
  vec_observed <- NULL
  for (i in testCLS){
    
    clicks <- paste(i, collapse = ",")
    pattern <- paste(startPattern@sequence, collapse = ",")
    pattern <- as.character(pattern)
    patternchar <- paste(pattern, ",", sep = "")
    totalPages <- patternMatch(patternchar, clicks)
    vec_totalPages <- append(vec_totalPages, totalPages)
    expectedSeq <- paste(pattern, predicted, sep = ",")
    obs <- patternMatch(expectedSeq, clicks)
    vec_observed <- append(vec_observed, obs)
  }
  TotalExpectded <- sum(vec_totalPages) * expec
  result = new(
    "Result", expected = TotalExpectded, observed = sum(vec_observed),
    predictedNextClick = predicted, patternSequence = pattern, totalclicks = sum(vec_totalPages),
    probability = expec
  )
  return(result)
}

#' evaluates all next page clicks in a clickstream training data set 
#' against the test data on the basis of a 1st Order Markov chain
#' and produces a report of observed and expected values in a matrix
#'
#' @export
#' @param mc a markovchain object that corresponds to a list of clusters
#' @param trainingCLS clickstream object with training data (this should be the data used to build the markov chain object)
#' @param testCLS clickstream object with test data 
#' @param mc the markov chain against which to compare the clickstream data
#' @examples
#' training <- c("User1,h,c,c,p,c,h,c,p,p,c,p,p,o",
#'               "User2,i,c,i,c,c,c,d",
#'               "User3,h,i,c,i,c,p,c,c,p,c,c,i,d",
#'               "User4,c,c,p,c,d")
#' 
#' test <- c("User1,h,h,h,h,c,c,p,p,h,c,p,p,c,p,p,o",
#'           "User2,i,c,i,c,c,c,d",
#'           "User4,c,c,c,c,d,c,c,c,c,c,c,c,c,c,c,c,c,c,c")
#' 
#' csf <- tempfile()
#' writeLines(training, csf)
#' trainingCLS <- readClickstreams(csf, header = TRUE)
#' 
#' csf <- tempfile()
#' writeLines(test, csf)
#' testCLS <- readClickstreams(csf, header = TRUE)
#' 
#' mc <- fitMarkovChain(trainingCLS, order = 1)
#' res <- mcEvaluateAll(mc, testCLS, trainingCLS)
#' write.csv(res, file = "results.csv")


mcEvaluateAll = function(mc, testCLS, trainingCLS, ord=1){
  results <- data.frame( "totalclicks" = character(), "observed" = character(), "expected" = character(), "predictedNextClick" = character(), "patternSequence" = character(), "probability" = character(),  stringsAsFactors=FALSE)
  vec<-unlist(trainingCLS)
  dedupe <- vec[which(!duplicated(vec))]
  for (d in dedupe){
    if(d !="Defer"){
      value <- d[[1]]
      startPattern <- new("Pattern", sequence = c(value)) 
      res <- mcEvaluate(mc, startPattern, testCLS)
      if (res@totalclicks != 0 && res@expected > 0){
        vec_results <- c(res@totalclicks, res@observed, res@expected, res@patternSequence, res@predictedNextClick,res@probability)
        results[nrow(results) + 1, ] <- c(res@totalclicks, res@observed, res@expected, res@predictedNextClick, res@patternSequence, res@probability)
      }
    }
  }
  return(results)
}


#' evaluates all next page clicks in a clickstream training data set 
#' against the test data on the basis of a set of pre-computed 1st Order Markov chains
#' and corresponding clusters, and produces a report of observed and expected values in a matrix
#'
#' @export
#' @param markovchains a list of pre-computed markovchain objects that correspond to a list of clusters
#' @param clusters a clist of clusters
#' @param trainingCLS clickstream object with training data (this should be the data used to build the markov chains objects)
#' @param testCLS clickstream object with test data 
#' @param mc the markov chain against which to compare the clickstream data
#' @examples
#' training <- c("User1,h,c,c,p,c,h,c,p,p,c,p,p,o",
#'               "User2,i,c,i,c,c,c,d",
#'               "User3,h,i,c,i,c,p,c,c,p,c,c,i,d",
#'               "User4,c,c,p,c,d")
#' 
#' test <- c("User1,h,c,c,p,p,h,c,p,p,c,p,p,o",
#'           "User2,i,c,i,c,c,c,d",
#'           "User4,c,c,c,c,d")
#' 
#' csf <- tempfile()
#' writeLines(training, csf)
#' trainingCLS <- readClickstreams(csf, header = TRUE)
#' 
#' csf <- tempfile()
#' writeLines(test, csf)
#' testCLS <- readClickstreams(csf, header = TRUE)
#' 
#' clusters <- clusterClickstreams(trainingCLS, centers = 2)
#' markovchains <- fitMarkovChains(clusters, order = 1)
#' res <- mcEvaluateAllClusters(markovchains, clusters, testCLS, trainingCLS)
#' write.csv(res, file = "results.csv")

mcEvaluateAllClusters = function(markovchains, clusters, testCLS, trainingCLS){
  results <- data.frame( "totalclicks" = character(), "observed" = character(), "expected" = character(), "predictedNextClick" = character(), "patternSequence" = character(), "probability" = character(),  stringsAsFactors=FALSE)
  vec<-unlist(trainingCLS)
  dedupe <- vec[which(!duplicated(vec))]
  for (d in dedupe){
    if(d !="Defer"){
      value <- d[[1]]
      startPattern <- new("Pattern", sequence = c(value)) 
      mc <- getOptimalMarkovChain(startPattern,markovchains,clusters)
      res <- mcEvaluate(mc, startPattern, testCLS)
      if (res@totalclicks != 0){
        vec_results <- c(res@totalclicks, res@observed, res@expected, res@patternSequence, res@predictedNextClick,res@probability)
        results[nrow(results) + 1, ] <- c(res@totalclicks, res@observed, res@expected, res@predictedNextClick, res@patternSequence, res@probability)
      }
    }
  }
  return(results)
}


getSequenceTotal = function(startPattern, nextSequence, testCLS){
  setClass(
    "Result",
    representation(
      totalclicks = "numeric",
      observed = "numeric",
      nextSequence = "character",
      startPattern = "character"
    )
  )
  patternMatch <- function(p,s) {
    gg <- gregexpr(paste0("(?=",p,")"),s,perl=TRUE)[[1]]
    if (length(gg)==1 && gg==-1) 0 else length(gg)
  }
  nextseq <- nextSequence
  vec_totalPages <- NULL
  vec_observed <- NULL
  for (i in testCLS){
    clicks <- paste(i, collapse = ",")
    clicks <- paste(clicks, ",", sep = "")
    pattern <- paste(startPattern@sequence, collapse = ",")
    pattern <- as.character(pattern)
    patternchar <- paste(pattern, ",", sep = "")
    totalPages <- patternMatch(patternchar, clicks)    
    vec_totalPages <- append(vec_totalPages, totalPages)
    expectedSeq <- paste(pattern, nextseq, sep = ",")
    obs <- patternMatch(expectedSeq, clicks)
    vec_observed <- append(vec_observed, obs)
  }
  result = new(
    "Result",  observed = sum(vec_observed),
    nextSequence = nextseq, startPattern = pattern, totalclicks = sum(vec_totalPages)
  )
  return(result)
}

getSequenceMatrix = function(mc, testCLS){
  empiricalMatrix <- mc
  cols <- colnames(empiricalMatrix)
  rows <- rownames(empiricalMatrix)
  for (c in rows){
    for(r in cols){
      value <- empiricalMatrix[[c,r]]
      startPattern <- new("Pattern", sequence = c(c))
      nextSequence <- r
      value <- getSequenceTotal(startPattern, nextSequence, testCLS)
      empiricalMatrix[[c,r]] <- value@observed
    }
  }
  return(empiricalMatrix)
}


#' calculates the Chi-Square Statistic, p-value, and degrees of freedom, for a transition matrix compared with observed state changes. 
#' 
#'
#' @export
#' @param cls The clickstream
#' @param mc the markov chain against which to compare the clickstream data
#' @examples
#' clickstreams <- c("User1,h,c,c,p,c,h,c,p,p,c,p,p,o",
#'                  "User2,i,c,i,c,c,c,d",
#'                  "User3,h,i,c,i,c,p,c,c,p,c,c,i,d",
#'                  "User4,c,c,p,c,d")
#'
#' csf <- tempfile()
#' writeLines(clickstreams, csf)
#' cls <- readClickstreams(csf, header = TRUE)
#' mc <- fitMarkovChain(cls)
#' chiSquareTest(cls,mc, verbose = TRUE)

chiSquareTest <- function(cls, mc, printValues = TRUE) {
  object <- as.data.frame(t(mc@transitions[[1]]))
  data <- cls
  data <- getSequenceMatrix(object, data)
  data <- data[match(rownames(data),names(object)),]
  data <- data[,match(colnames(data),names(object))]
  cols <- colSums(data)
  statistic <- 0
  for (i in 1:length(object)) {
    for (j in 1:length(object)) {
      if (data[i, j]>0&object[i, j]>0) statistic <- statistic + data[i, j]*log(data[i, j]/(cols[i]*object[i, j]))
    }
  }
  statistic <- statistic * 2
  null_elements <- sum(object == 0)
  dof <- length(object) * (length(object) - 1) - null_elements
  p.value <- 1 - pchisq(q = statistic,df = dof)
  if(printValues==TRUE){
    cat("Test data:\n");print(data);cat("Transition matrix: \n");print(object);print("...")
    cat("Chi-Square Statistic:",statistic,"\n")
    cat("degrees of freedom:",dof,"\n")
    cat("The p-value:",p.value,"\n")  
  }
  out <- list(statistic = statistic, dof = dof,pvalue = p.value)
  return(out)
}
