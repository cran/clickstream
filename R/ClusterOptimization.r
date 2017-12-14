#' Generates a list of markov chains from a given set of clusters
#'
#' @description The purpose of this function is to generate pre-computed markov chain objects from clusters of clickstreams.
#' @param clusters The clusters from which to generate markov chain objects.
#' @param order The order for the markov chain.
#' @author Theo van Kraay \email{theo.vankraay@@hotmail.com}
#' @examples 
#' 
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
#' unlink(csf)
#' 
#' csf <- tempfile()
#' writeLines(test, csf)
#' testCLS <- readClickstreams(csf, header = TRUE)
#' unlink(csf)
#' 
#' clusters <- clusterClickstreams(trainingCLS, centers = 2)
#' markovchains <- fitMarkovChains(clusters, order = 1)
#' @export fitMarkovChains
fitMarkovChains = function(clusters, order=1) {
    markovchains <- vector()
    for (i in clusters[[1]]){
        mc <- fitMarkovChain(i, order = order) 
        markovchains <- c(markovchains, mc)
    }
    return(markovchains)
}


#' Generates the optimal markov chains from a list of markov chains and corresponding clusters
#'
#' @export
#' @description The purpose of this function is to predict from a pattern using pre-computed markov chains and corresponding clusters. The markov chain corresponding with the cluster that is the best fit to the prediction value is used.
#' @param startPattern The pattern object to be used.
#' @param markovchains The pre-computed markov chains generated from a set of clusters.
#' @param clusters The corresponding clusters (should be in the corresponding order as the markov chains).
#' @author Theo van Kraay \email{theo.vankraay@@hotmail.com}
#' @examples 
#' 
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
#' unlink(csf)
#' 
#' csf <- tempfile()
#' writeLines(test, csf)
#' testCLS <- readClickstreams(csf, header = TRUE)
#' unlink(csf)
#' 
#' clusters <- clusterClickstreams(trainingCLS, centers = 2)
#' markovchains <- fitMarkovChains(clusters, order = 1)
#' startPattern <- new("Pattern", sequence = c("c")) 
#' mc <- getOptimalMarkovChain(startPattern, markovchains, clusters)
#' predict(mc, startPattern)
getOptimalMarkovChain =function(startPattern, markovchains, clusters) {
    markovchainIndex <- predict(clusters, startPattern)
    optimalPreComputedChain <- markovchains[[markovchainIndex]]
    return(optimalPreComputedChain)
}

#' Generates an optimal set of clusters for a clickstream based on certain constraints
#'
#' @description This is an experimental function for a consensus clustering algorithm based on targeting a range of average next state probabilities derived when fitting each cluster to a markov chain. 
#' @param trainingCLS Clickstream object with training data (this should be the data used to build the markov chain object).
#' @param testCLS Clickstream object with test data.
#' @param maxIterations Number of times to iterate (repeat) through the k-means clustering.
#' @param optimalProbMean The target average probability of each next page click prediction in a 1st order markov chain.
#' @param range The range above the optimal probability to target. 
#' @param centresMin The minimum cluster centres to evaluate.
#' @param clusterCentresRange the additional cluster centres to evaluate.
#' @param order The order for markov chains that will be used to evaluate each cluster.
#' @param takeHighest determines whether to default to the highest mean next click probability, or error if the target is not reached after the given number of k-means iterations. 
#' @param verbose Should this function report extra information on progress?
#' @author Theo van Kraay \email{theo.vankraay@@hotmail.com}
#' @examples
#' training <- c("User1,h,c,c,p,c,h,c,p,p,c,p,p,o",
#'               "User2,i,c,i,c,c,c,d",
#'               "User3,h,i,c,i,c,p,c,c,p,c,c,i,d",
#'               "User4,h,c,c,p,p,c,p,p,p,i,p,o",
#'               "User5,i,h,c,c,p,p,c,p,c,d",
#'               "User6,i,h,c,c,p,p,c,p,c,o",
#'               "User7,i,h,c,c,p,p,c,p,c,d",
#'               "User8,i,h,c,c,p,p,c,p,c,d,o")
#' 
#' test <- c(
#'     "User1,h,c,c,p,c,h,c,p,p,c,p,p,o",
#'     "User2,i,c,i,c,c,c,d",
#'     "User3,h,i,c,i,c,p,c,c,p,c,c,i,d"
#' )
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
#' clusters <- getConsensusClusters(trainingCLS, testCLS, maxIterations=5, 
#'                          optimalProbMean=0.40, range = 0.70, centresMin = 2, 
#'                          clusterCentresRange = 0, order = 1, takeHighest = FALSE, 
#'                          verbose = FALSE)
#' markovchains <- fitMarkovChains(clusters)
#' startPattern <- new("Pattern", sequence = c("i", "h", "c", "p"))
#' mc <- getOptimalMarkovChain(startPattern, markovchains, clusters)
#' predict(mc, startPattern)
#' @export getConsensusClusters
getConsensusClusters = function(trainingCLS, testCLS, maxIterations=10, optimalProbMean=0.50, range=0.30, 
                                centresMin=2, clusterCentresRange=0, order=1, takeHighest=FALSE, verbose=FALSE){
    cls <- trainingCLS
    vec <- unlist(cls)
    dedupe <- vec[which(!duplicated(vec))]
    centresMax <- centresMin + clusterCentresRange
    listOfClusters <- list()
    clusterCentres <- centresMin:centresMax
    iterations <- 1:maxIterations
    vectorOfAllProbsMeans <- vector()
    limit <- optimalProbMean + range
    for (i in iterations) {
        for (cc in clusterCentres) {
            clusters <- clusterClickstreams(cls, centers = cc) 
            markovchains <- fitMarkovChains(clusters, order = order)
            absorbingStates <- unique(unlist(lapply(markovchains, FUN = function(x) return(x@absorbingStates))))
            vectorOfProbs <- vector()
            if (verbose)
                cat("Starting next page probability aggregation....\n")
            for (d in dedupe){
                if(!d %in% absorbingStates){
                    value <- d[[1]]
                    startPattern <- new("Pattern", sequence = c(value)) 
                    mc <- getOptimalMarkovChain(startPattern, markovchains, clusters)
                    prob <- predict(mc, startPattern)
                    vectorOfProbs <- c(vectorOfProbs, prob@probability)
                }
            }
            vectorOfAllProbsMeans <- c(vectorOfAllProbsMeans, mean(vectorOfProbs))
            listOfClusters <- c(listOfClusters, list(clusters))
        }
        candidates <- which(vectorOfAllProbsMeans > optimalProbMean & vectorOfAllProbsMeans < limit)
        if (verbose)
            cat("Candidates are:", candidates, "\n")
    }
    if (takeHighest != TRUE){
        if (length(candidates) > 0){
            #get the candidate clusters into a vector
            candidateClusters <- list()
            for (i in candidates){
                clusters <- listOfClusters[[i]]
                candidateClusters <- c(candidateClusters, list(clusters))
            }
            if (verbose)
                cat("Evaluating candidates....\n")
            vec_variances <- vector()
            for(c in candidateClusters){
                markovchains <- fitMarkovChains(c) 
                variance <- mcEvaluateAllClusters(markovchains, c, testCLS, trainingCLS, returnChiSquareOnly = TRUE)
                if (verbose)
                    cat("Variance is....", variance, "\n")
                vec_variances <- c(vec_variances, variance)
            }
            if (verbose)
                cat("Vector of variances is:", vec_variances, "\n")
            winner <- which.min(vec_variances)
            if (verbose)
                cat("Winner is:", winner, "\n")
            return(candidateClusters[[winner]])
        }
        else{
            stop(("Target range was not reached with the given number of iterations."))
        }
    }
    else{
        if (length(candidates) == 0){
            warning("Target prediction accuracy was not reached with the given number of iterations. Taking highest probability mean.")
        }
        candidates <- which(vectorOfAllProbsMeans == max(vectorOfAllProbsMeans))
        return(listOfClusters[[candidates]])
    }
}

.getParallelClusterSets = function(trainingCLS, maxIterations, centres, cores){
    mkWorker <- function(centres) {
        fitMarkovChains =function(clusters, order=1) {
            markovchains <- vector()
            for (i in clusters[[1]]){
                mc <- fitMarkovChain(i, order = order) 
                markovchains <- c(markovchains, mc)
            }
            return(markovchains)
        }
        force(centres)
        worker <- function(cls) {
            clusterChainPair <- list()
            clusters <- clusterClickstreams(clickstreamList = cls, centers = centres)
            clusterChainPair <- c(clusterChainPair, list(clusters))
            mc <- fitMarkovChains(clusters)
            clusterChainPair <- c(clusterChainPair, list(mc))
            return (clusterChainPair)
        }
        return(worker)
    }
    
    listOfClickstreams <- list()
    for (i in maxIterations){
        listOfClickstreams <- c(listOfClickstreams, list(trainingCLS))
    }
    parallelCluster <- parallel::makeCluster(cores)
    clusterEvalQ(parallelCluster, {
        library(plyr) 
        library(methods) 
        library(stats)
        library(linprog)
    })
    clusterExport(parallelCluster, c("fitMarkovChain", "fitMarkovChains", "clusterClickstreams"))
    setOfClusterSets <- list()
    print(centres)
    for (c in centres){
        clusters <- parallel::parLapply(parallelCluster, listOfClickstreams, mkWorker(c))
        setOfClusterSets <- c(setOfClusterSets, list(clusters))
    }
    if(!is.null(parallelCluster)) {
        parallel::stopCluster(parallelCluster)
        parallelCluster <- c()
    }
    return (setOfClusterSets)
}

#' Generates an optimal set of clusters for a clickstream based on certain constraints and with parallel computation
#' 
#' @description This is an experimental function for a consensus clustering algorithm based on targeting a range of average next state probabilities derived when fitting each cluster to a markov chain. This function parallelizes k-means and fitToMarkovChain operations across computer cores, and depends on the parallel package to function.
#' @param trainingCLS Clickstream object with training data (this should be the data used to build the markov chain object).
#' @param testCLS Clickstream object with test data.
#' @param maxIterations Number of times to iterate (repeat) through the k-means clustering.
#' @param optimalProbMean The target average probability of each next page click prediction in a 1st order markov chain.
#' @param range The range above the optimal probability to target. 
#' @param centresMin The minimum cluster centres to evaluate.
#' @param clusterCentresRange the additional cluster centres to evaluate.
#' @param order The order for markov chains that will be used to evaluate each cluster.
#' @param cores Number of cores used for clustering.
#' @param takeHighest determines whether to default to the highest mean next click probability, or error if the target is not reached after the given number of k-means iterations. 
#' @param verbose Should this function report extra information on progress?
#' @author Theo van Kraay \email{theo.vankraay@@hotmail.com}
#' @examples
#' training <- c("User1,h,c,c,p,c,h,c,p,p,c,p,p,o",
#'               "User2,i,c,i,c,c,c,d",
#'               "User3,h,i,c,i,c,p,c,c,p,c,c,i,d",
#'               "User4,h,c,c,p,p,c,p,p,p,i,p,o",
#'               "User5,i,h,c,c,p,p,c,p,c,d",
#'               "User6,i,h,c,c,p,p,c,p,c,o",
#'               "User7,i,h,c,c,p,p,c,p,c,d",
#'               "User8,i,h,c,c,p,p,c,p,c,d,o")
#' 
#' test <- c(
#'     "User1,h,c,c,p,c,h,c,p,p,c,p,p,o",
#'     "User2,i,c,i,c,c,c,d",
#'     "User3,h,i,c,i,c,p,c,c,p,c,c,i,d"
#' )
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
#' clusters <- getConsensusClustersParallel(trainingCLS, testCLS, maxIterations=5, 
#'                                  optimalProbMean=0.40, range = 0.70, centresMin = 2, 
#'                                  clusterCentresRange = 0, order = 1, cores = 2,
#'                                  takeHighest = FALSE, verbose = FALSE)
#' markovchains <- fitMarkovChains(clusters)
#' startPattern <- new("Pattern", sequence = c("i", "h", "c", "p"))
#' mc <- getOptimalMarkovChain(startPattern, markovchains, clusters)
#' predict(mc, startPattern)
#' @export getConsensusClustersParallel
getConsensusClustersParallel = function(trainingCLS, testCLS, maxIterations=10, optimalProbMean=0.50, range=0.30, 
                                        centresMin=2, clusterCentresRange=0, order=1, cores=2, takeHighest=FALSE, verbose=FALSE){
    cls <- trainingCLS
    vec <- unlist(cls)
    dedupe <- vec[which(!duplicated(vec))]
    centresMax <- centresMin + clusterCentresRange
    listOfClusters <- list()
    clusterCentres <- centresMin:centresMax
    iterations <- 1:maxIterations
    vectorOfAllProbsMeans <- vector()
    limit <- optimalProbMean + range
    
    if (verbose)
        cat("Getting cluster sets in parallel....\n")
    clusterSets <- .getParallelClusterSets(trainingCLS, iterations, centres=clusterCentres, cores = cores)
    
    if (verbose)
        cat("Starting next page probability aggregation....\n")
    for (i in clusterSets){
        for (c in i){
            clusters <- c[[1]]
            markovchains <- c[[2]]
            absorbingStates <- unique(unlist(lapply(markovchains, FUN = function(x) return(x@absorbingStates))))
            vectorOfProbs <- vector()
            for (d in dedupe){
                if(!d %in% absorbingStates){
                    value <- d[[1]]
                    startPattern <- new("Pattern", sequence = c(value)) 
                    mc <- getOptimalMarkovChain(startPattern, markovchains, clusters)
                    prob <- predict(mc, startPattern)
                    vectorOfProbs <- c(vectorOfProbs, prob@probability)
                }
            }
            vectorOfAllProbsMeans <- c(vectorOfAllProbsMeans, mean(vectorOfProbs))
            listOfClusters <- c(listOfClusters, list(clusters))
        }
        print(vectorOfAllProbsMeans)
        candidates <- which(vectorOfAllProbsMeans>optimalProbMean & vectorOfAllProbsMeans < limit)
        if (verbose)
            cat("Candidates are: ", candidates, "\n")
        # Shutdown cluster neatly
    }
    if (verbose)
        cat("Finished next page probability aggregation....\n")
    if (takeHighest != TRUE){
        if (length(candidates) > 0){
            #get the candidate clusters into a vector
            candidateClusters <- list()
            for (i in candidates) {
                clusters <- listOfClusters[[i]]
                candidateClusters <- c(candidateClusters, list(clusters))
            }
            if (verbose)
                cat("Evaluating candidates....\n")
            vec_variances <- vector()
            for(c in candidateClusters){
                markovchains <- fitMarkovChains(c) 
                variance <- mcEvaluateAllClusters(markovchains, c, testCLS, trainingCLS, returnChiSquareOnly = TRUE)
                if (verbose)
                    cat("Variance is....", variance,"\n")
                vec_variances <- c(vec_variances,variance)
            }
            if (verbose)
                cat("Vector of variances is:", vec_variances,"\n")
            #winner <- which(vec_variances==min(vec_variances))
            winner <- which.min(vec_variances)
            if (verbose)
                cat("Winner is:", winner, "\n")
            return(candidateClusters[[winner]])
        }
        else{
            stop(("Target range was not reached with the given number of iterations."))
        }
    }
    else{
        if (length(candidates) == 0) {
            warning("Target prediction accuracy was not reached with the given number of iterations. Taking highest probability mean.")
        }
        candidates <- which(vectorOfAllProbsMeans ==  max(vectorOfAllProbsMeans))
        return(listOfClusters[[candidates]])
    }
}
