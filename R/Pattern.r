setClass("Pattern", 
         representation(sequence="character", 
                        probability="numeric",
                        absorbingProbabilities="data.frame")
)

setMethod("initialize",
          signature(.Object = "Pattern"),
          function (.Object, sequence, probability, absorbingProbabilities, ...) {
              if(missing(probability)) probability=1 
              if(missing(absorbingProbabilities)) 
                  absorbingProbabilities=data.frame(None=0)
              callNextMethod(.Object, sequence=sequence, probability=probability, 
                             absorbingProbabilities=absorbingProbabilities,...)
          }
) 

setMethod("show", "Pattern",
          function(object) {
              cat("Sequence: ")
              cat(object@sequence)
              cat("\n")
              cat("Probability: ")
              cat(object@probability)
              cat("\n")
              cat("Absorbing Probabilities: \n")
              print(object@absorbingProbabilities)
              cat("\n\n")
          }
)

setMethod("[", signature(x="Pattern", i="ANY"),
          function(x) {
              out=x@sequence[i]
              return(out)
          }
)   

setMethod("+", c("Pattern", "Pattern"),
          function(e1, e2) {
              sequence=c(e1@sequence, e2@sequence)
              probability=e1@probability*e2@probability              
              pattern=new("Pattern", sequence=sequence, probability=probability)
              return(pattern)
          }
)  