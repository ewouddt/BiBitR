#' @title A biclustering algorithm for extracting bit-patterns from binary datasets
#'
#' @description BiBitR is a simple R wrapper which directly calls the original Java code for applying the BiBit algorithm.
#' The original Java code can be found at \url{http://eps.upo.es/bigs/BiBit.html} by Domingo S. Rodriguez-Baena, Antonia J. Perez-Pulido and Jesus S. Aguilar-Ruiz.
#' 
#' The BiBitR package also includes the following functions and/or workflows:
#' \itemize{
#' \item A slightly adapted version of the original BiBit algorithm which now allows allows noise when adding rows to the bicluster (\code{\link{bibit2}}).
#' \item A function which accepts a pattern and, using the BiBit algorithm, will find biclusters fully or partly fitting the given pattern (\code{\link{bibit3}}).
#' \item A workflow which can discover larger patterns (and their biclusters) using BiBit and classic hierarchical clustering approaches (\code{\link{BiBitWorkflow}}).
#' }
#' 
#' 
#' @references Domingo S. Rodriguez-Baena, Antonia J. Perez-Pulido and Jesus S. Aguilar-Ruiz (2011), "A biclustering algorithm for extracting bit-patterns from binary datasets", \emph{Bioinformatics}
#'
#' @docType package
#' @name BiBitR
NULL