#' getMutexGroup function
#' 
#' Given a binary matrix and its corresponding probability matrix pij, compute the Poisson Binomial
#' method to estimate mutual exclusive events.
#'
#' @param A The binary matrix
#' @param PM The corresponding probability matrix of A. It can be computed using function getPM. By default equal to getPM(A)
#' @param type one of Coverage, Exclusivity or Impurity. By default is Impurity
#' @param lower.tail True if mutually exclusive test. False for co-ocurrence. By default is TRUE.
#'
#' @return A symmetric matrix with the p.value of the corresponding test.
#'
#' @examples 
#'
#'   #This first example is a basic 
#'   #example of how to perform getMutexGroup
#'   
#'   data("A_example")
#'   A2 <- A_example[,1:30]
#'   A2[1,1:10] <- 1
#'   A2[2,1:10] <- 0
#'   A2[3,1:10] <- 0
#'   A2[1,11:20] <- 0
#'   A2[2,11:20] <- 1
#'   A2[3,11:20] <- 0
#'   A2[1,21:30] <- 0
#'   A2[2,21:30] <- 0
#'   A2[3,21:30] <- 1
#'   PM2 <- getPM(A2)
#'   A <- A2[1:3,]
#'   PM <- PM2[1:3,]
#'   
#'   getMutexGroup(A, PM, "Impurity")
#'   getMutexGroup(A, PM, "Coverage")
#'   getMutexGroup(A, PM, "Exclusivity")
#'   
#'
#' @import Matrix
#' @import parallel
#' @importFrom speedglm control
#' @importFrom PoissonBinomial ppbinom dpbinom
#' @export

getMutexGroup <- function(A = NULL, PM = NULL, type ="Impurity", lower.tail = TRUE){
  
  if(is.null(A)){
    stop("not input matrix A")
  }
  
  if(!is(A,"matrix") & !is(A,"Matrix")){
    stop("input A must be a Matrix or a matrix class")
  }
  
  if(is.null(PM)){
    stop("not probability matrix PM")
  }
  
  if(!is(PM,"matrix") & !is(PM,"Matrix")){
    stop("input PM must be a Matrix or a matrix class")
  }
  
  if(nrow(A)==0 | ncol(A) == 0){
    stop("input A must have at least 1 row and 1 column")
  }
  
  if(nrow(PM)==0 | ncol(PM) == 0){
    stop("input PM must have at least 1 row and 1 column")
  }
  
  if(max(A)>1){
    stop("input A must be binary")
  }
  
  if(max(PM)>1){
    stop("input PM must be binary")
  }

  if (type == "Impurity") {
    success <- sum((colSums(A) >= 2))
    pp <- 1- apply(PM, 2, function(PM) {ppbinom(1,PM)})
    p.value <- ppbinom(success, pp, lower.tail = T)
  }
  
  if (type == "Coverage") {
    success <- sum((colSums(A) >= 1))
    pp <- 1- apply(PM, 2, function(PM) {ppbinom(0,PM)})
    p.value <- ppbinom(success, pp, lower.tail = F)
  }
  
  if (type == "Exclusivity") {
    success <- sum((colSums(A) == 1))
    pp <- apply(PM, 2, function(PM) {dpbinom(1,PM)})
    p.value <- ppbinom(success, pp, lower.tail = F)
  }
  
  return(p.value)
}


# # Test script for the function. To be removed in the future.
# 
# data("A_example")
# A2 <- A_example[,1:100]
# A2[1,1:20] <- 1
# A2[2,1:20] <- 0
# A2[2,21:40] <- 1
# A2[1,21:40] <- 0
# PM2 <- getPM(A2)
# A <- A2[1:2,]
# PM <- PM2[1:2,]
# image(A) # These two genes are mutually exclusive (to a certain extent)
# 
# 
# getMutexGroup(A, PM, "Impurity")
# getMutex(A, PM)[1,2]
# 
# getMutexGroup(A, PM, "Coverage")
# getMutexGroup(A, PM, "Exclusivity")
