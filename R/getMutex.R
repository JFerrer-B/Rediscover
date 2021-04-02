#' getMutex function
#' 
#' Given a binary matrix and its corresponding probability matrix pij, compute the Poisson Binomial
#' method to estimate mutual exclusive events.
#'
#' @param A The binary matrix
#' @param PM The corresponding probability matrix of A. It can be computed using function getPM. By default equal to getPM(A)
#' @param lower.tail True if mutually exclusive test. False for co-ocurrence. By default is TRUE.
#' @param mixed option to compute lower p.values with an exact method. By default TRUE
#' @param th upper threshold of p.value to apply the exact method.
#' @param verbose The verbosity of the output
#' @param parallel If the exact method is executed with a parallel process.
#'
#' @return A symmetric matrix with the p.value of the corresponding test.
#'
#' @examples 
#' 
#'   #This first example is a basic 
#'   #example of how to perform getMutex. 
#'   
#'   data("A_example")
#'   PMA <- getPM(A_example)
#'   mismutex <- getMutex(A=A_example,PM=PMA)
#'   
#'   \donttest{
#'   #The next example, is the same as the first one but,
#'   # using a matrix of class Matrix. 
#'   
#'   data("A_Matrix")
#'   PMA_Matrix <- getPM(A_Matrix)
#'   mismutex <- getMutex(A=A_Matrix,PM=PMA_Matrix)
#'   
#'   #Finally, the last example, shows a real 
#'   #example of how to perform this function when using 
#'   #data from TCGA, Colon Adenocarcinoma in this case. 
#'   
#'   data("TCGA_COAD")
#'   data("PM_COAD")
#'   
#'   PM_COAD <- getMutex(TCGA_COAD, PM_COAD)
#'   }
#'
#' @import Matrix
#' @import parallel
#' @importFrom stats pnorm dnorm
#' @importFrom speedglm control
#' @importFrom PoissonBinomial ppbinom
#' @export

getMutex <- function(A = NULL, PM = getPM(A), lower.tail = TRUE, 
                     mixed = FALSE,
                     th = 1e-2, verbose = FALSE, parallel = FALSE){
  
  if(verbose){
    message("checking inputs...")
  }
  if(is.null(A)){
    stop("not input matrix A")
  }
  
  if(!is(A,"matrix") & !is(A,"Matrix")){
    stop("input A must be a Matrix or a matrix class")
  }
  
  if(nrow(A)==0 | ncol(A) == 0){
    stop("input A must have at least 1 row and 1 column")
  }
  
  if(max(A)>1){
    stop("input A must be binary")
  }
  
  if(verbose){
    message("Building model...")
  }
  Mevents <- Matrix(A)
  Mevents <- tcrossprod(Mevents) 
  PM <- as.matrix(PM)
  MeanEst <- tcrossprod(PM) # expected means
  varEst <- MeanEst - tcrossprod(PM*PM) # expected variance
  gammEst <- varEst - 2*(MeanEst-varEst) + 2*tcrossprod(PM * PM * PM) # 3rd order correction
  sqrtVarEst <- sqrt(varEst) # expected standard deviations
  
  kk1 <- (Mevents + 0.5 - MeanEst)/sqrtVarEst
  kk1 <- as(kk1, "dspMatrix")
  ind = gammEst/(6 * sqrtVarEst^3)
  ind <- as(ind, "dspMatrix")
  
  pvals <- kk1
  if (lower.tail) {
    ppp <- pnorm(kk1@x, lower.tail = TRUE) + ind@x * (1-(kk1@x)^2)*dnorm(kk1@x)
  } else {
    ppp <- pnorm(kk1@x, lower.tail = FALSE) - ind@x * (1-(kk1@x)^2)*dnorm(kk1@x)
  }
  ppp[ppp < 0] <- 0
  ppp[ppp > 1] <- 1
  
  pvals@x <- ppp
  pvals <- as.matrix(pvals)
  if (mixed) {
    if(verbose){
      message("Performing exact method...")
    }
      # Mixed method: use approximation and, if the p.value is small, compute the exact p.value
      II <- which(pvals < th, arr.ind = TRUE)
      II <- II[II[,2] > II[,1],,drop=F] # Remove half of them
      if (nrow(II) > 0){
        if(parallel) {
          no_cores <- detectCores(logical = FALSE) - 1
          cl <- makeCluster(no_cores)
          if(verbose){
            message("Creating cluster...")
          }
          clusterExport(cl, c("ppbinom"))
          pair <- 1:nrow(II)
          pvalue <- parSapply(cl, pair, function (pair, II, PM, A) {
            genei <- II[pair,1]
            genej <- II[pair,2]
            pp <- PM[genei,] * PM[genej,]
            pvalue <- ppbinom(sum(A[genei,]*A[genej,]), pp, method = "DivideFFT", lower.tail = lower.tail)
          }, II, PM, A)
          pvals[II] <- pvalue
          pvals[II[,c(2,1)]] <- pvalue # To keep symmetry
          
          stopCluster(cl)
        } else {
          pair <- 1:nrow(II)
          pvalue <- sapply(pair, function (pair, II, PM, A) {
            genei <- II[pair,1]
            genej <- II[pair,2]
            pp <- PM[genei,] * PM[genej,]
            pvalue <- ppbinom(sum(A[genei,]*A[genej,]), pp, method = "DivideFFT", lower.tail = lower.tail)
          }, II, PM, A)
          pvals[II] <- pvalue
          pvals[II[,c(2,1),drop=F]] <- pvalue # To keep symmetry
        }
      }
  }
  if(verbose){
    message("Building output...")
  }
  pvals <- as(pvals, "dspMatrix")
  return(pvals)
}


