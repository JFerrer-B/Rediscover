#' getMutexAB function
#' 
#' Given two binary matrices and its corresponding probability matrices PAij and PBij, compute the Poisson Binomial
#' method to estimate mutual exclusive events between A and B
#' 
#' @param A The binary matrix of events A
#' @param PMA The corresponding probability matrix of A. It can be computed using function getPM. By default equal to getPM(A)
#' @param B The binary matrix of events B
#' @param PMB The corresponding probability matrix of B. It can be computed using function getPM. By default equal to getPM(B)
#' @param lower.tail True if mutually exclusive test. False for co-ocurrence. By default is TRUE.
#' @param mixed option to compute lower p.values with an exact method. By default TRUE
#' @param th upper threshold of p.value to apply the exact method.
#' @param verbose The verbosity of the output
#' @param parallel If the exact method is executed with a parallel process.
#'
#' @return A  matrix with the p.value of the corresponding test.
#'
#' @examples 
#'   
#'   #This first example is a basic 
#'   #example of how to perform getMutexAB. 
#'   
#'   data("A_example")
#'   data("B_example")
#'   PMA <- getPM(A_example)
#'   PMB <- getPM(B_example)
#'   mismutex <- getMutexAB(A=A_example, PM=PMA, B=B_example, PMB = PMB)
#'   
#'   \donttest{
#'   #The next example, is the same as the first
#'   # one but, using a matrix of class Matrix. 
#'   
#'   data("A_Matrix")
#'   data("B_Matrix")
#'   PMA <- getPM(A_Matrix)
#'   PMB <- getPM(B_Matrix)
#'   mismutex <- getMutexAB(A=A_Matrix, PM=PMA, B=B_Matrix, PMB = PMB)
#'   
#'   #Finally, the last example, shows a 
#'   #real example of how to perform this function
#'   # when using data from TCGA, Colon Adenocarcinoma in this case. 
#'   
#'   data("TCGA_COAD_AMP")
#'   data("AMP_COAD")
#'   data("PM_TCGA_COAD_AMP")
#'   data("PM_AMP_COAD")
#'   
#'   mismutex <- getMutexAB(A=TCGA_COAD_AMP, 
#'                          PMA=PM_TCGA_COAD_AMP,
#'                          B=AMP_COAD,
#'                          PMB = PM_AMP_COAD)
#'  }
#'
#' @import Matrix
#' @importFrom stats dnorm
#' @importFrom speedglm control
#' @importFrom PoissonBinomial ppbinom
#' @export

getMutexAB <- function(A, PMA = getPM(A), B, PMB = getPM(B), lower.tail = TRUE, 
                     mixed = FALSE,
                     th = 1e-2, verbose = FALSE, parallel = FALSE){
  
  if(verbose){
    message("checking inputs...")
  }
  
  if(is.null(A)){
    stop("not input matrix A")
  }
  
  if(is.null(B)){
    stop("not input matrix B")
  }
  
  if(!is(A,"matrix") & !is(A,"Matrix")){
    stop("input A must be a Matrix or a matrix class")
  }
  
  if(!is(B,"matrix") & !is(B,"Matrix")){
    stop("input B must be a Matrix or a matrix class")
  }
  
  if(nrow(A)==0 | ncol(A) == 0){
    stop("input A must have at least 1 row and 1 column")
  }
  
  if(nrow(B)==0 | ncol(B) == 0){
    stop("input B must have at least 1 row and 1 column")
  }
  
  if(max(A)>1){
    stop("input A must be binary")
  }
  
  if(max(B)>1){
    stop("input B must be binary")
  }
  
  if(verbose){
    message("Building model...")
  }
  
  Mevents <- A %*% t(B) 
  PMA <- as.matrix(PMA)
  PMB <- as.matrix(PMB)
  MeanEst <- PMA %*% t(PMB) # expected means
  varEst <- MeanEst - ( (PMA*PMA) %*% t(PMB*PMB) ) # expected variance
  gammEst <- varEst - 2*(MeanEst-varEst) + 2*((PMA * PMA * PMA) %*% t(PMB * PMB * PMB)) # 3rd order correction
  sqrtVarEst <- sqrt(varEst) # expected standard deviations
  
  kk1 <- (Mevents + 0.5 - MeanEst)/sqrtVarEst
  kk1 <- as(kk1, "Matrix")
  ind <- gammEst/(6 * sqrtVarEst^3)
  ind <- as(ind, "Matrix")
  
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
      if (nrow(II) > 0){
        if (parallel) {
          no_cores <- detectCores(logical = FALSE) - 1
          cl <- makeCluster(no_cores)
          if(verbose){
            message("Creating cluster...")
          }
          clusterExport(cl, c("ppbinom"))
          pair <- 1:nrow(II)
          pvalue <- parSapply(cl, pair, function (pair, II, A, PMA, B, PMB) {
            genei <- II[pair,1]
            genej <- II[pair,2]
            pp <- PMA[genei,] * PMB[genej,]
            pvalue <- ppbinom(sum(A[genei,]*B[genej,]), pp, 
                              method = "DivideFFT", lower.tail = lower.tail)
          },II, A, PMA, B, PMB)
          pvals[II] <- pvalue
          stopCluster(cl)
        } else {
          pair <- 1:nrow(II)
          pvalue <- sapply(pair, function (pair, II, A, PMA, B, PMB) {
            genei <- II[pair,1]
            genej <- II[pair,2]
            pp <- PMA[genei,] * PMB[genej,]
            pvalue <- ppbinom(sum(A[genei,]*B[genej,]), pp, 
                              method = "DivideFFT", lower.tail = lower.tail)
          },II, A, PMA, B, PMB)
          pvals[II] <- pvalue
        }
      }
  }
  # pvals <- as(pvals, "dspMatrix")
  return(pvals)
}


