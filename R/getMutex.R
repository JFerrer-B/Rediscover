#' getMutex function
#' 
#' Given a binary matrix and its corresponding probability matrix pij, compute the Poisson Binomial
#' method to estimate mutual exclusive events.
#'
#' @param A The binary matrix
#' @param PM The corresponding probability matrix of A. It can be computed using function getPM. By default equal to getPM(A)
#' @param lower.tail True if mutually exclusive test. False for co-ocurrence. By default is TRUE.
#' @param method one of the following: "ShiftedBinomial" (default),"Exact", "Binomial", and "RefinedNormal".
#' @param mixed option to compute lower p-values with an exact method. By default TRUE
#' @param th upper threshold of p.value to apply the exact method.
#' @param verbose The verbosity of the output
#' @param parallel If the exact method is executed with a parallel process.
#' @param no_cores number of cores. If not stated number of cores of the CPU - 1
#' 
#' @details we  implemented three different approximations of the Poison-Binomial distribution function:
#' \itemize{
#'  \item "ShiftedBinomial" (by default) that correspond to a shifted Binomial with three parameters (Peköz, Shwartz, Christiansen, & Berlowitz, 2010).
#'  \item"Exact" that use the exact formula using the `PoissonBinomial` Rpackage based on the work from (Biscarri, Zhao, & Brunner, 2018).
#'  \item"Binomial" with two parameters (Cam, 1960).
#'  \item"RefinedNormal" that is based on the work from  (Volkova, 1996).
#' }
#'  If `mixed` option is selected (by default is FALSE), the "Exact" method is computed for p-values lower than a threshold
#'   (`th` parameter, that by default is 0.05). When the exact method is computed, it is possible to parallelize the process by
#'   selecting the option `parallel` (by default FALSE) and setting the number of cores (`no_cores` parameter)
#'
#' @return A symmetric matrix with the p-values of the corresponding test.
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
#' @importFrom stats pnorm dnorm pbeta
#' @importFrom speedglm control
#' @importFrom PoissonBinomial ppbinom
#' @export

getMutex <- function(A = NULL, PM = getPM(A), lower.tail = TRUE, 
                            method = "ShiftedBinomial", mixed = TRUE,
                            th = 5e-2, verbose = FALSE, parallel = FALSE,no_cores=NULL){
  

  
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
  
  if(!method %in% c("Exact", "RefinedNormal", "Binomial", "ShiftedBinomial")){
    stop('method must be "Exact", "RefinedNormal", "Binomial", "ShiftedBinomial"')
  }
  
  if(verbose){
    message("Building model...")
  }
  Mevents <- Matrix(A)
  Mevents <- as.matrix(tcrossprod(Mevents) )
  PM <- as.matrix(PM)
  
  if (method == "ShiftedBinomial") {
    l1 <- tcrossprod(PM) # expected means
    l2 <- tcrossprod(PM*PM)
    l3 <- tcrossprod(PM * PM * PM)
    
    pstar <- (l2-l3)/(l1-l2)
    nstar <- (l1 -l2)/(pstar*(1-pstar))
    sstar <- l1 - nstar*pstar
    rm(l1);rm(l2);rm(l3)
    
    # TODO: since the matrix is symmetric and takes a lot of time, use lower.tri
    # to compute only half of the values. Check the RefinedNormal code.
    II <- lower.tri(pstar, diag = T)
    ppp <- pbeta(pstar[II], pmax(Mevents[II]+1-sstar[II],0),pmax(nstar[II]-Mevents[II]+sstar[II],0), lower.tail = !lower.tail)
    
    # Previous line equivalent to the following (round needed for pbinom) -symmetric trick not used
    # n <- round(nstar)
    # s <- round(sstar)
    # p <- ((nstar * pstar) + (sstar-s))/n
    # ppp <- pbinom(Mevents - s, n,p, lower.tail = lower.tail)
    pvals <- pstar
    pvals[II] <- ppp
    pvals <- as(forceSymmetric(pvals, "L"),"dspMatrix")
  }
  if (method == "Binomial") {
    l1 <- tcrossprod(PM)
    l2 <- tcrossprod(PM*PM)
    
    p <- l2/l1
    n <- l1 / p
    
    II <- lower.tri(p, diag = T)
    
    ppp <- pbeta(p[II], Mevents[II]+1,pmax(n[II]-Mevents[II],0), lower.tail = !lower.tail)
    
    # Previous line equivalent to the following (round needed for pbinom)
    # n <- round(n)
    # p <- A/n
    # ppp <- pbinom(Mevents, n,p, lower.tail = lower.tail)
    
    pvals <- p
    pvals[II] <- ppp
    pvals <- as(forceSymmetric(pvals, "L"),"dspMatrix")
  }
  if (method == "RefinedNormal") {
    MeanEst <- tcrossprod(PM) # expected means
    varEst <- MeanEst - tcrossprod(PM*PM) # expected variance
    gammEst <- varEst - 2*(MeanEst-varEst) + 2*tcrossprod(PM * PM * PM) # 3rd order correction
    sqrtVarEst <- sqrt(varEst) # expected standard deviations
    
    kk1 <- (Mevents + 0.5 - MeanEst)/sqrtVarEst
    kk1 <- as(kk1, "dspMatrix")
    ind <- gammEst/(6 * sqrtVarEst^3)
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
  }
  if (method == "Exact"){
    pvals <- matrix(0, nrow = nrow(A), ncol=nrow(A))
    mixed <- T
    th <- 1.01 # Just to be sure :)
    
  }
  if (mixed) {
    if(verbose){
      message("Improving low p-values with exact method...")
    }
    # Mixed method: use approximation and, if the p.value is small, compute the exact p.value
    II <- which(pvals < th, arr.ind = TRUE)
    II <- II[II[,2] > II[,1],,drop=F] # Remove half of them
    if (nrow(II) > 0){
      if(parallel) {
        if(is.null(no_cores)){
          no_cores <- detectCores(logical = FALSE) - 1          
        }
        cl <- makeCluster(no_cores)
        if(verbose){
          message("Creating cluster...")
        }
        clusterExport(cl, c("ppbinom"))
        pair <- 1:nrow(II)
        pvalue <- parSapply(cl, pair, function (pair, II, PM,Mevents) {
          genei <- II[pair,1]
          genej <- II[pair,2]
          pp <- PM[genei,] * PM[genej,]
          pvalue <- ppbinom(Mevents[genei,genej], pp, method = "DivideFFT", lower.tail = lower.tail)
        }, II, PM,Mevents)
        pvals[II] <- pvalue
        pvals[II[,c(2,1),drop=F]] <- pvalue # To keep symmetry
        
        stopCluster(cl)
      } else {
        pair <- 1:nrow(II)
        pvalue <- sapply(pair, function (pair, II, PM, Mevents) {
          genei <- II[pair,1]
          genej <- II[pair,2]
          pp <- PM[genei,] * PM[genej,]
          pvalue <- ppbinom(Mevents[genei,genej], pp, method = "DivideFFT", lower.tail = lower.tail)
        }, II, PM,Mevents)
        pvals[II] <- pvalue
        pvals[II[,c(2,1),drop=F]] <- pvalue # To keep symmetry
      }
    }
    pvals <- as(pvals, "dspMatrix")
  }
  if(verbose){
    message("Building output...")
  }
  pvals <- as(pvals, "dspMatrix")
  return(pvals)
}


