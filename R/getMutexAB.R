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
#' @param method one of the following: "ShiftedBinomial" (default),"Exact", "RefinedNormal", and "Binomial".
#' @param mixed option to compute lower p-values with an exact method. By default TRUE
#' @param th upper threshold of p-value to apply the exact method.
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
#' @return A  matrix with the p-values of the corresponding test.
#'
#' @examples 
#' 
#'   \donttest{   
#'   
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
#' @import parallel
#' @importFrom utils combn
#' @importFrom stats pnorm dnorm pbeta
#' @importFrom speedglm control
#' @importFrom PoissonBinomial ppbinom
#' @importFrom matrixStats rowProds
#' @importFrom stats rnorm
#' @importFrom ShiftConvolvePoibin ppoisbin
#' @export

getMutexAB <- function(A, PMA = getPM(A), B, PMB = getPM(B), lower.tail = TRUE, 
                       method = "ShiftedBinomial",mixed = TRUE,
                     th = 5e-2, verbose = FALSE, parallel = FALSE,no_cores=NULL){
  
  normal <- FALSE
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
  
  if(!method %in% c("Exact", "RefinedNormal", "Binomial", "ShiftedBinomial")){
    stop('method must be "Exact", "RefinedNormal", "Binomial", "ShiftedBinomial"')
  }
  
  if(verbose){
    message("Building model...")
  }
  
  Mevents <- A %*% t(B)
  Mevents <- Matrix(Mevents)
  
  
  if(ncol(A)<3800){
    ppoisbin_method <- "DC"  
  }else{
    ppoisbin_method <- "ShiftConvolve"
  }
  
  
  ###### RefinedNormal ######
  if (method == "RefinedNormal"){
    PMA_2 <- as.matrix(PMA)
    PMB_2 <- as.matrix(PMB)
    MeanEst <- PMA_2 %*% t(PMB_2) # expected means
    varEst <- MeanEst - ( (PMA_2*PMA_2) %*% t(PMB_2*PMB_2) ) # expected variance
    gammEst <- varEst - 2*(MeanEst-varEst) + 2*((PMA_2 * PMA_2 * PMA_2) %*% t(PMB_2 * PMB_2 * PMB_2)) # 3rd order correction
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
    
  }
  ###### ShiftedBinomial ######
  if (method == "ShiftedBinomial"){
    PMA_2 <- as.matrix(PMA)
    PMB_2 <- as.matrix(PMB)
    l1 <- PMA_2 %*% t(PMB_2) # expected means
    l2 <- (PMA_2*PMA_2) %*% t(PMB_2*PMB_2)
    l3 <- (PMA_2 * PMA_2 * PMA_2) %*% t(PMB_2 * PMB_2 * PMB_2)
    pstar <- (l2-l3)/(l1-l2)
    nstar <- (l1 -l2)/(pstar*(1-pstar))
    sstar <- l1 - nstar*pstar
    rm(l1);rm(l2);rm(l3)
    
    pvals <- Mevents
    # TODO: since the matrix is symmetric and takes a lot of time, use lower.tri
    # to compute only half of the values. Check the RefinedNormal code.
    # II <- lower.tri(pstar, diag = T)
    # ppp <- pbeta(pstar[II], pmax(Mevents[II]+1-sstar[II],0),pmax(nstar[II]-Mevents[II]+sstar[II],0), lower.tail = !lower.tail)
    ppp <- pbeta(as.vector(pstar),
                 as.vector(pmax(Mevents+1-sstar,0)),
                 as.vector(pmax(nstar-Mevents+sstar,0)), lower.tail = !lower.tail)
    # Previous line equivalent to the following (round needed for pbinom) -symmetric trick not used
    # n <- round(nstar)
    # s <- round(sstar)
    # p <- ((nstar * pstar) + (sstar-s))/n
    # ppp <- pbinom(Mevents - s, n,p, lower.tail = lower.tail)
    # pvals <- pstar
    pvals@x <- ppp
    
    
    
    
    # pvals <- as(forceSymmetric(pvals, "L"),"dspMatrix")
  }
  ###### Binomial ######
  if (method == "Binomial"){
    PMA_2 <- as.matrix(PMA)
    PMB_2 <- as.matrix(PMB)
    l1 <- PMA_2 %*% t(PMB_2)
    l2 <- (PMA_2*PMA_2) %*% t(PMB_2*PMB_2)
    
    p <- l2/l1
    n <- l1 / p
    
    # II <- lower.tri(p, diag = T)
    
    ppp <- pbeta(as.vector(p),
                 as.vector(Mevents+1),
                 as.vector(pmax(n-Mevents,0)),
                 lower.tail = !lower.tail)
    
    # Previous line equivalent to the following (round needed for pbinom)
    # n <- round(n)
    # p <- A/n
    # ppp <- pbinom(Mevents, n,p, lower.tail = lower.tail)
    
    pvals <- Mevents
    pvals@x <- ppp
    
    
    
    
    # pvals <- as(forceSymmetric(pvals, "L"),"dspMatrix")
  }
  ###### exact ######
  if (method == "Exact"){
    pvals <- matrix(0, nrow = nrow(A), ncol=nrow(B))
    mixed <- FALSE
    # th <- 1.01 # Just to be sure :)
    genes_factor_A <- factor(PMA@rowExps)
    Idx_A <- as.matrix(sparseMatrix(i=as.numeric(genes_factor_A),j = 1:nrow(PMA),x = 1))
    miniPM_A <- as.matrix(PMA[match(levels(genes_factor_A),PMA@rowExps),])
    
    genes_factor_B <- factor(PMB@rowExps)
    Idx_B <- as.matrix(sparseMatrix(i=as.numeric(genes_factor_B),j = 1:nrow(PMB),x = 1))
    miniPM_B <- as.matrix(PMB[match(levels(genes_factor_B),PMB@rowExps),])
    
    llx <- expand.grid_fast(1:nrow(miniPM_A),1:nrow(miniPM_B))
    # miniPM_2 <- miniPM_A[llx[,1],]*miniPM_B[llx[,2],]
    
    if(parallel){
      if(is.null(no_cores)){
        no_cores <- detectCores(logical = FALSE) - 1          
      }
      my_os <- get_os()
      if(my_os=="linux" | my_os=="osx"){
        type_cl <- "FORK"
      }else{
        type_cl <- "PSOCK"
      }
      cl <- makeCluster(no_cores,type = type_cl)
      if(verbose){
        message("Creating cluster...")
      }
      # clusterExport(cl, c("Matrix","ppbinom","ppoisbin","expand.grid_fast"))
      clusterExport(cl, c("Matrix","ppbinom","ppoisbin"))
      i <- 1:nrow(llx)
      pvalue <- parLapply(cl, i, function (i, Idx_A,Idx_B,llx, miniPM_A,miniPM_B,Mevents) {
        # idx_kk <- expand.grid_fast(which(Idx_A[llx[i,1],]==1),which(Idx_B[llx[i,2],]==1))
        a <- which(Idx_A[llx[i,1],]==1)
        b <- which(Idx_B[llx[i,2],]==1)
        idx_kk <- cbind(rep(a,each=length(b)),b)
        
        # mi_pp <- miniPM_2[i,]
        mi_pp <- miniPM_A[llx[i,1],]*miniPM_B[llx[i,2],]
        
        pvals <- ppoisbin(Mevents[idx_kk], mi_pp, method = "DC", lower.tail = TRUE)
        if(any(Mevents[idx_kk]==0)){
          oox <- which(Mevents[idx_kk]==0)
          pvalue <- prod(1-mi_pp)
          if(lower.tail==FALSE){
            pvalue <- 1-pvalue
          }
          pvals[oox] <- pvalue
        }
        if(any(Mevents[idx_kk]==1)){
          oox <- which(Mevents[idx_kk]==1)
          pvalue <- prod(1-mi_pp) * (1+sum( mi_pp/(1-mi_pp)))
          if(lower.tail==FALSE){
            pvalue <- 1-pvalue
          }
          pvals[oox] <- pvalue
        }
        return(cbind(idx_kk,pvals))
        
      }, Idx_A,Idx_B,llx, miniPM_A,miniPM_B,Mevents)
      stopCluster(cl)
      pvalue <- do.call(rbind,pvalue)
      pvals[cbind(pvalue[,1],pvalue[,2])] <- pvalue[,3]
      pvals[cbind(pvalue[,2],pvalue[,1])] <- pvalue[,3]
    }else{
      for(i in 1:nrow(llx)){
        # i <- 1
        # idx_kk <- as.matrix(expand.grid(which(Idx[llx[1,i],]==1),which(Idx[llx[2,i],]==1)))
        idx_kk <- expand.grid_fast(which(Idx_A[llx[i,1],]==1),which(Idx_B[llx[i,2],]==1))
        
        # mi_pp <- miniPM_2[i,]
        mi_pp <- miniPM_A[llx[i,1],]*miniPM_B[llx[i,2],]
        
        pvals[idx_kk] <- ppoisbin(Mevents[idx_kk], mi_pp, method = ppoisbin_method, lower.tail = lower.tail)
        if(any(Mevents[idx_kk]==0)){
          oox <- which(Mevents[idx_kk]==0)
          pvalue <- prod(1-mi_pp)
          if(lower.tail==FALSE){
            pvalue <- 1-pvalue
          }
          pvals[idx_kk][oox] <- pvalue
        }
        if(any(Mevents[idx_kk]==1)){
          oox <- which(Mevents[idx_kk]==1)
          pvalue <- prod(1-mi_pp) * (1+sum( mi_pp/(1-mi_pp)))
          if(lower.tail==FALSE){
            pvalue <- 1-pvalue
          }
          pvals[idx_kk][oox] <- pvalue
        }
        pvals[idx_kk[,c(2,1),drop=F]] <- pvals[idx_kk] # To keep symmetry
      }
    }
    
    # diag(pvals) <- 0
  }
  ####### mixed #######
  if (mixed) {
    if(verbose){
      message("Performing exact method...")
    }
    # Mixed method: use approximation and, if the p-value is small, compute the exact p-value
    # browser()
    II <- which(pvals < th, arr.ind = TRUE)
    if (nrow(II) > 0){
      if(parallel) {
        if(is.null(no_cores)){
          no_cores <- detectCores(logical = FALSE) - 1          
        }
        
        my_os <- get_os()
        if(my_os=="linux" | my_os=="osx"){
          type_cl <- "FORK"
        }else{
          type_cl <- "PSOCK"
        }
        cl <- makeCluster(no_cores,type = type_cl)
        
        if(verbose){
          message("Creating cluster...")
        }
        clusterExport(cl, c("Matrix","ppbinom","ppoisbin"))
        
        genes_factor_A <- factor(PMA@rowExps)
        Idx_A <- as.matrix(sparseMatrix(i=as.numeric(genes_factor_A),j = 1:nrow(PMA),x = 1))
        miniPM_A <- as.matrix(PMA[match(levels(genes_factor_A),PMA@rowExps),])
        
        genes_factor_B <- factor(PMB@rowExps)
        Idx_B <- as.matrix(sparseMatrix(i=as.numeric(genes_factor_B),j = 1:nrow(PMB),x = 1))
        miniPM_B <- as.matrix(PMB[match(levels(genes_factor_B),PMB@rowExps),])
        
        llx <- expand.grid_fast(1:nrow(miniPM_A),1:nrow(miniPM_B))
        
        II_2 <- cbind(which(Idx_A[,II[,1]] == 1,arr.ind = T)[,1],
                      which(Idx_B[,II[,2]] == 1,arr.ind = T)[,1])
        
        pos <- factor(II_2 %*% rnorm(2))
        pvals_2 <- as.matrix(pvals)
        
        i <- 1:length(levels(pos))
        p_to_ad <- parLapply(cl, i, function (i, II_2,pos,miniPM_A,miniPM_B,II,pvals_2,th,Mevents) {
          ix <- II_2[match(levels(pos)[i],pos),]
          mi_pp <- miniPM_A[ix[1],]*miniPM_B[ix[2],]
          idx_kk <- II[which(pos == levels(pos)[i]),,drop=F]
          
          # idx_kk <- idx_kk[which(pvals[idx_kk]<th),,drop=F]
          
          idx_kk <- idx_kk[which(pvals_2[idx_kk]<th),,drop=F]
          
          miskk <- Mevents[idx_kk]
          mispvalues <- vector(mode="numeric",length=length(miskk))
          # mispvalues <- ppoisbin(miskk, mi_pp, method = ppoisbin_method, lower.tail = lower.tail)
          oox_0 <- c()
          oox_1 <- c()
          if(any(Mevents[idx_kk]==0)){
            oox_0 <- which(Mevents[idx_kk]==0)
            pvalue <- prod(1-mi_pp)
            if(lower.tail==FALSE){
              pvalue <- 1-pvalue
            }
            # pvals[idx_kk][oox] <- pvalue
            mispvalues[oox_0] <- pvalue
            if(length(oox_0) == length(mispvalues)){
              return(cbind(idx_kk,mispvalues))
            }
          }
          if(any(Mevents[idx_kk]==1)){
            oox_1 <- which(Mevents[idx_kk]==1)
            pvalue <- prod(1-mi_pp) * (1+sum( mi_pp/(1-mi_pp)))
            if(lower.tail==FALSE){
              pvalue <- 1-pvalue
            }
            # pvals[idx_kk][oox] <- pvalue
            mispvalues[oox_1] <- pvalue
            if(length(oox_1) == length(mispvalues)){
              return(cbind(idx_kk,mispvalues))
            }
          }
          if(length(oox_0) > 0 | length(oox_1) > 0){
            if(length(c(oox_0,oox_1)) == length(mispvalues)){
              return(cbind(idx_kk,mispvalues))
            }else{
              mispvalues[-c(oox_0,oox_1)] <- ppoisbin(miskk[-c(oox_0,oox_1)], mi_pp, method = ppoisbin_method, lower.tail = lower.tail) 
            }
          }else{
            mispvalues <- ppoisbin(miskk, mi_pp, method = ppoisbin_method, lower.tail = lower.tail)
          }
          return(cbind(idx_kk,mispvalues))
        },II_2,pos,miniPM_A,miniPM_B,II,pvals_2,th,Mevents)
        stopCluster(cl)
        p_to_ad <- do.call(rbind,p_to_ad)  
        pvals[p_to_ad[,c(1,2),drop=F]] <- p_to_ad[,3]
        
      }else{
        if(normal){
          
          pair <- 1:nrow(II)
          pvalue <- sapply(pair, function (pair, II, PMA, PMB,Mevents) {
            genei <- II[pair,1]
            genej <- II[pair,2]
            pp <- PMA[genei,] * PMB[genej,]
            if(Mevents[genei,genej]==0){
              pvalue <- prod(1-pp)
              if(lower.tail==FALSE){
                pvalue <- 1-pvalue
              }
            }else if(Mevents[genei,genej]==1){
              pvalue <- prod(1-pp) * (1+sum( pp/(1-pp)))
              if(lower.tail==FALSE){
                pvalue <- 1-pvalue
              }
            }else{
              # pvalue <- ppbinom(sum(A[genei,]*B[genej,]), pp, method = "DivideFFT", lower.tail = lower.tail)  
              pvalue <- ppoisbin(Mevents[genei,genej], pp, method = ppoisbin_method, lower.tail = lower.tail)    
            }
            # return(pvalue)
            
          },II, PMA, PMB,Mevents)
          pvals[II] <- pvalue
          
        }else{
          genes_factor_A <- factor(PMA@rowExps)
          Idx_A <- as.matrix(sparseMatrix(i=as.numeric(genes_factor_A),j = 1:nrow(PMA),x = 1))
          miniPM_A <- as.matrix(PMA[match(levels(genes_factor_A),PMA@rowExps),])
          
          genes_factor_B <- factor(PMB@rowExps)
          Idx_B <- as.matrix(sparseMatrix(i=as.numeric(genes_factor_B),j = 1:nrow(PMB),x = 1))
          miniPM_B <- as.matrix(PMB[match(levels(genes_factor_B),PMB@rowExps),])
          
          llx <- expand.grid_fast(1:nrow(miniPM_A),1:nrow(miniPM_B))
          
          II_2 <- cbind(which(Idx_A[,II[,1]] == 1,arr.ind = T)[,1],
                        which(Idx_B[,II[,2]] == 1,arr.ind = T)[,1])
          
          pos <- factor(II_2 %*% rnorm(2))
          pvals_2 <- as.matrix(pvals)
          
          p_to_ad <- lapply(1:length(levels(pos)),function(i){
            ix <- II_2[match(levels(pos)[i],pos),]
            mi_pp <- miniPM_A[ix[1],]*miniPM_B[ix[2],]
            idx_kk <- II[which(pos == levels(pos)[i]),,drop=F]
            
            # idx_kk <- idx_kk[which(pvals[idx_kk]<th),,drop=F]
            
            idx_kk <- idx_kk[which(pvals_2[idx_kk]<th),,drop=F]
            
            miskk <- Mevents[idx_kk]
            mispvalues <- vector(mode="numeric",length=length(miskk))
            # mispvalues <- ppoisbin(miskk, mi_pp, method = ppoisbin_method, lower.tail = lower.tail)
            oox_0 <- c()
            oox_1 <- c()
            if(any(Mevents[idx_kk]==0)){
              oox_0 <- which(Mevents[idx_kk]==0)
              pvalue <- prod(1-mi_pp)
              if(lower.tail==FALSE){
                pvalue <- 1-pvalue
              }
              # pvals[idx_kk][oox] <- pvalue
              mispvalues[oox_0] <- pvalue
              if(length(oox_0) == length(mispvalues)){
                return(cbind(idx_kk,mispvalues))
              }
            }
            if(any(Mevents[idx_kk]==1)){
              oox_1 <- which(Mevents[idx_kk]==1)
              pvalue <- prod(1-mi_pp) * (1+sum( mi_pp/(1-mi_pp)))
              if(lower.tail==FALSE){
                pvalue <- 1-pvalue
              }
              # pvals[idx_kk][oox] <- pvalue
              mispvalues[oox_1] <- pvalue
              if(length(oox_1) == length(mispvalues)){
                return(cbind(idx_kk,mispvalues))
              }
            }
            if(length(oox_0) > 0 | length(oox_1) > 0){
              if(length(c(oox_0,oox_1)) == length(mispvalues)){
                return(cbind(idx_kk,mispvalues))
              }else{
                mispvalues[-c(oox_0,oox_1)] <- ppoisbin(miskk[-c(oox_0,oox_1)], mi_pp, method = ppoisbin_method, lower.tail = lower.tail) 
              }
            }else{
              mispvalues <- ppoisbin(miskk, mi_pp, method = ppoisbin_method, lower.tail = lower.tail)
            }
            return(cbind(idx_kk,mispvalues))
          })
          
          p_to_ad <- do.call(rbind,p_to_ad) 
          pvals[p_to_ad[,c(1,2),drop=F]] <- p_to_ad[,3]
        }
      }
      pvals <- as(pvals,"dgCMatrix")
    }
  }
  
  # pvals <- as(pvals, "dspMatrix")
  return(pvals)
}


