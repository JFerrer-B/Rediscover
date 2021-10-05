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
#' @importFrom utils combn
#' @importFrom stats pnorm dnorm pbeta
#' @importFrom speedglm control
#' @importFrom PoissonBinomial ppbinom
#' @importFrom ShiftConvolvePoibin ppoisbin
#' @importFrom matrixStats rowProds
#' @importFrom stats rnorm
#' @export

getMutex <- function(A = NULL, PM = getPM(A), lower.tail = TRUE, 
                            method = "ShiftedBinomial", mixed = TRUE,
                            th = 5e-2, verbose = FALSE, parallel = FALSE,no_cores=NULL){
  
  normal <- FALSE #for mutex option, use enhanced algorithm
  
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
  
  
  if(ncol(A)<3800){
    ppoisbin_method <- "DC"  
  }else{
    ppoisbin_method <- "ShiftConvolve"
  }
  
  
  ###### ShiftedBinomial ########
  if (method == "ShiftedBinomial") {
    PM_2 <- as.matrix(PM)
    l1 <- tcrossprod(PM_2) # expected means
    l2 <- tcrossprod(PM_2*PM_2)
    l3 <- tcrossprod(PM_2 * PM_2 * PM_2)
    
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
  ###### Binomial ########
  if (method == "Binomial") {
    PM_2 <- as.matrix(PM)
    l1 <- tcrossprod(PM_2)
    l2 <- tcrossprod(PM_2*PM_2)
    
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
  ####### RefinedNormal #########
  if (method == "RefinedNormal") {
    PM_2 <- as.matrix(PM)
    MeanEst <- tcrossprod(PM_2) # expected means
    varEst <- MeanEst - tcrossprod(PM_2*PM_2) # expected variance
    gammEst <- varEst - 2*(MeanEst-varEst) + 2*tcrossprod(PM_2 * PM_2 * PM_2) # 3rd order correction
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
  ####### exact ######
  if (method == "Exact"){
    
    pvals <- matrix(0, nrow = nrow(A), ncol=nrow(A))
    mixed <- FALSE
    # th <- 1.01 # Just to be sure :)
    genes_factor <- factor(PM@rowExps)
    Idx <- as.matrix(sparseMatrix(i=as.numeric(genes_factor),j = 1:nrow(PM),x = 1))
    miniPM <- as.matrix(PM[match(levels(genes_factor),PM@rowExps),])
    llx <- combn(nrow(miniPM),2,)
    llx <- cbind(llx,rbind(c(1:nrow(miniPM)),c(1:nrow(miniPM))))
    # miniPM_2 <- miniPM[llx[1,],]*miniPM[llx[2,],]
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
      # clusterExport(cl, c("ppbinom","ppoisbin","expand.grid_fast"))
      clusterExport(cl, c("ppbinom","ppoisbin"))
      
      i <- 1:ncol(llx)
      pvalue <- parSapply(cl, i, function (i, Idx,llx, miniPM,Mevents) {
        
        # idx_kk <- expand.grid_fast(which(Idx[llx[1,i],]==1),which(Idx[llx[2,i],]==1))
        a <- which(Idx[llx[1,i],]==1)
        b <- which(Idx[llx[2,i],]==1)
        idx_kk <- cbind(rep(a,each=length(b)),b)
        
        # mi_pp <- miniPM_2[i,]
        mi_pp <- miniPM[llx[1,i],]*miniPM[llx[2,i],]
        
        pvals <- ppoisbin(Mevents[idx_kk], mi_pp, method = ppoisbin_method, lower.tail = lower.tail)
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
          pvalue <- prod(1-mi_pp) * (1+sum(mi_pp/(1-mi_pp)))
          if(lower.tail==FALSE){
            pvalue <- 1-pvalue
          }
          pvals[oox] <- pvalue
        }
        return(cbind(idx_kk,pvals))
      }, Idx,llx, miniPM,Mevents)
      stopCluster(cl)
      pvalue <- do.call(rbind,pvalue)
      pvals[cbind(pvalue[,1],pvalue[,2])] <- pvalue[,3]
      pvals[cbind(pvalue[,2],pvalue[,1])] <- pvalue[,3]
    }else{
      for(i in 1:ncol(llx)){
        # i <- 1
        # idx_kk <- as.matrix(expand.grid(which(Idx[llx[1,i],]==1),which(Idx[llx[2,i],]==1)))
        # idx_kk <- expand.grid_fast(which(Idx[llx[1,i],]==1),which(Idx[llx[2,i],]==1))
        a <- which(Idx[llx[1,i],]==1)
        b <- which(Idx[llx[2,i],]==1)
        idx_kk <- cbind(rep(a,each=length(b)),b)
        # mi_pp <- miniPM_2[i,]
        mi_pp <- miniPM[llx[1,i],]*miniPM[llx[2,i],]
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
    
    diag(pvals) <- 0
  }
  ####### mixed #######
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
        clusterExport(cl, c("ppbinom","ppoisbin"))
        if(normal){
          pair <- 1:nrow(II)
          pvalue <- parSapply(cl, pair, function (pair, II, PM,Mevents) {
            genei <- II[pair,1]
            genej <- II[pair,2]
            pp <- PM_2[genei,] * PM_2[genej,]
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
              # pvalue <- ppbinom(Mevents[genei,genej], pp, method = "DivideFFT", lower.tail = lower.tail)
              pvalue <- ppoisbin(Mevents[genei,genej], pp, method = ppoisbin_method, lower.tail = lower.tail)    
            }
            return(pvalue)
          }, II, PM,Mevents)
          pvals[II] <- pvalue
          pvals[II[,c(2,1),drop=F]] <- pvalue # To keep symmetry
          
          stopCluster(cl)
        }else{
          genes_factor <- factor(PM@rowExps)
          Idx <- as.matrix(sparseMatrix(i=as.numeric(genes_factor),j = 1:nrow(PM),x = 1))
          miniPM <- as.matrix(PM[match(levels(genes_factor),PM@rowExps),])
          llx <- combn(nrow(miniPM),2)
          llx <- cbind(llx,rbind(c(1:nrow(miniPM)),c(1:nrow(miniPM))))
          II_2 <- cbind(which(Idx[,II[,1]] == 1,arr.ind = T)[,1],
                        which(Idx[,II[,2]] == 1,arr.ind = T)[,1])
          pos <- factor(II_2 %*% rnorm(2))
          
          # cat("length_pos=",length(levels(pos)),"\nII_2 = ",nrow(II_2),"\n")
          pvals_2 <- as.matrix(pvals)
          
          # p_to_ad <- sapply(1:length(levels(pos)),function(i){
          i <- 1:length(levels(pos))
          p_to_ad <- parLapply(cl, i, function (i, II_2,pos,miniPM,II,pvals_2,th,Mevents) {
            ix <- II_2[match(levels(pos)[i],pos),]
            mi_pp <- miniPM[ix[1],]*miniPM[ix[2],]
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
          },II_2,pos,miniPM,II,pvals_2,th,Mevents)
          stopCluster(cl)
          if(is.list(p_to_ad)){
            p_to_ad <- do.call(rbind,p_to_ad)  
          }else{
            p_to_ad <- t(p_to_ad)
          }
          pvals[p_to_ad[,c(1,2),drop=F]] <- p_to_ad[,3] # To keep symmetry
          pvals[p_to_ad[,c(2,1),drop=F]] <- p_to_ad[,3] # To keep symmetry
          
          
        }
        
      } else {
        if(normal){
          pair <- 1:nrow(II)
          PM <- as.matrix(PM)
          pvalue <- sapply(pair, function (pair, II, PM, Mevents) {
            genei <- II[pair,1]
            genej <- II[pair,2]
            pp <- PM[genei,] * PM[genej,]
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
              # pvalue <- ppbinom(Mevents[genei,genej], pp, method = "DivideFFT", lower.tail = lower.tail)
              pvalue <- ppoisbin(Mevents[genei,genej], pp, method = ppoisbin_method, lower.tail = lower.tail)    
              
            }
            return(pvalue)
          }, II, PM,Mevents)
          pvals[II] <- pvalue
          pvals[II[,c(2,1),drop=F]] <- pvalue # To keep symmetry
        }else{
          genes_factor <- factor(PM@rowExps)
          Idx <- as.matrix(sparseMatrix(i=as.numeric(genes_factor),j = 1:nrow(PM),x = 1))
          miniPM <- as.matrix(PM[match(levels(genes_factor),PM@rowExps),])
          llx <- combn(nrow(miniPM),2)
          llx <- cbind(llx,rbind(c(1:nrow(miniPM)),c(1:nrow(miniPM))))
          II_2 <- cbind(which(Idx[,II[,1]] == 1,arr.ind = T)[,1],
                        which(Idx[,II[,2]] == 1,arr.ind = T)[,1])
          pos <- factor(II_2 %*% rnorm(2))
          
          # cat("length_pos=",length(levels(pos)),"\nII_2 = ",nrow(II_2),"\n")
          pvals_2 <- as.matrix(pvals)
          # for(i in 1:length(levels(pos))){
          p_to_ad <- lapply(1:length(levels(pos)),function(i){
            ix <- II_2[match(levels(pos)[i],pos),]
            mi_pp <- miniPM[ix[1],]*miniPM[ix[2],]
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
          if(is.list(p_to_ad)){
            p_to_ad <- do.call(rbind,p_to_ad)  
          }else{
            p_to_ad <- t(p_to_ad)
          }
          pvals[p_to_ad[,c(1,2),drop=F]] <- p_to_ad[,3] # To keep symmetry
          pvals[p_to_ad[,c(2,1),drop=F]] <- p_to_ad[,3] # To keep symmetry
          
          # pvals[idx_kk[,c(2,1),drop=F]] <- pvals[idx_kk] # To keep symmetry
          # }
        }
        
      }
    }
    diag(pvals) <- 0
    pvals <- as(pvals, "dspMatrix")
  }
  if(verbose){
    message("Building output...")
  }
  pvals <- as(pvals, "dspMatrix")
  return(pvals)
}


