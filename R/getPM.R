#' getPM function
#' 
#' Given a binary matrix estimates the corresponding probability  matrix pij.
#'
#' @param A The binary matrix
#'
#' @return A `PMatrix` object with the corresponding probability estimations.
#' This `PMatrix` object stored the corresponding coefficients of the 
#' logistic regression computed. With this coefficients it is possible to build
#' the complete matrix of probabilities.
#'
#' @examples 
#' 
#' 
#'   #This first example is a basic example of how to perform getPM: 
#'   
#'   data("A_example")
#'   PMA <- getPM(A_example)
#'   
#'   #The next example, is the same as the first one but, 
#'   #using a matrix of class Matrix: 
#'   
#'   data("A_Matrix")
#'   PMA_Matrix <- getPM(A_Matrix)
#'   
#'   #Finally, the last example, shows a real example 
#'   #of how to perform this function when when using
#'   #data from TCGA, Colon Adenocarcinoma in this case: 
#'   data("TCGA_COAD")
#'   PM_COAD <- getPM(TCGA_COAD)
#'
#'
#' @import Matrix
#' @importFrom stats coef binomial gaussian dnorm pnorm
#' @importFrom speedglm control is.sparse
#' @export


getPM <- function(A){
  
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
  
  # rSc <- rowSums(A)
  rSc <- as.numeric(A %*% matrix(1,nrow = ncol(A)))
  
  cSc <- colSums(A)
  rSums <- factor(rSc)
  cSums <- factor(cSc)

  Xbig2small <- sparseMatrix(i = 1:(nrow(A)+ncol(A)),
                             j = c(as.numeric(rSums), as.numeric(cSums)+length(levels(rSums))),
                             x =1)

  rM <- sparseMatrix(i = as.numeric(rSums),
                     j = 1:nrow(A),
                     x =1)
  cM <- sparseMatrix(i = 1:ncol(A),
                     j = as.numeric(cSums),
                     x =1)

  A1 <- rM %*% (A %*% cM) 
  #A0 <- rM %*% (1-A) %*% cM 
  # A0 <- rM %*% (1-A) %*% cM -A1 
  a <- rowSums(rM)
  b <- colSums(cM)
  A0 <- matrix(a,,1) %*% matrix(b,1) - A1

  X1 <- sparseMatrix(i = 1:(nrow(A1)*ncol(A1)),
                     j = rep(1:nrow(A1), ncol(A1)),
                     x = 1,
                     dims = c((nrow(A1)*ncol(A1)), (nrow(A1)+ncol(A1))))
  
  X2 <- sparseMatrix(i = 1:(nrow(A1)*ncol(A1)),
                     j = nrow(A1)+rep(1:ncol(A1), each = nrow(A1)),
                     x = 1,
                     dims = c((nrow(A1)*ncol(A1)), (nrow(A1)+ncol(A1))))
  Xtxiki <- X1 + X2 # Matriz de diseño
  Xtxiki <- rbind(Xtxiki, c(rep(1,nrow(A1)), rep(0,ncol(A1))))
  Salida5 <- speedglm.wfit2(X=(Xtxiki),
                            y=cbind(c(as.vector(A1),round(mean(A1))), c(as.vector(A0),round(mean(A0)))),
                            family = binomial(), sparse=TRUE)
  
  cS5A <- coef(Salida5)
  # cS5 <- cS5A %*% t(Xbig2small)
  # elog5 <- matrix(exp(cS5[1:nrow(A)]), ncol = 1) %*% matrix(exp(cS5[nrow(A) + 1:ncol(A)]), nrow = 1)
  # p5 <- elog5 /(1+elog5)
  
  expcS5A <- exp(cS5A)
  expcS5 <- expcS5A %*% t(Xbig2small)
  
  p5 <- PMatrix(rowExps = expcS5[1:nrow(A)], colExps = expcS5[nrow(A) + (1:ncol(A))])
  return(p5)
}

