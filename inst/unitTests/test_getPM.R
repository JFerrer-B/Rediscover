test_getPM <- function() {
  
  obs <- tryCatch(getPM(A = NULL), 
                  error=conditionMessage)
  checkIdentical("not input matrix A", obs)
  
  obs <- tryCatch(getPM(A = 12), 
                  error=conditionMessage)
  checkIdentical("input A must be a Matrix or a matrix class", obs)
  
  obs <- tryCatch(getPM(A = matrix(2,nrow = 10,ncol = 5)), 
                  error=conditionMessage)
  checkIdentical("input A must be binary", obs)
  
  obs <- tryCatch(getPM(A = matrix(NA,nrow = 10,ncol = 0)), 
                  error=conditionMessage)
  checkIdentical("input A must have at least 1 row and 1 column", obs)

}