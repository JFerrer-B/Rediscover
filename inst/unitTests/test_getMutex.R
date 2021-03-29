test_getMutex <- function() {
  
  obs <- tryCatch(getMutex(A = NULL, 
                           PM = getPM(A), 
                           lower.tail = TRUE, 
                           mixed = FALSE,
                           th = 1e-2, verbose = FALSE, parallel = FALSE), error=conditionMessage)
  checkIdentical("not input matrix A", obs)
  
  obs <- tryCatch(getMutex(A = 12, 
                           PM = getPM(A), 
                           lower.tail = TRUE, 
                           mixed = FALSE,
                           th = 1e-2, verbose = FALSE, parallel = FALSE), error=conditionMessage)
  checkIdentical("input A must be a Matrix or a matrix class", obs)
  
  obs <- tryCatch(getMutex(A = matrix(2,nrow = 10,ncol = 5), 
                           PM = getPM(A), 
                           lower.tail = TRUE, 
                           mixed = FALSE,
                           th = 1e-2, verbose = FALSE, parallel = FALSE), error=conditionMessage)
  checkIdentical("input A must be binary", obs)
  
  obs <- tryCatch(getMutex(A = matrix(NA,nrow = 10,ncol = 0), 
                           PM = getPM(A), 
                           lower.tail = TRUE, 
                           mixed = FALSE,
                           th = 1e-2, verbose = FALSE, parallel = FALSE), error=conditionMessage)
  checkIdentical("input A must have at least 1 row and 1 column", obs)
  
}