test_getMutexGroup <- function() {
  
  obs <- tryCatch(getMutexGroup(A = NULL, 
                           PM = matrix(sample(0:1, 25, replace=TRUE), ncol=5), 
                           type = "Impurity", 
                           lower.tail = TRUE), error=conditionMessage)
  checkIdentical("not input matrix A", obs)
  
  obs <- tryCatch(getMutexGroup(A = 12, 
                                PM = matrix(sample(0:1, 25, replace=TRUE), ncol=5), 
                                type = "Impurity", 
                                lower.tail = TRUE), error=conditionMessage)
  checkIdentical("input A must be a Matrix or a matrix class", obs)
  
  obs <- tryCatch(getMutexGroup(A = matrix(sample(0:1, 25, replace=TRUE), ncol=5), 
                                PM = NULL, 
                                type = "Impurity", 
                                lower.tail = TRUE), error=conditionMessage)
  checkIdentical("not probability matrix PM", obs)
  
  obs <- tryCatch(getMutexGroup(A = matrix(sample(0:1, 25, replace=TRUE), ncol=5), 
                                PM = 12, 
                                type = "Impurity", 
                                lower.tail = TRUE), error=conditionMessage)
  checkIdentical("input PM must be a Matrix or a matrix class", obs)
  
  obs <- tryCatch(getMutexGroup(A = matrix(NA,nrow = 10,ncol = 0), 
                           PM = matrix(sample(0:1, 25, replace=TRUE), ncol=5), 
                           type = "Impurity", 
                           lower.tail = TRUE), error=conditionMessage)
  checkIdentical("input A must have at least 1 row and 1 column", obs)
  
  obs <- tryCatch(getMutexGroup(A = matrix(sample(0:1, 25, replace=TRUE), ncol=5), 
                                PM = matrix(NA,nrow = 10,ncol = 0), 
                                type = "Impurity", 
                                lower.tail = TRUE), error=conditionMessage)
  checkIdentical("input PM must have at least 1 row and 1 column", obs)
  
  obs <- tryCatch(getMutexGroup(A = matrix(2,nrow = 10,ncol = 5), 
                                PM = matrix(sample(0:1, 25, replace=TRUE), ncol=5), 
                                type = "Impurity", 
                                lower.tail = TRUE), error=conditionMessage)
  checkIdentical("input A must be binary", obs)
  
  obs <- tryCatch(getMutexGroup(A = matrix(sample(0:1, 25, replace=TRUE), ncol=5), 
                                PM = matrix(2,nrow = 10,ncol = 5), 
                                type = "Impurity", 
                                lower.tail = TRUE), error=conditionMessage)
  checkIdentical("input PM must be binary", obs)
  
}

