test_getMutexAB <- function() {

  obs <- tryCatch(getMutexAB(A = NULL, 
                             PMA = getPM(A),
                             B = matrix(sample(0:1, 10, replace = TRUE), ncol = 5),
                             PMB = getPM(B),
                             lower.tail = TRUE, 
                             mixed = FALSE,
                             th = 1e-2, verbose = FALSE, parallel = FALSE), error=conditionMessage)
  checkIdentical("not input matrix A", obs)
  
  obs <- tryCatch(getMutexAB(A = matrix(sample(0:1, 10, replace = TRUE), ncol = 5), 
                             PMA = getPM(A),
                             B = NULL,
                             PMB = getPM(B),
                             lower.tail = TRUE, 
                             mixed = FALSE,
                             th = 1e-2, verbose = FALSE, parallel = FALSE), error=conditionMessage)
  checkIdentical("not input matrix B", obs)
  
  obs <- tryCatch(getMutexAB(A = 12, 
                             PMA = getPM(A),
                             B = matrix(sample(0:1, 10, replace = TRUE), ncol = 5),
                             PMB = getPM(B),
                             lower.tail = TRUE, 
                             mixed = FALSE,
                             th = 1e-2, verbose = FALSE, parallel = FALSE), error=conditionMessage)
  checkIdentical("input A must be a Matrix or a matrix class", obs)
  
  obs <- tryCatch(getMutexAB(A = matrix(sample(0:1, 10, replace = TRUE), ncol = 5), 
                             PMA = getPM(A),
                             B = 30,
                             PMB = getPM(B),
                             lower.tail = TRUE, 
                             mixed = FALSE,
                             th = 1e-2, verbose = FALSE, parallel = FALSE), error=conditionMessage)
  checkIdentical("input B must be a Matrix or a matrix class", obs)
  
  obs <- tryCatch(getMutexAB(A = matrix(NA,nrow = 10,ncol = 0), 
                             PMA = getPM(A),
                             B = matrix(sample(0:1, 10, replace = TRUE), ncol = 5),
                             PMB = getPM(B),
                             lower.tail = TRUE, 
                             mixed = FALSE,
                             th = 1e-2, verbose = FALSE, parallel = FALSE), error=conditionMessage)
  checkIdentical("input A must have at least 1 row and 1 column", obs)
  
  obs <- tryCatch(getMutexAB(A = matrix(sample(0:1, 10, replace = TRUE), ncol = 5), 
                             PMA = getPM(A),
                             B = matrix(NA,nrow = 10,ncol = 0),
                             PMB = getPM(B),
                             lower.tail = TRUE, 
                             mixed = FALSE,
                             th = 1e-2, verbose = FALSE, parallel = FALSE), error=conditionMessage)
  checkIdentical("input B must have at least 1 row and 1 column", obs)
  
  obs <- tryCatch(getMutexAB(A = matrix(2,nrow = 10,ncol = 5), 
                             PMA = getPM(A),
                             B = matrix(sample(0:1, 10, replace = TRUE), ncol = 5),
                             PMB = getPM(B),
                             lower.tail = TRUE, 
                             mixed = FALSE,
                             th = 1e-2, verbose = FALSE, parallel = FALSE), error=conditionMessage)
  checkIdentical("input A must be binary", obs)
  
  obs <- tryCatch(getMutexAB(A = matrix(sample(0:1, 10, replace = TRUE), ncol = 5), 
                             PMA = getPM(A),
                             B = matrix(2,nrow = 10,ncol = 5),
                             PMB = getPM(B),
                             lower.tail = TRUE, 
                             mixed = FALSE,
                             th = 1e-2, verbose = FALSE, parallel = FALSE), error=conditionMessage)
  checkIdentical("input B must be binary", obs)
  
}

