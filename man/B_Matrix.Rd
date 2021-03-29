\name{B_Matrix}
\alias{B_Matrix}
\docType{data}
\title{B_Matrix data}
\description{
A binary dgCMatrix matrix used as toy example in \code{getPM} and \code{getMutex} and \code{getMutexAB} and \code{getMutexGroup}. 
}
\usage{data("B_Matrix")}
\format{
  The format is:
  
Formal class 'dgCMatrix' [package "Matrix"] with 6 slots

  ..@ i       : int [1:249526] 1 2 4 8 9 11 13 15 18 20 ...
  
  ..@ p       : int [1:501] 0 498 1014 1527 2048 2558 3036 3511 4035 4537 ...
  
  ..@ Dim     : int [1:2] 1000 500
  
  ..@ Dimnames:List of 2
  
  .. ..$ : NULL
  
  .. ..$ : NULL
  
  ..@ x       : num [1:249526] 1 1 1 1 1 1 1 1 1 1 ...
  
  ..@ factors : list()
}
\examples{
data(B_Matrix)
}
\keyword{datasets}
