\name{A_Matrix}
\alias{A_Matrix}
\docType{data}
\title{A_Matrix data}
\description{
A binary dgCMatrix matrix used as toy example in \code{getPM} and \code{getMutex} and \code{getMutexAB} and \code{getMutexGroup}
}
\usage{data("A_Matrix")}
\format{
  The format is:
  
Formal class 'dgCMatrix' [package "Matrix"] with 6 slots

  ..@ i       : int [1:249838] 5 9 10 11 13 14 18 20 23 24 ...
  
  ..@ p       : int [1:501] 0 503 1010 1506 1995 2497 2981 3488 4002 4474 ...
  
  ..@ Dim     : int [1:2] 1000 500
  
  ..@ Dimnames:List of 2
  
  .. ..$ : NULL
  
  .. ..$ : NULL
  
  ..@ x       : num [1:249838] 1 1 1 1 1 1 1 1 1 1 ...
  
  ..@ factors : list()
}
\examples{
data(A_Matrix)
}
\keyword{datasets}
