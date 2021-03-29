\name{PM_COAD}
\alias{PM_COAD}
\docType{data}
\title{PM_COAD data}
\description{
Probability matrix, with information of genes being mutated in samples in Colon Adenocarcinoma, created by TCGA_COAD.rda applying \code{getPM} and used as real example in \code{getMutex} and \code{getMutexAB} and \code{getMutexGroup}.
}
\usage{data("PM_COAD")}
\format{
  The format is:
  
Formal class 'PMatrix' [package "Rediscover"] with 2 slots

  ..@ rowExps: num [1:399] 13.1 1.02 7.43 3.26 0.4 ...
  
  ..@ colExps: num [1:17616] 2.54 1.78 1.76 1.35 0.6 ...
}
\examples{
data(PM_COAD)
}
\keyword{datasets}
