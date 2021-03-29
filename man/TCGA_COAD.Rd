\name{TCGA_COAD}
\alias{TCGA_COAD}
\docType{data}
\title{TCGA_COAD}
\description{
A binary matrix, with information about genes mutations in Colon Adenocarcinoma, created by applying \code{maftools} to .maf file and used as real example in \code{getPM} and \code{getMutex} and \code{getMutexAB} and \code{getMutexGroup}.
}
\usage{data("TCGA_COAD")}
\format{
  The format is:
  
 num [1:399, 1:17616] 1 1 1 1 1 1 1 1 1 1 ...
 
 - attr(*, "dimnames")=List of 2
 
  ..$ : chr [1:399] "TCGA-CA-6718" "TCGA-D5-6931" "TCGA-AZ-6601" "TCGA-G4-6320" ...
  
  ..$ : chr [1:17616] "APC" "TP53" "TTN" "KRAS" ...
}
\examples{
data(TCGA_COAD)
}
\keyword{datasets}
