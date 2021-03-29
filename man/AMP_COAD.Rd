\name{AMP_COAD}
\alias{AMP_COAD}
\docType{data}
\title{AMP_COAD data}
\description{
A binary matrix, with information about amplifications in Colon Adenocarcinoma, created by applying \code{GDCquery} and used as real example in \code{getMutexAB}.
}
\usage{data("AMP_COAD")}
\format{
  The format is:
  
 num [1:1000, 1:391] 0 0 0 0 0 0 0 0 0 0 ...
 
 - attr(*, "dimnames")=List of 2
 
  ..$ : chr [1:1000] "ENSG00000212993.4" "ENSG00000279524.1" "ENSG00000136997.13" "ENSG00000101294.15" ...
  
  ..$ : chr [1:391] "TCGA-CA-6718" "TCGA-D5-6931" "TCGA-AZ-6601" "TCGA-G4-6320" ...
}
\examples{
data(AMP_COAD)
## maybe str(AMP_COAD) ; plot(AMP_COAD) ...
}
\keyword{datasets}
