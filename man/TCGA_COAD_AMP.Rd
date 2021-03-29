\name{TCGA_COAD_AMP}
\alias{TCGA_COAD_AMP}
\docType{data}
\title{TCGA_COAD_AMP data}
\description{
A binary matrix, with information about genes mutations in Colon Adenocarcinoma, created from TCGA_COAD.rda by applying \code{GDCquery} in order to obtain information about the amplifications, and used as real example in \code{getMutexAB}.
}
\usage{data("TCGA_COAD_AMP")}
\format{
  The format is:
  
 num [1:1000, 1:391] 1 1 1 1 1 1 1 1 1 1 ...
 
 - attr(*, "dimnames")=List of 2
 
  ..$ : chr [1:1000] "APC" "TP53" "TTN" "KRAS" ...
  
  ..$ : chr [1:391] "TCGA-CA-6718" "TCGA-D5-6931" "TCGA-AZ-6601" "TCGA-G4-6320" ...
}
\examples{
data(TCGA_COAD_AMP)
}
\keyword{datasets}
