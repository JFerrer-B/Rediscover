\docType{class}
\name{PMatrix-class}
\alias{PMatrix-class}
\alias{PMatrix}
\alias{[,PMatrix,matrix,missing,ANY-method}
\alias{[,PMatrix,missing,numeric,ANY-method}
\alias{[,PMatrix,numeric,missing,ANY-method}
\alias{[,PMatrix,numeric,numeric,ANY-method}
\alias{dim,PMatrix-method}
\title{An S4 class to store the probabilities}
\description{
An S4 class to store the probabilities of gene i being mutated in sample j
}
\section{Slots}{
\describe{
\item{\code{rowExps}}{Sample depending estimated coefficients obtained from the logistic regression}
\item{\code{colExps}}{gene depending estimated coefficients obtained from the logistic regression}
}}