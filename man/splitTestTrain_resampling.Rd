% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/splitTestTrain_partition.R
\name{splitTestTrain_resampling}
\alias{splitTestTrain_resampling}
\title{Assign train/test labels over several resamplings of the data.}
\usage{
splitTestTrain_resampling(pheno_DF, nFold = 3L, predClass, verbose = FALSE)
}
\arguments{
\item{pheno_DF}{(data.frame) table with patient ID and status.
Must contain columns for Patient ID (named 'ID') and class
(named 'STATUS'). Status should be a char; value of predictor class 
should be specified in \code{predClass} param; 
all other values are considered non-predictor class
Expects rows with unique IDs
Rows with duplicate IDs will be excluded.}

\item{nFold}{(integer) number of resamplings. Each sample will be a test
sample exactly once.}

\item{predClass}{(char) name of predictor class}

\item{verbose}{(logical) print messages}
}
\value{
(list) of length nFold, each with char vector of length 
nrow(pheno_DF). Values of 'TRAIN' or 'TEST'
}
\description{
Assign train/test labels over several resamplings of the data.
}
\details{
This function is useful when feature selection needs to 
occur over multiple resamplings of the data, as a strategy to reduce
overfitting. Each sample serves as a test for exactly one resampilng, 
and as a training sample for the others. The method is provided with the
positive label and splits the samples so that an even number of positive
and negative classes are represented in all the resamplings (i.e. it
avoids the situation where one resampling has too many positives and 
another has too few).
}
\examples{
data(pheno) 
x <- splitTestTrain_resampling(pheno,predClass='LumA')
}
