\name{predRes}
\alias{predRes}
\docType{data}
\title{
	Example output of getPatientRankings, used to call labels for test patients.
}
\description{
List of lists. First level is a list of size 4, with one key entry for each tumour type in example medulloblastoma dataset (WNT,SHH,Group3,Group4). Each list in the second level is of length 8, with structure corresponding to the output of getPatientRankings().
}
\usage{data(predRes)}
\examples{
data(predRes)
summary(predRes)
summary(predRes[[1]])
}
\keyword{datasets}
