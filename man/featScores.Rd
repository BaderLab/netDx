\name{featScores}
\alias{featScores}
\docType{data}
\title{
Demo feature-level scores from running feature selection on two-class problem
}
\description{
List with one entry per patient label ("SURVIVEYES" and "SURVIVENO"). Each entry contains scores obtained through feature-selection acros 100 train/test splits. Scores range from 0 to 10. Scores in data.frame format, with rows corresponding to features and columns to a particular train/test split.
}
\usage{data(featScores)}
\examples{
data(featScores)
head(featScores)
}
\keyword{datasets}
