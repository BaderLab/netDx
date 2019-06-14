\name{silh}
\alias{silh}
\docType{data}
\title{
	Toy network.
}
\description{
List with two entries. 
net: Network specification. "X" and "Y" are source and target columns respectively. "DIST" specifies weights.
groups: Node labsls. A data.frame with columns "ID" and "GROUP"
}
\usage{data(silh)}
\examples{
data(silh)
summary(silh)
silh$net
silh$groups
}
\keyword{datasets}
