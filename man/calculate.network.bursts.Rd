\name{calculate.network.bursts}
\alias{calculate.network.bursts}

\title{
Compute network bursts for a list of MEA recordings.
}
\description{
For a list of MEA recordings, ususally from the same plate at different time point,This function detects and report network burst features at the well level.
}
\usage{
calculate.network.bursts(s,Sigma, min_electrodes, local_region_min_nAE)
}

\arguments{
  \item{s}{ A list of MEA recordings, typically from the same MEA plate at different time point.}
  \item{Sigma}{The window size used to generate network bursts.}
  \item{min_electrodes}{ Minimum number of electrodes to call a network burst}
  \item{local_region_min_nAE}{Indicates if an adaptive threthold method should be used.}
}

\value{
Returns an object containing all the network features for all wells.
}
\references{
Add reference to Yi-Fan Lu's paper when it is in press.
}
\author{
Quanli Wang
}


\keyword{ network bursts }
