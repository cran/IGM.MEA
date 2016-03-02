\name{axion.spikestodf}
\alias{axion.spikestodf}
\title{
axion.spikestodf
}

\description{
Diana needs to add description for this.
}
\usage{
axion.spikestodf(spikes)
}
\arguments{
  \item{spikes}{A list of spike trains, named by their electrode name. see example
}
}

\value{
returns a data frame with columns 'elec' and 'time'. One row for each spike in input train, with electode name and time of spike recorded. 
}

\seealso{
axion.elec.name.to.xy
}
\examples{
data("S") 
g<-axion.spikestodf(S$spikes )

}
\keyword{ axion }
\keyword{ electrode }% __ONLY ONE__ keyword per line
