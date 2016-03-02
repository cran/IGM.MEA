\name{plateinfo}
\alias{plateinfo}
\title{
Retrieve plate information}
\description{
Retrieves layout parameters from internal variable according to whether Axion plate contains 12 or 48 wells.
}
\usage{
plateinfo(arrayname)
}
\arguments{
  \item{arrayname}{
arrayname must be either:
  
`Axion 12 well`

`Axion 48 well`

according to whether plate is 12 or 48 wells.
}
}
\value{
\item{n.well }{# wells on plate}
\item{wells }{well names}
\item{n.well.r }{# rows of wells per plate}
\item{n.well.c}{# columns of wells per plate}
\item{layout }{vector of # rows by # columns of wells}
\item{n.elec.r }{# rows of electrodes per well }
\item{n.elec.c}{# columns of electrodes per well }
}
\author{
Diana Hall
}
\note{
axion.guess.well.number can be used to automatically generate the probable # wells on a plate. see axion.guess.well.number function documentation for further details
}

\seealso{
See Also axion.guess.well.number

}
\examples{
data("S") # load data
wells<-axion.guess.well.number(S$channels) #get # wells
array <- paste0("Axion ", wells, " well") #create array input
plateinfo <- plateinfo(array) # get plateinfo

}
\keyword{ plate }
