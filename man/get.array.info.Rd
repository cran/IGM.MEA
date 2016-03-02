\name{get.array.info}
\alias{get.array.info}
\title{
get.array.info
}
\description{
Retrieves data previously created in spike object related to layout of plate.
}
\usage{
get.array.info(data)
}
\arguments{
  \item{data}{
A spike list object that must have a field named 'array', & a field named 'epos'.
}
}
\value{
creates a object of class "mealayout" 
\item{layout }{created by platelayout function}
\item{xlim }{created by axion.elec.name.to.xy function}
\item{ylim }{created by axion.elec.name.to.xy function }
\item{spacing }{space between electrode grid in plot layout}
\item{pos }{created by axion.elec.name.to.xy function }
\item{pos }{created by axion.elec.name.to.xy function }
\item{corr.breaks }{ inherited }
}

\seealso{
axion.guess.well.number
axion.elec.name.to.xy
plateinfo
}
\examples{

data("S") 
wells <- axion.guess.well.number(S$channels)
array <- paste0("Axion ", wells, " well") #make array info
S$array<-array #store in S object
plateinfo <- plateinfo(S$array)
epos <- axion.elec.name.to.xy(S$channels, plateinfo) #create spacing
S$epos<-epos
arrayinfo <- get.array.info(S)
}
\keyword{ s-object }

