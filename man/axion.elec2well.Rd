\name{axion.elec2well}
\alias{axion.elec2well}
\title{
Extract well name from electrode name
}
\description{
Takes the first two characters of electrode name and creates well name.
"E7_22" becomes "E7"
}
\usage{
axion.elec2well(elec)
}
\arguments{
  \item{elec}{
A vector of electrode names  e.g. c("E7_22", "E6_31")
}
}
\examples{
data("S") 
elec<-S$channels
wells<-axion.elec2well(elec)
wells[1:5] 
}
