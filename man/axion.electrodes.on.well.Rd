\name{axion.electrodes.on.well}
\alias{axion.electrodes.on.well}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
axion.electrodes.on.well
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Retrieve the set of channels belonging to the input well and list of channels.
}
\usage{
axion.electrodes.on.well(well, electrodes)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{well}{
%%     ~~Describe \code{well} here~~
Well name. e.g. "B2" or "C4"
}
  \item{electrodes}{
%%     ~~Describe \code{electrodes} here~~
A vector of electrode names. 
}
}

\value{
Returns a vector of channel names
}
\author{
%%  ~~who you are~~
Diana Hall
}

\seealso{

axion.elec2well
}
\examples{
data("S") 
well<-"B2"
axion.electrodes.on.well(well, S$channels)

}

