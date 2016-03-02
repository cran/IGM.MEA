\name{axion.elec.name.to.xy}
\alias{axion.elec.name.to.xy}

\title{
%%  ~~function to do ... ~~
axion.elec.name.to.xy
}

\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Convert 1 or multple electrode names into the (x,y) position where it would be found on a plate.
}
\usage{
axion.elec.name.to.xy(name, plateinfo)
}
\arguments{
  \item{name}{one or multiple channel names e.g. "E7_22".  see example}
  \item{plateinfo}{
  a list with\describe{ 
  \item{n.well}{# wells in recording}
  \item{wells}{ well names}
  \item{n.well.r}{ # well rows in plate layout}
  \item{n.well.c}{# well columns in plate layout}
  \item{layout}{ vector, first entry # columns, 2nd entry # rows}
  \item{n.elec.r}{ # rows of electrodes/ well}
  \item{n.elec.c}{ # columns of electrodes/ well}
  }
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
Used to create coordinates for each electrode on a plate for plotting purposes.

}
\value{

returns a matrix with number of rows equal to number of channels in data and 2 columns.
Each row corresponds to an 'x' and 'y' coordinates of the channel location on a plate.

}

\author{
%%  ~~who you are~~
Diana Hall
}
\note{
%%  ~~further notes~~
Useful for 12-well and 48-well plates.

}


\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
plateinfo
axion.guess.well.number
}
\examples{
# load data
data('S') 
wells <- axion.guess.well.number(S$channels)
array <- paste0("Axion ", wells, " well")
plateinfo <- plateinfo(array)
epos <- axion.elec.name.to.xy(S$channels, plateinfo)
##plot(epos)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ compute.ns }

