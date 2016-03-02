\name{make.raster}
\alias{make.raster}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
make.raster
}
\description{
Creates a pdf raster plot of selected user selected well from an 's' object. Options include verticle lines showing network spike times, vertical bars showing bursts as well as # showing count of spikes in burst and network spikes.
}
\usage{
make.raster(RobjectFile = NULL, outputdir = NULL, well.for.raster
                 = NULL, interval.for.raster = NULL, show.bursts = F,
                 show.burst.number = F, show.networkspikes = F,
                 show.ns.number = F)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{RobjectFile}{
Default value is NULL, in which case tcltk pop-up file chooser will prompt user to select an 's' object.  Otherwise, provide a full path to to a .RData 's' object that contains burst and network data.
}

  \item{outputdir}{
   A directory (character string in quotes) where pdf is to be saved.
   Default is NULL, in which case the plot will be saved in the directory RobjectFile location.
}
  \item{well.for.raster}{
   A well name, character string, from plate. e.g. well.for.raster="A3".
   Default is NULL, in which case first well in plate will appear in plot.
}
  \item{interval.for.raster}{
     A vector of min and max time (s) for raster marks. e.g. interval.for.raster=c(30,60)
     Default is NULL, in which case the whole recording interval will be used.
}
  \item{show.bursts}{
  A boolean value sets whether bursts are indicated by red horizontal line (TRUE/FALSE)
  e.g. show.bursts=FALSE
  Default=FALSE
}
  \item{show.burst.number}{
 A boolean value sets whether # spikes/bursts are indicated (TRUE/FALSE).
 show.bursts must be set to true in order that show.burst.number=T
 e.g. show.burst.number=FALSE
 Default=FALSE
}
  \item{show.networkspikes}{
A boolean value sets whether network spikes are indicated by green vertical line (TRUE/FALSE)
  e.g. show.networkspikes=FALSE
   Default=FALSE
}
  \item{show.ns.number}{
A boolean value sets whether # electrodes in network spikes are indicated (TRUE/FALSE)
  e.g. show.ns.number=FALSE
   Default=FALSE
}

}
\value{
A pdf raster plot will be displayed in system viewer.
}
\author{
%%  ~~who you are~~
Diana Hall
}
\examples{
##make.raster(RobjectFile=NULL,
#                      well.for.raster=NULL, 
#                      interval.for.raster=NULL,
#                      show.bursts=F, 
#                      show.burst.number=F, 
#                      show.networkspikes=F,
#                      show.ns.number=F )

}
\keyword{ utility }
\keyword{ spikes }% __ONLY ONE__ keyword per line
