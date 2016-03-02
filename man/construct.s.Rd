\name{construct.s}
\alias{construct.s}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
construct.s
}
\description{
Creates an 's' object of class 'mm.s', an RData list, from list of spike trains and channels names.

}
\usage{
construct.s(spikes, ids, time.interval, beg, end, corr.breaks, layout, filename)
}
\arguments{
  \item{spikes}{
named list of spike trains, names are given by electrodes e.g. "E4_33"
}
  \item{ids}{
  NULL: all input spike trains in 'spikes' will appear in output 's' object
  vector of channel names: select subset of channels to appear in output 's' object
}
  \item{time.interval}{
bin length (s) to partition the firing rate. e.g. time.interval=2 results in a 20 second recording having 10 bins of firing rate data, occuring at time 0.01, 2.01, 4.01s etc
}
  \item{beg}{
NULL: include all spikes
a number n: filter out spikes that occured before n from all spike trains

}
  \item{end}{
NULL: include all spikes
a number n: filter out spikes that occured after n from all spike trains

}
  \item{corr.breaks}{
0
}
  \item{layout}{
info about plate layout as compiled by get.array.info. See get.array.info for further details.
}
  \item{filename}{
name of recording file
}
}
\value{
\item{channels }{Description of 'comp2'}
\item{spikes }{ list of spike trains, named by channel}
\item{nspikes }{# spikes on each train}
\item{NCells }{# spike trains}
\item{layout}{ the input argument 'layout' }
\item{meanfiringrate }{mean firing rate for each channel}
\item{file }{recording file name }
\item{rates }{a list that included binned firing rate for each spike train}
\item{rec.time }{ vector c('beg', 'end') according to input parameters}

}

\seealso{
get.file.basename
plateinfo
axion.elec.name.to.xy
get.array.info
}
\examples{

data("S") 
wells <- axion.guess.well.number(S$channels)
array <- paste0("Axion ", wells, " well") #make array info
S$array<-array #store in S object
plateinfo <- plateinfo(S$array)

epos <- axion.elec.name.to.xy(S$channels, plateinfo) #create spacing
S$epos<-epos
S$layout <- get.array.info(S)

b<-construct.s(spikes=S$spikes, 
  ids=NULL, 
  time.interval=1,
  beg=NULL, end=NULL,
  corr.breaks=0,
  layout=S$layout,
  filename=get.file.basename(S$file) )
  
  
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ spike-list }
