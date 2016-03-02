\name{axion.spikesum}
\alias{axion.spikesum}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
axion.spikesum
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Provides a brief description of spike trainr.
}
\usage{
axion.spikesum(spikes)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{spikes}{
%%     ~~Describe \code{spikes} here~~
a spike train. see example
}
}

\value{
Returns # electodes, # spikes, first and last spike.
}
\author{
%%  ~~who you are~~
Diana Hall
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
axion.spikesum2
plateinfo
axion.guess.well.number
}
\examples{
data("S")
axion.spikesum(S$spikes)
}


