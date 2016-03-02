\name{axion.spikesum2}
\alias{axion.spikesum2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
axion.spikesum2
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Returns a vector containing summary information about spike train input.
}
\usage{
axion.spikesum2(spikes)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{spikes}{
%%     ~~Describe \code{spikes} here~~
A spike train, or a list of vectors of spike trains
}
}
\author{
%%  ~~who you are~~
Diana Hall
}

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
axion.spikesum
}
\examples{
data("S") 
axion.spikesum2(S$spikes)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ spikes }

