\name{axion.guess.well.number}
\alias{axion.guess.well.number}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
axion.guess.well.number
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Given a plate recording, guess whether the plate contained 12 or 48 wells.
}
\usage{
axion.guess.well.number(channels)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{channels}{A list of channels names. e.g. c("A3_12","B4_14")}
}

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
axion.elec2well
plateinfo
}
\examples{

data("S") 
well.number <- axion.guess.well.number(S$channels)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ layout }

