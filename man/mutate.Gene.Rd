%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  
% 
% on Wed Feb 08 14:37:29 2006.
% 
% Generator was the Rdoc class, which is part of the R.oo package written
% by Henrik Bengtsson, 2001-2004.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{mutate.Gene}
\alias{mutate.Gene}
\alias{Gene.mutate}
\alias{mutate.Gene}
\alias{mutate,Gene-method}
\alias{mutate}
\keyword{methods}
\keyword{internal}


\title{Mutates a gene}

\description{
 Mutate a gene. This method is the proper way to call generateRandom method. For better description \code{\link[galgo:generateRandom.Gene]{*generateRandom}()}.
}

\usage{\method{mutate}{Gene}(.O, ...)}

\arguments{
  \item{n}{Number of random values.}
}

\value{
 Returns random values.
}

\examples{
  ge <- Gene(shape1=1, shape2=100)
  ge
  mutate(ge)
  mutate(ge)
  mutate(ge)
}

\references{Goldberg, David E. 1989 \emph{Genetic Algorithms in Search, Optimization and Machine Learning}. Addison-Wesley Pub. Co. ISBN: 0201157675}

\author{Victor Trevino. Francesco Falciani Group. University of Birmingham, U.K. http://www.bip.bham.ac.uk/bioinf}

\seealso{
  For more information see \code{\link{Gene}}.
  \code{\link[galgo:generateRandom.Gene]{*generateRandom}()}.
}

\keyword{methods}
