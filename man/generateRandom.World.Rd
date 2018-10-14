%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  
% 
% on Wed Feb 08 14:37:48 2006.
% 
% Generator was the Rdoc class, which is part of the R.oo package written
% by Henrik Bengtsson, 2001-2004.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{generateRandom.World}
\alias{generateRandom.World}
\alias{World.generateRandom}
\alias{generateRandom.World}
\alias{generateRandom,World-method}
\keyword{methods}
\keyword{internal}


\title{Generates random values for all niches in the world}

\description{
 Generates random values for all niches in the world. It only pass the message \code{generateRandom} to all its niches.
}

\usage{\method{generateRandom}{World}(.O, ...)}

\value{
 Returns nothing.
}

\examples{
  cr <- Chromosome(genes=newCollection(Gene(shape1=1, shape2=100),5))
  ni <- Niche(chromosomes=newRandomCollection(cr, 10))
  wo <- World(niches=newRandomCollection(ni,2))
  wo
  wo2 <- wo
  generateRandom(wo2)
  wo2
  wo			# wo and wo2 are the very same object
  wo3 <- clone(wo2)
  generateRandom(wo3)
  wo3
  wo2			# now wo2 is different to wo3
  wo			# but wo2 is still the same than wo
}

\references{Goldberg, David E. 1989 \emph{Genetic Algorithms in Search, Optimization and Machine Learning}. Addison-Wesley Pub. Co. ISBN: 0201157675}

\author{Victor Trevino. Francesco Falciani Group. University of Birmingham, U.K. http://www.bip.bham.ac.uk/bioinf}

\seealso{
  For more information see \code{\link{World}}.
  \code{\link{unObject}},
  \code{\link{as.list.Object}},
  \code{\link{newCollection.World}},
  \code{\link{newRandomCollection.World}}
}

\keyword{methods}
