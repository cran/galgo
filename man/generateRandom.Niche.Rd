%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  
% 
% on Wed Feb 08 14:37:38 2006.
% 
% Generator was the Rdoc class, which is part of the R.oo package written
% by Henrik Bengtsson, 2001-2004.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{generateRandom.Niche}
\alias{generateRandom.Niche}
\alias{Niche.generateRandom}
\alias{generateRandom.Niche}
\alias{generateRandom,Niche-method}
\keyword{methods}
\keyword{internal}


\title{Generates random values for all genes contained in all chromosomes in the niche}

\description{
 It only pass the message \code{generateRandom} to all its chromosomes.
}

\usage{\method{generateRandom}{Niche}(.O, ...)}

\value{
 Returns nothing.
}

\examples{
  cr <- Chromosome(genes=newCollection(Gene(shape1=1, shape2=100),5))
  ni <- Niche(chromosomes=newRandomCollection(cr, 10))
  ni
  ni2 <- ni
  generateRandom(ni2)
  ni2
  ni			# ni and ni2 are the very same object
  ni3 <- clone(ni2)
  generateRandom(ni3)
  ni3
  ni2			# now cr2 is different to cr3
  ni			# but cr2 is still the same than cr
}

\references{Goldberg, David E. 1989 \emph{Genetic Algorithms in Search, Optimization and Machine Learning}. Addison-Wesley Pub. Co. ISBN: 0201157675}

\author{Victor Trevino. Francesco Falciani Group. University of Birmingham, U.K. http://www.bip.bham.ac.uk/bioinf}

\seealso{
  For more information see \code{\link{Niche}}.
  \code{\link{unObject}},
  \code{\link{as.list.Object}},
  \code{\link{newCollection.Niche}},
  \code{\link{newRandomCollection.Niche}}.
}

\keyword{methods}
