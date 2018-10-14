%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  
% 
% on Wed Feb 08 14:37:42 2006.
% 
% Generator was the Rdoc class, which is part of the R.oo package written
% by Henrik Bengtsson, 2001-2004.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{newCollection.Niche}
\alias{newCollection.Niche}
\alias{Niche.newCollection}
\alias{newCollection.Niche}
\alias{newCollection,Niche-method}
\keyword{methods}
\keyword{internal}


\title{Generates a list of cloned niches}

\description{
 Generates a list of cloned niches. It only use the generic newCollection method.
}

\usage{\method{newCollection}{Niche}(.O, ...)}

\arguments{
  \item{n}{Number of object clones.}
}

\value{
 Returns a list with cloned objects. The names are build with the class and a consecutive number.
}

\examples{
  cr <- Chromosome(genes=newCollection(Gene(shape1=1, shape2=100),5))
  cr
  ni <- Niche(chromosomes=newRandomCollection(cr, 10))
  ni
  newCollection(ni, 2)                  # list of two new identical Niche objects
  newRandomCollection(ni, 2)            # list of two new different Niche objects
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
