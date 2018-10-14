%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  
% 
% on Wed Feb 08 14:37:34 2006.
% 
% Generator was the Rdoc class, which is part of the R.oo package written
% by Henrik Bengtsson, 2001-2004.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{generateRandom.Chromosome}
\alias{generateRandom.Chromosome}
\alias{Chromosome.generateRandom}
\alias{generateRandom.Chromosome}
\alias{generateRandom,Chromosome-method}
\keyword{methods}
\keyword{internal}


\title{Generates random values for all genes in the chromosome}

\description{
 Updates gene values calling the method \code{generateRandom} to all its genes.
}

\usage{\method{generateRandom}{Chromosome}(.O, ...)}

\value{
 Returns nothing.
}

\examples{
  cr <- Chromosome(genes=newCollection(Gene(shape1=1, shape2=100),5))
  cr
  #
  # from clone example
  cr2 <- cr
  generateRandom(cr2)
  cr2
  cr			# cr and cr2 are the very same object
  cr3 <- clone(cr2)
  generateRandom(cr3)
  cr3
  cr2			# now cr2 is different to cr3
  cr			# but cr2 is still the same than cr
}

\references{Goldberg, David E. 1989 \emph{Genetic Algorithms in Search, Optimization and Machine Learning}. Addison-Wesley Pub. Co. ISBN: 0201157675}

\author{Victor Trevino. Francesco Falciani Group. University of Birmingham, U.K. http://www.bip.bham.ac.uk/bioinf}

\seealso{
  For more information see \code{\link{Chromosome}}.
  \code{\link{unObject}},
  \code{\link{as.list.Object}},
  \code{\link{newCollection.Chromosome}},
  \code{\link{newRandomCollection.Chromosome}}
}

\keyword{methods}
