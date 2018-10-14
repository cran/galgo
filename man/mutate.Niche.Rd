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

 \name{mutate.Niche}
\alias{mutate.Niche}
\alias{Niche.mutate}
\alias{mutate.Niche}
\alias{mutate,Niche-method}
\keyword{methods}
\keyword{internal}


 \title{Mutates a niche calling mutate method for all chromosomes}

 \description{
  Mutates a niche calling mutate method for all chromosomes.
 }

 \usage{\method{mutate}{Niche}(ni, n=(ni$mutationsFunc)(ni), ...)}

 \arguments{
   \item{n}{Number of chromosomes to mutate. The default is the result of calling \code{mutationsFunc}.}
 }

 \details{
	This method update the gene values for random chromsomes. The number of chromosomes to mutate is normally obtained calling \code{mutationFunc}.
 }

 \value{
  This methods returns the chromosome indexes mutated.
 }

 \examples{
   cr <- Chromosome(genes=newCollection(Gene(shape1=1, shape2=100),5))
   cr
   ni <- Niche(chromosomes = newRandomCollection(cr, 10))
   ni
   mutate(ni, 3)
   ni
 }

 \references{Goldberg, David E. 1989 \emph{Genetic Algorithms in Search, Optimization and Machine Learning}. Addison-Wesley Pub. Co. ISBN: 0201157675}

 \author{Victor Trevino. Francesco Falciani Group. University of Birmingham, U.K. http://www.bip.bham.ac.uk/bioinf}

 \seealso{
   For more information see \code{\link{Niche}}.
   \code{\link{mutate.Chromosome}}(),
   \code{\link{mutate.Gene}}().
 }

 \keyword{methods}
