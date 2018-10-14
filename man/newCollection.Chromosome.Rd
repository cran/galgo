%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  
% 
% on Wed Feb 08 14:37:36 2006.
% 
% Generator was the Rdoc class, which is part of the R.oo package written
% by Henrik Bengtsson, 2001-2004.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{newCollection.Chromosome}
\alias{newCollection.Chromosome}
\alias{Chromosome.newCollection}
\alias{newCollection.Chromosome}
\alias{newCollection,Chromosome-method}
\keyword{methods}
\keyword{internal}


 \title{Generates a list of chromosomes cloning the original chromosome object}

 \description{
 Generates a list of chromosomes cloning the original chromosome object. It only use the generic newCollection method.
 }

 \usage{\method{newCollection}{Chromosome}(.O, ...)}

 \arguments{
   \item{n}{Number of object clones.}
 }

 \value{
  Returns a list with cloned objects. The names are build with the class and a consecutive number.
 }

 \examples{
   cr <- Chromosome(genes=newCollection(Gene(shape1=1, shape2=100),5))
   cr
   # list of two new identical Chromosome objects (different id)
   newCollection(cr, 2)                
   ni <- Niche(chromosomes = newCollection(cr, 2))
   ni # same genes values, different objects
   generateRandom(ni)
   ni # different genes values
 
   # creation and random generation at the same time
   ni <- Niche(chromosomes = newRandomCollection(cr, 2))
   ni
 }

 \references{Goldberg, David E. 1989 \emph{Genetic Algorithms in Search, Optimization and Machine Learning}. Addison-Wesley Pub. Co. ISBN: 0201157675}

 \author{Victor Trevino. Francesco Falciani Group. University of Birmingham, U.K. http://www.bip.bham.ac.uk/bioinf}

 \seealso{
   For more information see \code{\link{Chromosome}}.
   \code{\link{unObject}},
   \code{\link{as.list.Object}},
   \code{\link{newCollection.Chromosome}},
   \code{\link{newRandomCollection.Chromosome}},
   \code{\link{Niche}}.
 }

 \keyword{methods}
