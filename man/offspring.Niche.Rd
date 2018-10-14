%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  
% 
% on Wed Feb 08 14:37:46 2006.
% 
% Generator was the Rdoc class, which is part of the R.oo package written
% by Henrik Bengtsson, 2001-2004.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{offspring.Niche}
\alias{offspring.Niche}
\alias{Niche.offspring}
\alias{offspring.Niche}
\alias{offspring,Niche-method}
\alias{offspring}
\keyword{methods}
\keyword{internal}


 \title{Overwrites the new niche selecting a new population from the best chromosomes}

 \description{
  Overwrites the new niche selecting from the best chromosomes. This method is called in the \code{progeny} method. For more information see \code{\link{Niche}}. for complete description.
 }

 \usage{\method{offspring}{Niche}(ni, ...)}

 \details{
  The basic idea to generate a progeny is a selection biased toward
  the best chromosomes (see Goldberg). We implented this idea as a weighted
  probability (\code{\link[galgo:scaling.Niche]{*scaling}()}).

  \code{offspring} is part of \code{progeny} method.

  For related details For more information see \code{\link{Niche}}.
 }

 \value{
  Returns the selected chromosomes.
 }

 \examples{
   cr <- Chromosome(genes=newCollection(Gene(shape1=1, shape2=100),5))
   cr
   ni <- Niche(chromosomes = newRandomCollection(cr, 10))
   ni
   ni$fitness <- 1:10/10 # tricky fitness, only for showing purposes
   scaling(ni)
   offspring(ni)
   ni
 }

 \references{Goldberg, David E. 1989 \emph{Genetic Algorithms in Search, Optimization and Machine Learning}. Addison-Wesley Pub. Co. ISBN: 0201157675}

 \author{Victor Trevino. Francesco Falciani Group. University of Birmingham, U.K. http://www.bip.bham.ac.uk/bioinf}

 \seealso{
   For more information see \code{\link{Niche}}.
   \code{\link[galgo:progeny.Niche]{*progeny}()},
   \code{\link[galgo:scaling.Niche]{*scaling}()}.
 }

 \keyword{methods}
