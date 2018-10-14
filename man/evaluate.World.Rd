%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  
% 
% on Wed Feb 08 14:37:54 2006.
% 
% Generator was the Rdoc class, which is part of the R.oo package written
% by Henrik Bengtsson, 2001-2004.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{evaluate.World}
\alias{evaluate.World}
\alias{World.evaluate}
\alias{evaluate.World}
\alias{evaluate,World-method}
\keyword{methods}
\keyword{internal}


 \title{Evaluate all niches with a fitness function}

 \description{
  Evaluate all niches with a fitness function. The result of this evaluation is treated as the ``fitness'' value as defined by Goldberg (see references).
  The \code{Galgo} object call this method and store the resulted value in order to decide which chromosomes are better choices to be part of the next generation.
  The ``fitness function'' should returns a numeric value scaled from 0 to 1. As close to 1 as better chance it have to be part of the next generation.
 }

 \usage{\method{evaluate}{World}(.O, fn, parent, ...)}

 \arguments{
   \item{fn}{The ``fitness'' function to be called to evaluate all niches. It should follow the format \code{function(obj, parent) \{ ... \}}}
   \item{parent}{The original object calling for the evaluation. This is passed when the function is sensitive to data stored in parent object. Commonly it is a \code{BigBang} object (perhaps \code{Galgo} instead).}
 }

 \value{
  Returns nothing.
 }

 \examples{
   cr <- Chromosome(genes=newCollection(Gene(shape1=1, shape2=100),5))
   cr
   ni <- Niche(chromosomes = newRandomCollection(cr, 10))
   ni
   fn <- function(chr, parent) { sd(as.double(chr))/mean(as.double(chr)) }
   evaluate(ni, fn, parent)
   getFitness(ni) ## see results
   summary(ni)
   wo <- World(niches=newRandomCollection(ni,2))
   evaluate(wo, fn, parent)
 }

 \references{Goldberg, David E. 1989 \emph{Genetic Algorithms in Search, Optimization and Machine Learning}. Addison-Wesley Pub. Co. ISBN: 0201157675}

 \author{Victor Trevino. Francesco Falciani Group. University of Birmingham, U.K. http://www.bip.bham.ac.uk/bioinf}

 \seealso{
   For more information see \code{\link{World}}.
 }

 \keyword{methods}
