%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  
% 
% on Wed Feb 08 14:38:07 2006.
% 
% Generator was the Rdoc class, which is part of the R.oo package written
% by Henrik Bengtsson, 2001-2004.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{meanFitness.BigBang}
\alias{meanFitness.BigBang}
\alias{BigBang.meanFitness}
\alias{meanFitness.BigBang}
\alias{meanFitness,BigBang-method}
\alias{meanFitness}
\keyword{methods}
\keyword{internal}


 \title{Computes the ``mean'' fitness from several solutions}

 \description{
  Computes the ``mean'' fitness from several solutions.
 }

 \usage{\method{meanFitness}{BigBang}(o, filter="none", subset=TRUE, ...)}

 \arguments{
	\item{filter}{The \code{BigBang} object can save information about solutions that did not reach the \code{goalFitness}. \code{filter=="solutions"} ensures that only chromosomes that reach the \code{goalFitness} are considered. \code{fitlter=="none"} take all chromosomes. \code{filter=="nosolutions"} consider only no-solutions (for comparative purposes).}
	\item{subset}{Second level of filter. \code{subset} can be a vector specifying which filtered chromosomes are used. It can be a logical vector or a numeric vector (indexes in order given by \code{$bestChromosomes} in \code{BigBang} object variable). If it is a numeric vector length one, a positive value means take those top chromosomes sorted by fitness, a negative value take those at bottom.}
 }

 \details{
    The mean is built considering all solutions. For solutions that have finished earlier, the final fitness is used for futher genertions.
 }

 \value{
  A vector with the mean fitness in each generation.
 }

 \examples{\dontrun{
   #bb is a BigBang object
   geneRankStability(bb)
   }
 }

 \references{Goldberg, David E. 1989 \emph{Genetic Algorithms in Search, Optimization and Machine Learning}. Addison-Wesley Pub. Co. ISBN: 0201157675}

 \author{Victor Trevino. Francesco Falciani Group. University of Birmingham, U.K. http://www.bip.bham.ac.uk/bioinf}

 \seealso{
   For more information see \code{\link{BigBang}}.
 }

 \keyword{methods}
