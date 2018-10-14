%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  
% 
% on Wed Feb 08 14:38:09 2006.
% 
% Generator was the Rdoc class, which is part of the R.oo package written
% by Henrik Bengtsson, 2001-2004.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{geneCoverage.BigBang}
\alias{geneCoverage.BigBang}
\alias{BigBang.geneCoverage}
\alias{geneCoverage.BigBang}
\alias{geneCoverage,BigBang-method}
\alias{geneCoverage}
\keyword{methods}
\keyword{internal}


 \title{Computes the fraction of genes present in the top-rank from the total genes present in chromosomes}

 \description{
  Computes the fraction of genes present in the top-rank from the total genes present in chromosomes.
 }

 \usage{\method{geneCoverage}{BigBang}(o, filter="none", subset=TRUE, chromosomes=NULL, ...)}

 \arguments{
	\item{filter}{The \code{BigBang} object can save information about solutions that did not reach the \code{goalFitness}. \code{filter=="solutions"} ensures that only chromosomes that reach the \code{goalFitness} are considered. \code{fitlter=="none"} take all chromosomes. \code{filter=="nosolutions"} consider only no-solutions (for comparative purposes).}
	\item{subset}{Second level of filter. \code{subset} can be a vector specifying which filtered chromosomes are used. It can be a logical vector or a numeric vector (indexes in order given by \code{$bestChromosomes} in \code{BigBang} object variable). If it is a numeric vector length one, a positive value means take those top chromosomes sorted by fitness, a negative value take those at bottom.}
	\item{chromosomes}{The chromosomes to process. The default is using \code{filter} and \code{subset} to extract the chromosomes from the \code{BigBang} object.}
 }

 \value{
  A vector with the fraction of genes present in each rank from the total genes present in chromosomes.
 }

 \examples{\dontrun{
   #bb is a BigBang object
   gc <- geneCoverage(bb)
   gc
   }
 }

 \references{Goldberg, David E. 1989 \emph{Genetic Algorithms in Search, Optimization and Machine Learning}. Addison-Wesley Pub. Co. ISBN: 0201157675}

 \author{Victor Trevino. Francesco Falciani Group. University of Birmingham, U.K. http://www.bip.bham.ac.uk/bioinf}

 \seealso{
   For more information see \code{\link{BigBang}}.
   \code{\link[galgo:plot.BigBang]{*plot}()}.
 }

 \keyword{methods}
