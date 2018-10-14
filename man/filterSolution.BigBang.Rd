%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  
% 
% on Wed Feb 08 14:38:03 2006.
% 
% Generator was the Rdoc class, which is part of the R.oo package written
% by Henrik Bengtsson, 2001-2004.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{filterSolution.BigBang}
\alias{filterSolution.BigBang}
\alias{BigBang.filterSolution}
\alias{filterSolution.BigBang}
\alias{filterSolution,BigBang-method}
\alias{filterSolution}
\keyword{methods}
\keyword{internal}


 \title{Filters solutions}

 \description{
  Filters solutions.
 }

 \usage{\method{filterSolution}{BigBang}(o, filter=c("none", "solutions", "nosolutions"), subset, ...)}

 \arguments{
	\item{filter}{The \code{BigBang} object can save information about solutions that did not reach the \code{goalFitness}. \code{filter=="solutions"} ensures that only chromosomes that reach the \code{goalFitness} are considered. \code{fitlter=="none"} take all chromosomes. \code{filter=="nosolutions"} consider only no-solutions (for comparative purposes).}
	\item{subset}{Second level of filter. \code{subset} can be a vector specifying which filtered chromosomes are used. It can be a logical vector or a numeric vector (indexes in order given by \code{$bestChromosomes} in \code{BigBang} object variable).}
 }

 \value{
  Returns a logical vector indicating which chromosomes are valid for the given filter.
 }

 \examples{\dontrun{
   ok <- filterSolution(bb, filter="solutions", subset=1:100)
   sum(ok)
   }
 }

 \references{Goldberg, David E. 1989 \emph{Genetic Algorithms in Search, Optimization and Machine Learning}. Addison-Wesley Pub. Co. ISBN: 0201157675}

 \author{Victor Trevino. Francesco Falciani Group. University of Birmingham, U.K. http://www.bip.bham.ac.uk/bioinf}

 \seealso{
   For more information see \code{\link{BigBang}}.
 }

 \keyword{methods}
