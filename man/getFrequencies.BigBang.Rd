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

 \name{getFrequencies.BigBang}
\alias{getFrequencies.BigBang}
\alias{BigBang.getFrequencies}
\alias{getFrequencies.BigBang}
\alias{getFrequencies,BigBang-method}
\alias{getFrequencies}
\keyword{methods}
\keyword{internal}


 \title{Computes gene freqencies}

 \description{
  Computes gene freqencies.
 }

 \usage{\method{getFrequencies}{BigBang}(o, filter="none", subset=TRUE, ...)}

 \arguments{
	\item{filter}{The \code{BigBang} object can save information about solutions that did not reach the \code{goalFitness}. \code{filter=="solutions"} ensures that only chromosomes that reach the \code{goalFitness} are considered. \code{fitlter=="none"} take all chromosomes. \code{filter=="nosolutions"} consider only no-solutions (for comparative purposes).}
	\item{subset}{Second level of filter. \code{subset} can be a vector specifying which filtered chromosomes are used. It can be a logical vector or a numeric vector (indexes in order given by \code{$bestChromosomes} in \code{BigBang} object variable).}
 }

 \value{
  Returns a list of components.
  \item{new}{Current gene count.}
  \item{n}{Number of genes in new.}
  \item{nChr}{Number of chromosomes counted.}
  \item{ord}{Decreasing order.}
  \item{rnk}{Rank (starting with maximum).}
 }


 \references{Goldberg, David E. 1989 \emph{Genetic Algorithms in Search, Optimization and Machine Learning}. Addison-Wesley Pub. Co. ISBN: 0201157675}

 \author{Victor Trevino. Francesco Falciani Group. University of Birmingham, U.K. http://www.bip.bham.ac.uk/bioinf}

 \seealso{
   For more information see \code{\link{BigBang}}.
   \code{\link[galgo:geneFrequency.BigBang]{*geneFrequency}()}.
 }

 \keyword{methods}
