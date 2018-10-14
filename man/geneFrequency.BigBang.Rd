%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  
% 
% on Wed Feb 08 14:38:06 2006.
% 
% Generator was the Rdoc class, which is part of the R.oo package written
% by Henrik Bengtsson, 2001-2004.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{geneFrequency.BigBang}
\alias{geneFrequency.BigBang}
\alias{BigBang.geneFrequency}
\alias{geneFrequency.BigBang}
\alias{geneFrequency,BigBang-method}
\alias{geneFrequency}
\keyword{methods}
\keyword{internal}


 \title{Computes the frequency of genes based on chromosomes}

 \description{
  Computes the frequency of genes based on chromosomes. It really returns \code{getFrequiences} using the \code{$new} variable, adding gene names, and filtering for \code{cutoff}.
 }

 \usage{\method{geneFrequency}{BigBang}(o,
	filter="none",
	subset=TRUE,
	gene.names=TRUE,
	cutoff=-1,
	value=c("frequency", "indexes", "ranks"),
	...)}

 \arguments{
	\item{filter}{The \code{BigBang} object can save information about solutions that did not reach the \code{goalFitness}. \code{filter=="solutions"} ensures that only chromosomes that reach the \code{goalFitness} are considered. \code{fitlter=="none"} take all chromosomes. \code{filter=="nosolutions"} consider only no-solutions (for comparative purposes).}
	\item{subset}{Second level of filter. \code{subset} can be a vector specifying which filtered chromosomes are used. It can be a logical vector or a numeric vector (indexes in order given by \code{$bestChromosomes} in \code{BigBang} object variable). If it is a numeric vector length one, a positive value means take those top chromosomes sorted by fitness, a negative value take those at bottom.}
	\item{gene.names}{\code{TRUE} for naming the result with the stored \code{$geneNames} in oject \code{BigBang}. Other character vector to name-specific.}
	\item{cutoff}{Only genes whose frequency is greather than \code{cutoff} are repored.}
	\item{value}{The result. \code{"frequency","indexes","ranks"}}
 }

 \value{
  A table when \code{value=="frequency"}, otherwise, a vector.
 }

 \examples{\dontrun{
   #bb is a BigBang object
   geneFrequency(bb)
   }
 }

 \references{Goldberg, David E. 1989 \emph{Genetic Algorithms in Search, Optimization and Machine Learning}. Addison-Wesley Pub. Co. ISBN: 0201157675}

 \author{Victor Trevino. Francesco Falciani Group. University of Birmingham, U.K. http://www.bip.bham.ac.uk/bioinf}

 \seealso{
   For more information see \code{\link{BigBang}}.
 }

 \keyword{methods}
