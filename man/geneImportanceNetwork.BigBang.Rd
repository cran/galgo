%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  
% 
% on Wed Feb 08 14:38:12 2006.
% 
% Generator was the Rdoc class, which is part of the R.oo package written
% by Henrik Bengtsson, 2001-2004.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{geneImportanceNetwork.BigBang}
\alias{geneImportanceNetwork.BigBang}
\alias{BigBang.geneImportanceNetwork}
\alias{geneImportanceNetwork.BigBang}
\alias{geneImportanceNetwork,BigBang-method}
\alias{geneImportanceNetwork}
\keyword{methods}
\keyword{internal}


 \title{Computes the number of times a couple of top-ranked-genes are present in models}

 \description{
  Computes the number of times top-ranked-genes are present in models.
 }

 \usage{\method{geneImportanceNetwork}{BigBang}(o, 
 filter="none", 
 subset=TRUE, 
 mord=min(ncol(o$data$data),50), 
 inc.rank=FALSE, 
 inc.index=FALSE, ...)}

 \arguments{
	\item{filter}{The \code{BigBang} object can save information about solutions that did not reach the \code{goalFitness}. \code{filter=="solutions"} ensures that only chromosomes that reach the \code{goalFitness} are considered. \code{fitlter=="none"} take all chromosomes. \code{filter=="nosolutions"} consider only no-solutions (for comparative purposes).}
	\item{subset}{Second level of filter. \code{subset} can be a vector specifying which filtered chromosomes are used. It can be a logical vector or a numeric vector (indexes in order given by \code{$bestChromosomes} in \code{BigBang} object variable). If it is a numeric vector length one, a positive value means take those top chromosomes sorted by fitness, a negative value take those at bottom.}
	\item{mord}{The number of ``top-ranked-genes'' to highlight.}
	\item{inc.rank}{Incluye the gene rank in rownames and colnames.}
	\item{inc.index}{Incluye the gene index in rownames and colnames.}
 }

 \value{
	Returns a matrix with number of overlaps for every top-ranked-gene pairs. The order correspond to rank.
 }

 \examples{\dontrun{
   #bb is a BigBang object
   bb
   gin <- geneImportanceNetwork(bb)
   gin
   gin <- geneImportanceNetwork(bbm, mord=5)
   gin
   }
 }

 \references{Goldberg, David E. 1989 \emph{Genetic Algorithms in Search, Optimization and Machine Learning}. Addison-Wesley Pub. Co. ISBN: 0201157675}

 \author{Victor Trevino. Francesco Falciani Group. University of Birmingham, U.K. http://www.bip.bham.ac.uk/bioinf}

 \seealso{
   For more information see \code{\link{BigBang}}.
   \code{\link[galgo:distanceImportanceNetwork.BigBang]{*distanceImportanceNetwork}()}.
 }

 \keyword{methods}
