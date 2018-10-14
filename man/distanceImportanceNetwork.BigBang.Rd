%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  
% 
% on Wed Feb 08 14:38:13 2006.
% 
% Generator was the Rdoc class, which is part of the R.oo package written
% by Henrik Bengtsson, 2001-2004.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{distanceImportanceNetwork.BigBang}
\alias{distanceImportanceNetwork.BigBang}
\alias{BigBang.distanceImportanceNetwork}
\alias{distanceImportanceNetwork.BigBang}
\alias{distanceImportanceNetwork,BigBang-method}
\alias{distanceImportanceNetwork}
\keyword{methods}
\keyword{internal}


 \title{Converts geneImportanceNetwork matrix to distance matrix}

 \description{
  Converts geneImportanceNetwork matrix to distance matrix.
 }

 \usage{\method{distanceImportanceNetwork}{BigBang}(o, gin=NULL, ...)}

 \arguments{
	\item{gin}{The geneImportanceNetwork matrix. If it is not provided, a default call to geneImportanceNetwork will be performed.}
	\item{...}{Parameters for geneImportanceNetwork when gin is not provided.}
 }

 \value{
	Returns a matrix representing gene distances.
 }

 \examples{\dontrun{
   #bb is a BigBang object
   bb
   gin <- geneImportanceNetwork(bb)
   din <- distanceImportanceNetwork(bb, gin)
   }
 }

 \references{Goldberg, David E. 1989 \emph{Genetic Algorithms in Search, Optimization and Machine Learning}. Addison-Wesley Pub. Co. ISBN: 0201157675}

 \author{Victor Trevino. Francesco Falciani Group. University of Birmingham, U.K. http://www.bip.bham.ac.uk/bioinf}

 \seealso{
   For more information see \code{\link{BigBang}}.
   \code{\link[galgo:geneImportanceNetwork.BigBang]{*geneImportanceNetwork}()}.
 }

 \keyword{methods}
