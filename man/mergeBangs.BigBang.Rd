%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  
% 
% on Wed Feb 08 14:38:11 2006.
% 
% Generator was the Rdoc class, which is part of the R.oo package written
% by Henrik Bengtsson, 2001-2004.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{mergeBangs.BigBang}
\alias{mergeBangs.BigBang}
\alias{BigBang.mergeBangs}
\alias{mergeBangs.BigBang}
\alias{mergeBangs,BigBang-method}
\alias{mergeBangs}
\keyword{methods}
\keyword{internal}


 \title{Merges the information from other BigBang objects}

 \description{
  Merges the information from other BigBang objects. It is assumed that all \code{BigBang} objects compute the same results, thus they are ``mergeable''. Possibly from parallelization.
 }

 \usage{\method{mergeBangs}{BigBang}(mb,  ..., from=NULL, to=NULL, detect.duplicates=TRUE, clone=FALSE)}

 \arguments{
	\item{...}{List of \code{BigBang} objects.}
	\item{clone}{Logical. Specify if the original object must be cloned before merging.}
 }

 \value{
	Returns nothing. However, the original \code{BigBang} object has been modified adding the others \code{BigBang} objects.
 }

 \examples{\dontrun{
   #bb, bbmachine2, bbmachine3 are BigBang objects
   bb
   plot(bb)
   mergeBangs(bb, bbmachine2, bbmachine3)
   bb       # accumulated solutions
   plot(bb) # accumulated solutions
   }
 }

 \references{Goldberg, David E. 1989 \emph{Genetic Algorithms in Search, Optimization and Machine Learning}. Addison-Wesley Pub. Co. ISBN: 0201157675}

 \author{Victor Trevino. Francesco Falciani Group. University of Birmingham, U.K. http://www.bip.bham.ac.uk/bioinf}

 \seealso{
   For more information see \code{\link{BigBang}}.
 }

 \keyword{methods}
