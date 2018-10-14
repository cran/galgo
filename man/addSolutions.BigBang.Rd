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

 \name{addSolutions.BigBang}
\alias{addSolutions.BigBang}
\alias{BigBang.addSolutions}
\alias{addSolutions.BigBang}
\alias{addSolutions,BigBang-method}
\alias{addSolutions}
\keyword{methods}
\keyword{internal}


 \title{Adds user-built or external chromosomes as solutions}

 \description{
  Adds chromosomes as solutions. This method is used to build BigBang objects from an external source of chromosomes/solutions. These solutions were obtained from any other process. This can be used for simulations, verifications and comparisons. 
 }

 \usage{\method{addSolutions}{BigBang}(no, chrs, gen=NULL, sol=NULL, fit=NULL, ...)}

 \arguments{
  \item{no}{BigBang Object}
	\item{chrs}{Chromosomes to add. It is assumed to be a list. If data frame or matrix are provided, rows are considered as chromosomes and columns as genes/variables.}
	\item{gen}{Generation in which each solution was reached (numeric vector). Default to NULL, which is interpreted as 1.}
	\item{sol}{Vector of flags indicating whether each chromosome is a solution (1=TRUE, 0=FALSE). Default to NULL, which is interpreted as 1.}
	\item{fit}{Vector of fitness of each chromosome (between 0 and 1). Default to NULL, which is used to compute the fitness using internal configuration of the BigBang/Galgo objects.}
 }

 \value{
	Returns nothing. However, the original \code{BigBang} object has been modified adding the solutions.
 }

 \examples{
   #bb is a BigBang object
   \dontrun{
   bb
   plot(bb)
   addSolutions(bb, chrList)
   plot(bb) # accumulated solutions
   }
 }

 \references{Goldberg, David E. 1989 \emph{Genetic Algorithms in Search, Optimization and Machine Learning}. Addison-Wesley Pub. Co. ISBN: 0201157675}

 \author{Victor Trevino. Francesco Falciani Group. University of Birmingham, U.K. http://www.bip.bham.ac.uk/bioinf}

 \seealso{
   For more information see \code{\link{BigBang}}.
 }

 \keyword{methods}
