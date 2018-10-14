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

 \name{addRandomSolutions.BigBang}
\alias{addRandomSolutions.BigBang}
\alias{BigBang.addRandomSolutions}
\alias{addRandomSolutions.BigBang}
\alias{addRandomSolutions,BigBang-method}
\alias{addRandomSolutions}
\keyword{methods}
\keyword{internal}


 \title{Adds random pre-existed solutions}

 \description{
  Adds random pre-existed solutions. The purpose is to ``simulate'' new solutions in order to see gene stability and gene frequency. The results are biased but can be used to estimate the how many more solutions are needed to stabilize certain number of genes.
 }

 \usage{\method{addRandomSolutions}{BigBang}(no, n=length(no$bestChromosomes), ...)}

 \arguments{
	\item{n}{number of solutions to add}
 }

 \value{
	Returns nothing. However, the original \code{BigBang} object has been modified adding the random solutions.
 }

 \examples{
   #bb is a BigBang object
   \dontrun{
   bb
   plot(bb)
   addRandomSolutions(bb, 1000)
   plot(bb) # accumulated solutions
   }
 }

 \references{Goldberg, David E. 1989 \emph{Genetic Algorithms in Search, Optimization and Machine Learning}. Addison-Wesley Pub. Co. ISBN: 0201157675}

 \author{Victor Trevino. Francesco Falciani Group. University of Birmingham, U.K. http://www.bip.bham.ac.uk/bioinf}

 \seealso{
   For more information see \code{\link{BigBang}}.
 }

 \keyword{methods}
