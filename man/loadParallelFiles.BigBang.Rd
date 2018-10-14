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

 \name{loadParallelFiles.BigBang}
\alias{loadParallelFiles.BigBang}
\alias{BigBang.loadParallelFiles}
\alias{loadParallelFiles.BigBang}
\alias{loadParallelFiles,BigBang-method}
\alias{loadParallelFiles}
\keyword{methods}
\keyword{internal}


 \title{Load all files saved during the parallelization}

 \description{
  Load all files saved during the parallelization.
 }

 \usage{\method{loadParallelFiles}{BigBang}(bb, prefix="parallel-", ...)}

 \arguments{
	\item{prefix}{Prefix used in the parallel files for identification. Defaults to ``parallel-''.}
 }

 \value{
	Returns the file names loaded to the bigbang object.
 }

 \examples{\dontrun{

	#Initial process:
	#load data and configure initial objects, run once
	library(galgo)
	bb <- configBB.varSel(..., saveFile="bb.parallel.Rdata", ...)
	saveObject(bb)
	#

	#Parallel process:
	#run as many process as you want
	library(galgo)
	loadObject("bb.parallel.Rdata")
	assignParallelFile(bb)
	blast(bb)
	#


	#Analysis Process:
	library(galgo)
	loadObject("bb.parallel.Rdata")
	loadParallelFiles(bb)
	#
	}
 }

 \references{Goldberg, David E. 1989 \emph{Genetic Algorithms in Search, Optimization and Machine Learning}. Addison-Wesley Pub. Co. ISBN: 0201157675}

 \author{Victor Trevino. Francesco Falciani Group. University of Birmingham, U.K. http://www.bip.bham.ac.uk/bioinf}

 \seealso{
   For more information see \code{\link{BigBang}}.
 }

 \keyword{methods}
