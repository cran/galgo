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

 \name{assignParallelFile.BigBang}
\alias{assignParallelFile.BigBang}
\alias{BigBang.assignParallelFile}
\alias{assignParallelFile.BigBang}
\alias{assignParallelFile,BigBang-method}
\alias{assignParallelFile}
\keyword{methods}
\keyword{internal}


 \title{Assigns a different saveFile value for parallelization}

 \description{
  Assigns a different saveFile value for parallelization. The assignation is done looking for a existing filename with a consecutive number starting from 1. If a filename with a particular consecutive number is not found, it is assigned.
 }

 \usage{\method{assignParallelFile}{BigBang}(bb,
	save=TRUE,
	prefix="parallel-",
	use.random=TRUE,
	compatibilize=TRUE,
	...)}

 \arguments{
	\item{save}{Specify to save the file to ``lock'' the name and avoid other parallel process to use it. Defaults to true.}
	\item{prefix}{Prefix used in the parallel files for identification. Defaults to ``parallel-''.}
	\item{use.random}{Specify if an integer random value between 0 and 9999 is added to avoid duplicated filenames in parallel process started at the same time. Defaults to true.}
 }

 \value{
	Returns the file name assigned to the bigbang object for parallelization.
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
