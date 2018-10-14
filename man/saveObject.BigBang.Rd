%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  
% 
% on Wed Feb 08 14:38:01 2006.
% 
% Generator was the Rdoc class, which is part of the R.oo package written
% by Henrik Bengtsson, 2001-2004.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{saveObject.BigBang}
\alias{saveObject.BigBang}
\alias{BigBang.saveObject}
\alias{saveObject.BigBang}
\alias{saveObject,BigBang-method}
\alias{saveObject}
\keyword{methods}
\keyword{internal}


 \title{Saves the BigBang object into a file in a suitable format}

 \description{
  Saves the BigBang object into a file in a suitable format.
 }

 \usage{\method{saveObject}{BigBang}(.bbO, file=.bbO$saveFile, mode=.bbO$saveMode, ...)}

 \arguments{
   \item{file}{The file name where the data will be saved. The default is taking the \code{$saveFile} variable form the \code{BigBang} object.}
   \item{saveMode}{Character vector specifying the saving mode. The default is taking the \code{$saveMode} variable from the \code{BigBang} object. Any combinations of the two options \code{compress} and \code{unObject}. It can be character vector length 1 or larger. For example, \code{saveMode=="compress+unObject"} would call \code{unObject} and save the file using \code{compress=TRUE}. The vector \code{c("object","compress")} (or shorter \code{c("compress")}) would save the \code{BigBang} object and compressed. It is not recommended to save the crude object because the functions varibles are stuck to environments and R will try to save those environments together, the result can be a waste of disk space and saving time. We strongly recommend \code{saveMode="unObject+compress"}.}
 }

 \value{
  Returns nothing.
 }

 \examples{\dontrun{
   cr <- Chromosome(genes=newCollection(Gene(shape1=1, shape2=100),5))
   ni <- Niche(chromosomes=newRandomCollection(cr, 10))
   wo <- World(niches=newRandomCollection(ni,2))
   ga <- Galgo(populations=newRandomCollection(wo,1), goalFitness = 0.75,
               callBackFunc=plot,
               fitnessFunc=function(chr, parent) 5/sd(as.numeric(chr)))
 
   #evolve(ga) ## not needed here

   bb <- BigBang(galgo=ga, maxSolutions=10, maxBigBangs=10, saveGeneBreaks=1:100)
   blast(bb, 1)
   saveObject(bb, file="bb.Rdata", mode="unObject+compress")
   }
 }

 \references{Goldberg, David E. 1989 \emph{Genetic Algorithms in Search, Optimization and Machine Learning}. Addison-Wesley Pub. Co. ISBN: 0201157675}

 \author{Victor Trevino. Francesco Falciani Group. University of Birmingham, U.K. http://www.bip.bham.ac.uk/bioinf}

 \seealso{
   For more information see \code{\link{BigBang}}.
 }

 \keyword{methods}
