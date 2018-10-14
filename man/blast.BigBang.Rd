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

 \name{blast.BigBang}
\alias{blast.BigBang}
\alias{BigBang.blast}
\alias{blast.BigBang}
\alias{blast,BigBang-method}
\alias{blast}
\keyword{methods}
\keyword{internal}


 \title{Evolves Galgo objects saving the results for further analysis}

 \description{
  The basic process is as follows.\\n
  \\tab1. Clone \code{Galgo} and generate random chromosomes\\n
  \\tab2. Call \code{evolve} method\\n
  \\tab3. Save results in \code{BigBang} object\\n
  \\tab4. Verify stop rules\\n
  \\tab5. Goto 1\\n
 }

 \usage{\method{blast}{BigBang}(.bb, add=0, ...)}
 \arguments{
   \item{add}{Force to add a number to maxBigBangs and maxSolutions in order to search for more solutions.}
 }

 \value{
  Returns nothing. The results are saved in the the \code{BigBang} object.
 }

 \examples{
   cr <- Chromosome(genes=newCollection(Gene(shape1=1, shape2=100),5))
   ni <- Niche(chromosomes=newRandomCollection(cr, 10))
   wo <- World(niches=newRandomCollection(ni,2))
   ga <- Galgo(populations=newRandomCollection(wo,1), goalFitness = 0.75,
				callBackFunc=plot,
               fitnessFunc=function(chr, parent) 5/sd(as.numeric(chr)))
 
   #evolve(ga) ## not needed here

   bb <- BigBang(galgo=ga, maxSolutions=10, maxBigBangs=10, saveGeneBreaks=1:100)
   \dontrun{blast(bb)}
   \dontrun{plot(bb)}
   \dontrun{blast(bb, 1)}
   \dontrun{plot(bb)}
 }

 \references{Goldberg, David E. 1989 \emph{Genetic Algorithms in Search, Optimization and Machine Learning}. Addison-Wesley Pub. Co. ISBN: 0201157675}

 \author{Victor Trevino. Francesco Falciani Group. University of Birmingham, U.K. http://www.bip.bham.ac.uk/bioinf}

 \seealso{
   For more information see \code{\link{BigBang}}.
   \code{\link{evolve.Galgo}}().
 }

 \keyword{methods}
