%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  
% 
% on Wed Feb 08 14:37:35 2006.
% 
% Generator was the Rdoc class, which is part of the R.oo package written
% by Henrik Bengtsson, 2001-2004.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{as.double.Chromosome}
\alias{as.double.Chromosome}
\alias{Chromosome.as.double}
\alias{as.double.Chromosome}
\alias{as.double,Chromosome-method}
\keyword{methods}
\keyword{internal}


 \title{Converts the chromosome values (genes) to its numerical representation}

 \description{
  Converts the chromosome values (genes) to its numerical representation.
 }

 \usage{\method{as.double}{Chromosome}(x, ...)}

 \details{
	This function really calls \code{genes} method.
 }

 \value{
  Returns a vector containig the values for all genes in the chromosome.
 }

 \examples{
   cr <- Chromosome(genes=newCollection(Gene(shape1=1, shape2=100),5))
   cr
   as.double(cr) # the same
   as.numeric(cr) # the same
   as.vector(cr) # NA's is not the same
 }

 \references{Goldberg, David E. 1989 \emph{Genetic Algorithms in Search, Optimization and Machine Learning}. Addison-Wesley Pub. Co. ISBN: 0201157675}

 \author{Victor Trevino. Francesco Falciani Group. University of Birmingham, U.K. http://www.bip.bham.ac.uk/bioinf}

 \seealso{
   For more information see \code{\link{Chromosome}}.
   \code{\link[galgo:genes.Chromosome]{*genes}()}.
 }

 \keyword{methods}
