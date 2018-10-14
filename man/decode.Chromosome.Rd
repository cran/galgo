%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  
% 
% on Wed Feb 08 14:37:36 2006.
% 
% Generator was the Rdoc class, which is part of the R.oo package written
% by Henrik Bengtsson, 2001-2004.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{decode.Chromosome}
\alias{decode.Chromosome}
\alias{Chromosome.decode}
\alias{decode.Chromosome}
\alias{decode,Chromosome-method}
\alias{decode}
\keyword{methods}
\keyword{internal}


 \title{Converts the gene values to user-readable values}

 \description{
  Converts the gene values to user-readable values.
 }

 \usage{\method{decode}{Chromosome}(x, ...)}

 \details{
	This function really calls the function defined in the \code{$decode} variable in the \code{Chromosome} object.
 }

 \value{
  Returns the representation of the chromosome.
 }

 \examples{
   cr <- Chromosome(genes=newCollection(Gene(shape1=1, shape2=100),5))
   cr
   as.double(cr) # the same
   as.numeric(cr) # the same
   decode(cr) # the same
   as.vector(cr) # NA's is not the same
 }

 \references{Goldberg, David E. 1989 \emph{Genetic Algorithms in Search, Optimization and Machine Learning}. Addison-Wesley Pub. Co. ISBN: 0201157675}

 \author{Victor Trevino. Francesco Falciani Group. University of Birmingham, U.K. http://www.bip.bham.ac.uk/bioinf}

 \seealso{
   For more information see \code{\link{Chromosome}}.
   \code{\link[galgo:genes.Chromosome]{*genes}()}.
 }

 \keyword{methods}
