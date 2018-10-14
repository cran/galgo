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

\name{formatChromosome.BigBang}
\alias{formatChromosome.BigBang}
\alias{BigBang.formatChromosome}
\alias{formatChromosome.BigBang}
\alias{formatChromosome,BigBang-method}
\alias{formatChromosome}
\keyword{methods}
\keyword{internal}


\title{Converts chromosome for storage in BigBang object}

\description{
 Converts chromosome for storage in BigBang object. The default behaviour is convert the chromosome in a numeric format instead of the original Chromosome object to save memory and time during further calls to saveObject and loadObject methods.
}

\usage{\method{formatChromosome}{BigBang}(.bbO, chr, ...)}

\arguments{
  \item{chr}{Chromosome to format}
}

\value{
 Object to be stored in $bestChromosomes variable in BigBang object.
}

\references{Goldberg, David E. 1989 \emph{Genetic Algorithms in Search, Optimization and Machine Learning}. Addison-Wesley Pub. Co. ISBN: 0201157675}

\author{Victor Trevino. Francesco Falciani Group. University of Birmingham, U.K. http://www.bip.bham.ac.uk/bioinf}

\seealso{
  For more information see \code{\link{BigBang}}.
  \code{\link{saveObject.BigBang}},
  \code{\link{loadObject}}.
}

\keyword{methods}
