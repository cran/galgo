%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  
% 
% on Wed Feb 08 14:37:30 2006.
% 
% Generator was the Rdoc class, which is part of the R.oo package written
% by Henrik Bengtsson, 2001-2004.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{as.matrix.Gene}
\alias{as.matrix.Gene}
\alias{Gene.as.matrix}
\alias{as.matrix.Gene}
\alias{as.matrix,Gene-method}
\keyword{methods}
\keyword{internal}


\title{Converts the gene parameters (shape1, shape2) to matrix}

\description{
 Converts the gene parameters (shape1, shape2) to matrix. It is suitable to bind rows for many genes (like in a chromosome).
}

\usage{\method{as.matrix}{Gene}(x, ...)}

\value{
 Returns a one-row matrix containig id, shape1, and shape2.
}

\examples{
  ge <- Gene(shape1=1, shape2=100)
  ge
  as.matrix(ge)
}

\references{Goldberg, David E. 1989 \emph{Genetic Algorithms in Search, Optimization and Machine Learning}. Addison-Wesley Pub. Co. ISBN: 0201157675}

\author{Victor Trevino. Francesco Falciani Group. University of Birmingham, U.K. http://www.bip.bham.ac.uk/bioinf}

\seealso{
  For more information see \code{\link{Gene}}, \code{\link{Chromosome}}
}

\keyword{methods}
